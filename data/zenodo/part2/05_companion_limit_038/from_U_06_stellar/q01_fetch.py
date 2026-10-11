"""Audit 6 (stellar parameters), step 1: fresh catalogue queries for TIC 4206066 = Gaia DR3 3220388198192519424.

Every file written to data/ is logged in data/PROVENANCE.jsonl with service URL, query text, UTC time and sha256.
Services:
  Gaia archive TAP  https://gea.esac.esa.int/tap-server/tap  (gaiadr3.gaia_source, gaiadr3.astrophysical_parameters,
                    external.gaiaedr3_distance = Bailer-Jones+2021, 100-pc comparison sample)
  MAST TIC 8.2      astroquery.mast.Catalogs (catalog="TIC")
  VizieR            https://vizier.cds.unistra.fr  (II/246 2MASS PSC, II/328 AllWISE, II/365 CatWISE2020,
                    I/259 Tycho-2, II/336 APASS DR9, I/352 Bailer-Jones 2021)
"""
import hashlib, json, sys, time, datetime as dt
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
D = HERE / "data"; D.mkdir(exist_ok=True)
PROV = D / "PROVENANCE.jsonl"
SID = 3220388198192519424
TIC = 4206066


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def log(fn, service, query, note=""):
    rec = dict(file=str(Path(fn).name), service=service, query=query, utc=dt.datetime.now(dt.timezone.utc).isoformat(),
               sha256=sha(fn), note=note)
    with open(PROV, "a") as f:
        f.write(json.dumps(rec) + "\n")
    print("saved", fn, rec["sha256"][:12])


def tap(query, fn, note="", asynch=False):
    import pyvo
    svc = pyvo.dal.TAPService("https://gea.esac.esa.int/tap-server/tap")
    for attempt in range(3):
        try:
            if not asynch:
                t = svc.run_sync(query, maxrec=200000).to_table()
            else:
                job = svc.submit_job(query); job.run(); job.wait(phases=["COMPLETED", "ERROR", "ABORTED"], timeout=900)
                job.raise_if_error(); t = job.fetch_result().to_table()
            break
        except Exception as e:  # noqa
            print("TAP retry", attempt, e); time.sleep(10)
    else:
        raise RuntimeError("TAP failed")
    t.write(D / fn, format="ascii.ecsv", overwrite=True)
    log(D / fn, "Gaia archive TAP https://gea.esac.esa.int/tap-server/tap", query, note)
    return t


what = sys.argv[1:] or ["gaia", "bj", "tic", "vizier", "sample", "xp"]

if "gaia" in what:
    tap(f"SELECT * FROM gaiadr3.gaia_source WHERE source_id = {SID}", "gaia_dr3_source.ecsv")
    tap(f"SELECT * FROM gaiadr3.astrophysical_parameters WHERE source_id = {SID}", "gaia_dr3_ap.ecsv")
if "cone" in what:
    # neighbours within 30 arcsec (for blending of 2MASS/WISE apertures)
    tap("SELECT source_id, ra, dec, parallax, phot_g_mean_mag, bp_rp, DISTANCE(POINT(ra, dec), "
        "POINT(80.859827, -1.326964))*3600 AS sep_arcsec FROM gaiadr3.gaia_source WHERE "
        "1 = CONTAINS(POINT(ra, dec), CIRCLE(80.859827, -1.326964, 60.0/3600.)) ORDER BY sep_arcsec",
        "gaia_dr3_neighbours_60as.ecsv")

if "bj" in what:
    tap(f"SELECT * FROM external.gaiaedr3_distance WHERE source_id = {SID}", "bailerjones2021.ecsv",
        "Bailer-Jones et al. 2021 AJ 161, 147 (external.gaiaedr3_distance)")

if "tic" in what:
    from astroquery.mast import Catalogs
    t = Catalogs.query_criteria(catalog="TIC", ID=TIC)
    t.write(D / "tic82.ecsv", format="ascii.ecsv", overwrite=True)
    log(D / "tic82.ecsv", "MAST Catalogs API (astroquery.mast.Catalogs, catalog=TIC)", f"ID={TIC}")

if "vizier" in what:
    from astroquery.vizier import Vizier
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(80.859827 * u.deg, -1.326964 * u.deg)
    v = Vizier(columns=["**", "+_r"], row_limit=-1)
    for cat, r, fn in [("II/246/out", 10, "2mass_psc.ecsv"), ("II/328/allwise", 10, "allwise.ecsv"),
                       ("II/365/catwise", 10, "catwise2020.ecsv"), ("I/259/tyc2", 10, "tycho2.ecsv"),
                       ("II/336/apass9", 10, "apass9.ecsv"), ("I/352/gedr3dis", 5, "bj2021_vizier.ecsv")]:
        res = v.query_region(c, radius=r * u.arcsec, catalog=cat)
        if len(res) == 0:
            print("no rows", cat); continue
        t = res[0]
        t.write(D / fn, format="ascii.ecsv", overwrite=True)
        log(D / fn, "VizieR https://vizier.cds.unistra.fr (astroquery.vizier)", f"{cat} cone {r} arcsec around 80.859827 -1.326964")

if "sample" in what:
    # local main-sequence comparison sample (d < 100 pc, precise astrometry, clean photometry), colour window
    # around the target BP-RP = 1.739 (+-0.10). Quality cuts follow Gaia Collaboration, Smart et al. 2021 (GCNS)
    # in spirit: parallax_over_error > 20, RUWE < 1.4, corrected BP/RP excess |C*| < 3 sigma_C* (Riello+2021 eq. 6,18).
    q = ("SELECT source_id, parallax, parallax_error, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, bp_rp, "
         "ruwe, phot_bp_rp_excess_factor, ipd_frac_multi_peak, ipd_gof_harmonic_amplitude, non_single_star, "
         "duplicated_source, radial_velocity, radial_velocity_error, rv_chisq_pvalue, rv_nb_transits, rv_renormalised_gof, "
         "has_xp_continuous FROM gaiadr3.gaia_source WHERE parallax > 10 AND parallax_over_error > 20 "
         "AND bp_rp BETWEEN 1.54 AND 1.94 AND phot_bp_mean_mag IS NOT NULL AND phot_rp_mean_mag IS NOT NULL")
    tap(q, "gaia_100pc_bprp_1.54_1.94.ecsv", "local comparison sample for CMD overluminosity", asynch=False)

if "xp" in what:
    from astroquery.gaia import Gaia
    dl = Gaia.load_data(ids=[SID], data_release="Gaia DR3", retrieval_type="XP_SAMPLED", data_structure="INDIVIDUAL",
                        format="csv")
    for k, v in dl.items():
        tab = v[0].to_table() if hasattr(v[0], "to_table") else v[0]
        fn = D / "gaia_dr3_xp_sampled.ecsv"
        tab.write(fn, format="ascii.ecsv", overwrite=True)
        log(fn, "Gaia archive DataLink (astroquery.gaia.Gaia.load_data, XP_SAMPLED, Gaia DR3)", f"source_id={SID}", k)
        break

if "irsample" in what:
    # 100-pc stars with colours like the target, with 2MASS and AllWISE photometry via the Gaia DR3 best-neighbour
    # cross-matches (Marrese et al. 2019/2022), for an empirical photospheric Ks-W1..W4 locus.
    q = ("SELECT g.source_id, g.parallax, g.parallax_error, g.phot_g_mean_mag, g.bp_rp, g.ruwe, g.non_single_star, "
         "t.j_m, t.j_msigcom, t.h_m, t.h_msigcom, t.ks_m, t.ks_msigcom, t.ph_qual AS tm_qual, "
         "w.w1mpro, w.w1mpro_error, w.w2mpro, w.w2mpro_error, w.w3mpro, w.w3mpro_error, w.w4mpro, w.w4mpro_error, "
         "w.cc_flags, w.ext_flag, w.ph_qual AS w_qual, xw.angular_distance AS w_sep "
         "FROM gaiadr3.gaia_source AS g "
         "JOIN gaiadr3.allwise_best_neighbour AS xw ON xw.source_id = g.source_id "
         "JOIN gaiadr1.allwise_original_valid AS w ON w.allwise_oid = xw.allwise_oid "
         "JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS xt ON xt.source_id = g.source_id "
         "JOIN gaiadr3.tmass_psc_xsc_join AS tj ON tj.clean_tmass_psc_xsc_oid = xt.clean_tmass_psc_xsc_oid "
         "JOIN gaiadr1.tmass_original_valid AS t ON t.designation = tj.original_psc_source_id "
         "WHERE g.parallax > 10 AND g.parallax_over_error > 20 AND g.bp_rp BETWEEN 1.64 AND 1.84")
    tap(q, "gaia_100pc_bprp_1.64_1.84_2mass_allwise.ecsv", "empirical photospheric IR colour locus", asynch=False)
