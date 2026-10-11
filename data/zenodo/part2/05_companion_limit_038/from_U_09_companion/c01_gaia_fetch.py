"""Audit 9, step 1-2: fresh Gaia DR3 (gaia_source, NSS tables, neighbours <10"), RUWE comparison sample,
2MASS/AllWISE neighbours. All queries logged with URL, query text, UTC time, sha256 of saved file."""
import hashlib, json, datetime, sys
from pathlib import Path
import numpy as np
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord

HERE = Path(__file__).resolve().parent
D = HERE / "data"; D.mkdir(exist_ok=True)
PROV = D / "PROVENANCE.jsonl"
SID = 3220388198192519424
RA, DEC = 80.85983, -1.32696
Gaia.ROW_LIMIT = -1
import os
SYNC = os.environ.get("SYNC", "0") == "1"

def log(name, url, query, path):
    h = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    rec = dict(name=name, url=url, query=query, file=str(Path(path).name), sha256=h,
               utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"))
    with open(PROV, "a") as f: f.write(json.dumps(rec) + "\n")
    print(name, h[:16], rec["utc"])

def tap(name, q):
    job = Gaia.launch_job(q) if SYNC else Gaia.launch_job_async(q)
    t = job.get_results()
    p = D / f"{name}.ecsv"; t.write(p, format="ascii.ecsv", overwrite=True)
    log(name, "https://gea.esac.esa.int/tap-server/tap", q, p)
    return t

todo = sys.argv[1:] or ["src", "nss", "nb", "ruwe", "ir"]
if "src0" in todo:
    tap("gaia_dr3_source", f"SELECT * FROM gaiadr3.gaia_source WHERE source_id = {SID}")
if "ap" in todo:
    tap("gaia_dr3_astrophysical_parameters", f"SELECT * FROM gaiadr3.astrophysical_parameters WHERE source_id = {SID}")
if "nss" in todo:
    for tb in ["nss_two_body_orbit", "nss_acceleration_astro", "nss_non_linear_spectro", "nss_vim_fl"]:
        tap(f"gaia_dr3_{tb}", f"SELECT * FROM gaiadr3.{tb} WHERE source_id = {SID}")
    tap("gaia_dr3_vari_summary", f"SELECT * FROM gaiadr3.vari_summary WHERE source_id = {SID}")
if "nb" in todo:
    tap("gaia_dr3_neighbours_30as",
        f"SELECT source_id, ra, dec, parallax, parallax_error, pmra, pmdec, phot_g_mean_mag, phot_bp_mean_mag, "
        f"phot_rp_mean_mag, ruwe, astrometric_params_solved, "
        f"DISTANCE(POINT(ra, dec), POINT({RA}, {DEC}))*3600 AS sep_arcsec "
        f"FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra, dec), CIRCLE({RA}, {DEC}, {30/3600})) ORDER BY sep_arcsec")
if "ruwe" in todo:
    # comparison sample: nearby (parallax>10 mas) stars with similar G and BP-RP, 5-parameter or 6-parameter solutions
    tap("ruwe_comparison_sample",
        "SELECT source_id, phot_g_mean_mag, bp_rp, parallax, ruwe, astrometric_excess_noise, astrometric_excess_noise_sig, "
        "ipd_gof_harmonic_amplitude, ipd_frac_multi_peak, astrometric_params_solved, non_single_star, "
        "radial_velocity_error, rv_nb_transits, rv_amplitude_robust, rv_chisq_pvalue, rv_renormalised_gof, grvs_mag, rv_template_teff "
        "FROM gaiadr3.gaia_source WHERE phot_g_mean_mag BETWEEN 10.2 AND 10.8 AND bp_rp BETWEEN 1.55 AND 1.85 "
        "AND parallax > 10")
if "ir" in todo:
    c = SkyCoord(RA * u.deg, DEC * u.deg)
    for cat, nm in [("II/246/out", "2mass_psc_30as"), ("II/328/allwise", "allwise_30as"),
                    ("II/363/unwise", "unwise_30as"), ("I/355/gaiadr3", "vizier_gaiadr3_30as")]:
        v = Vizier(columns=["**", "_r"], row_limit=-1)
        r = v.query_region(c, radius=30 * u.arcsec, catalog=cat)
        if len(r) == 0:
            p = D / f"{nm}.ecsv"; p.write_text("# no rows\n")
        else:
            t = r[0]; p = D / f"{nm}.ecsv"; t.write(p, format="ascii.ecsv", overwrite=True)
        log(nm, f"https://vizier.cds.unistra.fr/ catalog {cat}", f"query_region r=30arcsec {RA} {DEC}", p)
