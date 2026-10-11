"""Step 3: fresh catalogue queries (Gaia DR3/DR2 via ESA TAP, TIC 8.2 via MAST, 2MASS/AllWISE/LAMOST/others via
VizieR TAP/ASU). Raw responses are stored with provenance; nothing is parsed from older files."""
import os
import sys
import urllib.parse

sys.path.insert(0, os.path.dirname(__file__))
from audit_common import AUD, fetch, utcnow  # noqa

GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
SID = 3220388198192519424
RA, DEC = 80.85982665567, -1.32696401034


def gaia(q, dest, note):
    data = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q}).encode()
    open(os.path.join(AUD, dest.replace(".csv", ".adql")), "w").write(q + "\n")
    return fetch(GAIA_TAP, dest, note=note, data=data)


qs = {
    "gaia_dr3_source.csv": f"SELECT * FROM gaiadr3.gaia_source WHERE source_id={SID}",
    "gaia_dr3_ap.csv": f"SELECT * FROM gaiadr3.astrophysical_parameters WHERE source_id={SID}",
    "gaia_dr3_neighbours_180as.csv": (
        "SELECT source_id, ra, dec, ref_epoch, pmra, pmdec, pmra_error, pmdec_error, parallax, parallax_error, "
        "phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, bp_rp, ruwe, ipd_frac_multi_peak, "
        "ipd_gof_harmonic_amplitude, duplicated_source, astrometric_params_solved, phot_variable_flag, "
        "non_single_star, radial_velocity, "
        f"DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {RA}, {DEC}))*3600. AS sep_arcsec "
        "FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, dec), "
        f"CIRCLE('ICRS', {RA}, {DEC}, 0.05)) ORDER BY sep_arcsec"),
    "gaia_dr3_nss_two_body.csv": f"SELECT * FROM gaiadr3.nss_two_body_orbit WHERE source_id={SID}",
    "gaia_dr3_nss_accel.csv": f"SELECT * FROM gaiadr3.nss_acceleration_astro WHERE source_id={SID}",
    "gaia_dr3_nss_vim.csv": f"SELECT * FROM gaiadr3.nss_vim_fl WHERE source_id={SID}",
    "gaia_dr3_nss_nonlin.csv": f"SELECT * FROM gaiadr3.nss_non_linear_spectro WHERE source_id={SID}",
    "gaia_dr3_vari_summary.csv": f"SELECT * FROM gaiadr3.vari_summary WHERE source_id={SID}",
    "gaia_dr3_vari_classifier.csv": f"SELECT * FROM gaiadr3.vari_classifier_result WHERE source_id={SID}",
    "gaia_dr3_vari_rotation.csv": f"SELECT * FROM gaiadr3.vari_rotation_modulation WHERE source_id={SID}",
    "gaia_dr3_vari_eb.csv": f"SELECT * FROM gaiadr3.vari_eclipsing_binary WHERE source_id={SID}",
    "gaia_dr3_vari_planetary_transit.csv": f"SELECT * FROM gaiadr3.vari_planetary_transit WHERE source_id={SID}",
    "gaia_dr2_source.csv": f"SELECT * FROM gaiadr2.gaia_source WHERE source_id={SID}",
    "gaia_dr2_ruwe.csv": f"SELECT * FROM gaiadr2.ruwe WHERE source_id={SID}",
    "gaia_dr3_dr2_neighbourhood.csv": f"SELECT * FROM gaiadr3.dr2_neighbourhood WHERE dr3_source_id={SID}",
    "gaia_dr3_comoving_1deg.csv": (
        "SELECT source_id, ra, dec, parallax, parallax_error, pmra, pmdec, phot_g_mean_mag, bp_rp, radial_velocity, "
        f"DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {RA}, {DEC}))*3600. AS sep_arcsec "
        "FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, dec), "
        f"CIRCLE('ICRS', {RA}, {DEC}, 1.0)) AND parallax BETWEEN 24 AND 32"),
}
for fn, q in qs.items():
    try:
        r = gaia(q, "data/gaia/" + fn, "Gaia archive TAP " + fn)
        print(fn, r["bytes"], r["sha256"][:12])
    except Exception as e:  # noqa
        print("FAILED", fn, e)

# TIC 8.2 cone (MAST catalogs API)
from astroquery.mast import Catalogs  # noqa

t0 = utcnow()
tic = Catalogs.query_region(f"{RA} {DEC}", radius=0.0834, catalog="TIC")  # 5 arcmin
p = os.path.join(AUD, "data/catalogs/tic82_cone5arcmin.csv")
tic.write(p, format="csv", overwrite=True)
from audit_common import log_prov, sha256_file  # noqa

log_prov(dict(url="astroquery.mast.Catalogs.query_region(catalog='TIC', radius=0.0834 deg)", dest=os.path.relpath(p, AUD),
              utc_start=t0, utc_end=utcnow(), bytes=os.path.getsize(p), sha256=sha256_file(p), note="TIC 8.2 cone 5'"))
print("TIC rows", len(tic))

# VizieR ASU (tab-separated) for photometric / spectroscopic catalogues
VZ = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?"
viz = {
    "2mass_psc.tsv": ("II/246/out", 10),
    "allwise.tsv": ("II/328/allwise", 10),
    "catwise2020.tsv": ("II/365/catwise", 10),
    "tycho2.tsv": ("I/259/tyc2", 10),
    "apass_dr9.tsv": ("II/336/apass9", 10),
    "lamost_dr5_lrs.tsv": ("V/164/dr5", 10),
    "lamost_dr7_mrs.tsv": ("V/156/dr7mrs", 10),
    "lamost_mrs_general.tsv": ("V/156/dr7lrs", 10),
    "tic82_vizier.tsv": ("IV/39/tic82", 10),
    "ctl_v8.tsv": ("IV/38/tic", 10),
    "lepine_gaidos2011.tsv": ("J/AJ/142/138/Mdwarfs", 10),
    "galex_gr6.tsv": ("II/335/galex_ais", 15),
    "rosat_2rxs.tsv": ("J/A+A/588/A103/cat2rxs", 60),
    "erosita_dr1_main.tsv": ("J/A+A/682/A34/erass1-m", 30),
    "vsx.tsv": ("B/vsx/vsx", 60),
    "wds.tsv": ("B/wds/wds", 60),
    "gaia_dr3_vizier.tsv": ("I/355/gaiadr3", 10),
    "bailer_jones2021.tsv": ("I/352/gedr3dis", 10),
}
for fn, (cat, rad) in viz.items():
    url = VZ + urllib.parse.urlencode({"-source": cat, "-c": f"{RA} {DEC}", "-c.rs": rad, "-out.all": "",
                                        "-out.max": "unlimited", "-out.form": "|"})
    try:
        r = fetch(url, "data/catalogs/" + fn, note=f"VizieR {cat} r={rad}arcsec")
        n = sum(1 for line in open(os.path.join(AUD, "data/catalogs/" + fn)) if line.strip() and not line.startswith("#"))
        print(fn, cat, r["bytes"], "lines", n)
    except Exception as e:  # noqa
        print("FAILED", fn, e)
