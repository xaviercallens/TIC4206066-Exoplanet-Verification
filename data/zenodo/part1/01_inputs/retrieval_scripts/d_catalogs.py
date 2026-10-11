"""(d)+(e) Gaia DR3 (ESA TAP), Bailer-Jones 2021, TIC 8.2 cone, VizieR photometry. Saves ECSV + CSV + query text."""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import astroquery
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
from astroquery.mast import Catalogs

OUT = os.path.join(ROOT, "catalogs")
os.makedirs(OUT, exist_ok=True)
c = SkyCoord(RA, DEC, unit="deg")
log = []

def save(tab, name, url, query, extra=None):
    p1 = os.path.join(OUT, name + ".ecsv"); p2 = os.path.join(OUT, name + ".csv")
    tab.write(p1, overwrite=True, format="ascii.ecsv")
    tab.write(p2, overwrite=True, format="csv")
    if query:
        open(os.path.join(OUT, name + ".query.txt"), "w").write(query + "\n")
    for p in (p1, p2):
        e = dict(relpath=os.path.relpath(p, ROOT), url=url, utc=utcnow(), bytes=os.path.getsize(p), sha256=sha256_file(p),
                 query=query, nrows=len(tab), ncols=len(tab.colnames), tool=f"astroquery {astroquery.__version__}", **(extra or {}))
        log_prov(e); log.append(e)
    print(name, len(tab), "rows", len(tab.colnames), "cols", flush=True)

GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap (astroquery.gaia GaiaClass sync job)"
def gaia(name, q):
    for k in range(4):
        try:
            job = Gaia.launch_job(q, dump_to_file=False)
            t = job.get_results()
            save(t, name, GAIA_TAP, q); return t
        except Exception as ex:
            print("gaia retry", name, repr(ex)[:200], flush=True); time.sleep(10 * (k + 1))
    log.append(dict(name=name, error="gaia query failed", query=q))

cone = (f"SELECT *, DISTANCE(POINT({RA},{DEC}), POINT(ra,dec))*3600. AS dist_arcsec_j2016pos "
        f"FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec), CIRCLE({RA},{DEC},{240/3600.:.10f})) ORDER BY dist_arcsec_j2016pos")
gaia("gaiadr3_gaia_source_cone240as", cone)
for tab in ["nss_two_body_orbit", "nss_acceleration_astro", "nss_non_linear_spectro", "nss_vim_fl", "vari_summary",
            "astrophysical_parameters"]:
    gaia(f"gaiadr3_{tab}_target", f"SELECT * FROM gaiadr3.{tab} WHERE source_id = {GAIA_DR3}")

# Bailer-Jones et al. 2021 (VizieR I/352/gedr3dis)
V = Vizier(columns=["**", "+_r"], row_limit=-1)
def viz(name, cat, radius_as, **kw):
    for k in range(4):
        try:
            if "constraints" in kw:
                r = V.query_constraints(catalog=cat, **kw["constraints"])
            else:
                r = V.query_region(c, radius=radius_as * u.arcsec, catalog=cat)
            break
        except Exception as ex:
            print("vizier retry", name, repr(ex)[:200], flush=True); time.sleep(10 * (k + 1))
    else:
        log.append(dict(name=name, error="vizier failed")); return
    if len(r) == 0:
        print(name, "NO ROWS", flush=True); log.append(dict(name=name, catalog=cat, nrows=0, utc=utcnow())); return
    for i, t in enumerate(r):
        q = f"VizieR {r.keys()[i]} " + (f"source_id constraint {kw['constraints']}" if "constraints" in kw else f"cone r={radius_as} arcsec around ({RA},{DEC}) J2000")
        save(t, name + ("" if len(r) == 1 else f"_{i}"), f"{Vizier.VIZIER_SERVER} (astroquery Vizier; table {r.keys()[i]})", q,
             extra=dict(vizier_table=r.keys()[i]))

viz("bailerjones2021_I352", "I/352/gedr3dis", None, constraints={"Source": str(GAIA_DR3)})
viz("2mass_psc_II246", "II/246/out", 30)
viz("allwise_II328", "II/328/allwise", 30)
viz("catwise2020_II365", "II/365/catwise", 30)
viz("tycho2_I259", "I/259/tyc2", 30)
viz("apass_dr9_II336", "II/336/apass9", 30)
viz("galex_gr67_ais_II335", "II/335/galex_ais", 30)

# TIC 8.2 cone 300 arcsec (MAST Catalogs, all columns)
for k in range(4):
    try:
        t = Catalogs.query_region(c, radius=300 * u.arcsec, catalog="TIC")
        break
    except Exception as ex:
        print("tic retry", repr(ex)[:200]); time.sleep(10 * (k + 1))
save(t, "tic82_cone300as", "https://mast.stsci.edu/api/v0/invoke (astroquery.mast.Catalogs, catalog=TIC)",
     f"Catalogs.query_region(({RA},{DEC}), radius=300 arcsec, catalog='TIC')")
json.dump(log, open(os.path.join(OUT, "catalog_query_log.json"), "w"), indent=1, default=str)
print("DONE")
