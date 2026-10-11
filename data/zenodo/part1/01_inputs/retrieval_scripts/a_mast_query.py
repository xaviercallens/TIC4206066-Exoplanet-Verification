"""(a) List ALL MAST observations/products for TIC 4206066 (incl. HLSPs)."""
import os, sys, json
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.mast import Observations
import astroquery

out = os.path.join(ROOT, "mast")
c = SkyCoord(RA, DEC, unit="deg")
meta = dict(utc=utcnow(), astroquery=astroquery.__version__)

# 1) cone search: all observations whose footprint contains the target, radius 10 arcsec
obs_cone = Observations.query_region(c, radius=10 * u.arcsec)
# 2) name-based: TIC id and Gaia DR3 id (T16/TGLC use Gaia ids) -- catch anything the footprint misses
names = ["4206066", "TIC 4206066", "TIC4206066", "3220388198192519424", "03220388198192519424"]
obs_name = Observations.query_criteria(target_name=names)
from astropy.table import vstack, unique
allobs = vstack([obs_cone, obs_name], metadata_conflicts="silent")
allobs = unique(allobs, keys="obsid")
meta["n_obs_cone10as"] = len(obs_cone)
meta["n_obs_name"] = len(obs_name)
meta["n_obs_union"] = len(allobs)
cols = ["obsid", "obs_collection", "provenance_name", "project", "instrument_name", "dataproduct_type", "target_name",
        "obs_id", "sequence_number", "t_exptime", "s_ra", "s_dec", "t_min", "t_max", "dataURL", "calib_level"]
cols = [k for k in cols if k in allobs.colnames]
allobs[cols].write(os.path.join(out, "mast_observations_all.csv"), overwrite=True)

# keep only time-series / TESS / HLSP; FFI sector-level image obs (huge footprints) are listed but not expanded
ts = allobs[(allobs["dataproduct_type"] == "timeseries")]
meta["n_obs_timeseries"] = len(ts)
prods = Observations.get_product_list(ts)
pc = [k for k in ["obsID", "obs_collection", "project", "provenance_name" if "provenance_name" in prods.colnames else "productGroupDescription",
                  "productType", "productSubGroupDescription", "description", "productFilename", "dataURI", "size", "dataRights"] if k in prods.colnames]
prods[pc].write(os.path.join(out, "mast_products_all.csv"), overwrite=True)
meta["n_products"] = len(prods)
# summary by collection/provenance
summ = {}
obs_prov = {str(r["obsid"]): (str(r["obs_collection"]), str(r["provenance_name"]), str(r["target_name"]), r["sequence_number"]) for r in allobs}
for r in allobs:
    key = f'{r["obs_collection"]}|{r["provenance_name"]}|{r["dataproduct_type"]}'
    summ.setdefault(key, []).append(dict(obsid=str(r["obsid"]), target=str(r["target_name"]), obs_id=str(r["obs_id"]),
                                          seq=None if np.ma.is_masked(r["sequence_number"]) else int(r["sequence_number"]),
                                          texp=float(r["t_exptime"]) if not np.ma.is_masked(r["t_exptime"]) else None))
meta["summary"] = summ
json.dump(meta, open(os.path.join(out, "mast_query_summary.json"), "w"), indent=1, default=str)
print(json.dumps({k: v for k, v in meta.items() if k != "summary"}, indent=1))
for k, v in sorted(summ.items()):
    print(k, len(v), [ (x["obs_id"], x["seq"], x["texp"]) for x in v][:12])
