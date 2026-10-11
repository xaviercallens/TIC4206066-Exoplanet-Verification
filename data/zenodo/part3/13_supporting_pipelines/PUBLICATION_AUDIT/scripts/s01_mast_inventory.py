"""Step 1: fresh MAST inventory of every TESS/HLSP product for TIC 4206066, then download the time-series
products used anywhere in the investigation, record sha256 and compare with the files used earlier."""
import json
import os
import sys

import numpy as np
from astroquery.mast import Observations

sys.path.insert(0, os.path.dirname(__file__))
from audit_common import AUD, OLD_DEEP, fetch, save_json, sha256_file, utcnow  # noqa

RA, DEC = 80.85982665567, -1.32696401034

t0 = utcnow()
obs = Observations.query_region(f"{RA} {DEC}", radius="0.005 deg")
obs = obs[np.isin(obs["obs_collection"], ["TESS", "HLSP"])]
rows = []
for r in obs:
    rows.append({k: (str(r[k]) if not isinstance(r[k], (int, float, np.integer, np.floating)) else float(r[k]))
                 for k in ["obs_collection", "provenance_name", "target_name", "obs_id", "sequence_number",
                           "t_min", "t_max", "t_exptime", "dataproduct_type", "s_ra", "s_dec", "dataRights", "obsid",
                           "t_obs_release"]})
save_json({"query_utc": t0, "center": [RA, DEC], "radius_deg": 0.005, "n": len(rows), "obs": rows},
          "data/mast/mast_observations.json")

# keep only observations of THIS star (target name = TIC id or Gaia id)
mine = [r for r in rows if ("4206066" in r["target_name"]) or ("3220388198192519424" in r["target_name"])
        or r["obs_collection"] == "TESS" and r["provenance_name"] in ("SPOC", "TESS-SPOC") and "4206066" in r["obs_id"]]
print("observations near star:", len(rows), " of this star:", len(mine))
for r in sorted(mine, key=lambda x: (x["provenance_name"], x["obs_id"])):
    print(f'{r["obs_collection"]:5s} {r["provenance_name"]:18s} {r["target_name"]:22s} {r["obs_id"][:90]:90s} '
          f'seq={r["sequence_number"]} exptime={r["t_exptime"]} release={r["t_obs_release"]}')

obsids = [r["obsid"] for r in mine]
prods = Observations.get_product_list([str(int(float(o))) for o in obsids])
plist = []
for p in prods:
    plist.append({k: str(p[k]) for k in ["obsID", "obs_collection", "productFilename", "dataURI", "size",
                                          "description", "productType", "productSubGroupDescription",
                                          "project", "prvversion", "dataRights"] if k in p.colnames})
save_json({"query_utc": utcnow(), "n": len(plist), "products": plist}, "data/mast/mast_products.json")
print("products:", len(plist))

# download light curves / target pixel files (science time series only)
want = []
for p in plist:
    fn = p["productFilename"]
    if not fn.endswith(".fits"):
        continue
    if any(s in fn for s in ["_lc.fits", "_tp.fits", "_llc.fits", "lc.fits", "_tp.fits"]) or fn.endswith("_lc.fits"):
        want.append(p)
    elif "dvt" in fn or "dvr" in fn or "dvm" in fn:
        want.append(p)
print("to download:", len(want))
old_sha = {}
dl = os.path.join(OLD_DEEP, "lc", "download_log.json")
for rec in json.load(open(dl)):
    old_sha[rec["file"]] = rec["sha256"]
res = []
for p in want:
    fn = p["productFilename"]
    url = "https://mast.stsci.edu/api/v0.1/Download/file?uri=" + p["dataURI"]
    dest = os.path.join("data", "mast", fn)
    if os.path.exists(os.path.join(AUD, dest)):
        rec = {"dest": dest, "sha256": sha256_file(os.path.join(AUD, dest)), "url": url, "cached": True}
    else:
        rec = fetch(url, dest, note="MAST product " + p["dataURI"])
    rec["old_sha256"] = old_sha.get(fn)
    rec["same_as_old"] = (rec["old_sha256"] == rec["sha256"]) if rec["old_sha256"] else None
    rec["productFilename"] = fn
    rec["dataURI"] = p["dataURI"]
    res.append(rec)
    print(f'{fn:95s} {rec["sha256"][:12]} old={str(rec["old_sha256"])[:12]} same={rec["same_as_old"]}')
save_json(res, "data/mast/download_check.json")
