"""(a) Download the required MAST time-series products (direct MAST Download API, one file per request)."""
import os, sys, json, csv
sys.path.insert(0, os.path.dirname(__file__))
from common import *

rows = list(csv.DictReader(open(os.path.join(ROOT, "mast", "mast_products_all.csv"))))
WANT = [  # required by the audit task
    "hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits",
    "hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_tp.fits",
    "tess2020324010417-s0032-0000000004206066-0200-s_lc.fits",
    "tess2020324010417-s0032-0000000004206066-0200-s_tp.fits",
    "hlsp_tess-spoc_tess_phot_0000000004206066-s0032_tess_v1_lc.fits",
    "hlsp_tess-spoc_tess_phot_0000000004206066-s0032_tess_v1_tp.fits",
    "hlsp_qlp_tess_ffi_s0006-0000000004206066_tess_v01_llc.fits",
    "hlsp_qlp_tess_ffi_s0032-0000000004206066_tess_v01_llc.fits",
    "hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits",
    # small optional HLSPs requested "if small"
    "hlsp_gsfc-eleanor-lite_tess_ffi_s0006-0000000004206066_tess_v1.0_lc.fits",
    "hlsp_tglc_tess_ffi_gaiaid-3220388198192519424-s0006-cam1-ccd1_tess_v1_llc.fits",
]
BASE = "https://mast.stsci.edu/api/v0.1/Download/file"
listed = []
for r in rows:
    fn = r["productFilename"]
    url = f"{BASE}?uri={r['dataURI']}"
    listed.append(dict(file=fn, collection=r["obs_collection"], project=r["project"], uri=r["dataURI"], url=url,
                       mast_size=int(r["size"]), downloaded=fn in WANT))
    if fn in WANT:
        e = http_get(BASE, os.path.join(ROOT, "mast", fn), params={"uri": r["dataURI"]}, note=f"MAST product {r['dataURI']}")
        e["mast_listed_size"] = int(r["size"])
        print(json.dumps(e), flush=True)
missing = [w for w in WANT if w not in [x["file"] for x in listed]]
json.dump(dict(listed=listed, missing_from_mast=missing), open(os.path.join(ROOT, "mast", "mast_products_listing_urls.json"), "w"), indent=1)
print("missing", missing)
