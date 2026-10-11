"""(c) TESS PRF model FITS (start_s0004 set, valid for Sector >= 4 per the MAST 00README.txt) for the needed camera/CCDs."""
import os, sys, re, json
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import requests
from astropy.io import fits

BASE = "https://archive.stsci.edu/missions/tess/models/prf_fitsfiles/"
need = sys.argv[1:] or ["cam1_ccd1", "cam1_ccd2"]
out = []
# READMEs / directory indexes as documentation snapshots
for rel in ["", "start_s0001/", "start_s0004/", "start_s0004/00README.txt", "start_s0001/00README.txt"]:
    fn = "prf_index_" + (rel.replace("/", "_") or "root") + (".html" if not rel.endswith(".txt") else "")
    e = http_get(BASE + rel, os.path.join(ROOT, "docs", "prf", fn), note="PRF archive doc/index snapshot")
    print(e.get("relpath"), e.get("http_status"), e.get("error"))
for cc in need:
    idx = requests.get(BASE + f"start_s0004/{cc}/", headers=UA, timeout=120).text
    files = sorted(set(re.findall(r'href="([^"]+-row\d{4}-col\d{4}\.fits)"', idx)))
    print(cc, len(files))
    for f in files:
        e = http_get(BASE + f"start_s0004/{cc}/{f}", os.path.join(ROOT, "prf", "start_s0004", cc, f), note="TESS PRF model FITS, start_s0004")
        h = fits.getheader(os.path.join(ROOT, e["relpath"]))
        e2 = dict(e, prf_header={k: h.get(k) for k in ("CAMERA", "CCD", "VERSION", "PRF_ROW", "PRF_COL", "CDELT1P", "CDELT2P", "NAXIS1", "NAXIS2") if k in h})
        out.append(e2)
json.dump(out, open(os.path.join(ROOT, "prf", "prf_download_log.json"), "w"), indent=1, default=str)
print("n files", len(out), "errors", sum(1 for e in out if "error" in e))
