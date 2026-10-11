"""(g) Official documentation snapshots (raw file + text extraction), with URL, UTC, sha256."""
import os, sys, json, subprocess, re, html
sys.path.insert(0, os.path.dirname(__file__))
from common import *

D = os.path.join(ROOT, "docs")
DOCS = [
    # SPOC data products / quality flags
    ("spoc/EXP-TESS-ARC-ICD-TM-0014-Rev-F.pdf", "https://archive.stsci.edu/files/live/sites/mast/files/home/missions-and-data/active-missions/tess/_documents/EXP-TESS-ARC-ICD-TM-0014-Rev-F.pdf", "TESS Science Data Products Description Document (SDPDD) Rev F: QUALITY bit definitions"),
    ("spoc/EXP-TESS-ARC-ICD-TM-0014-Rev-F-FINAL.pdf", "https://archive.stsci.edu/missions/tess/doc/EXP-TESS-ARC-ICD-TM-0014-Rev-F-FINAL.pdf", "SDPDD Rev F (archive mission doc directory copy)"),
    ("spoc/mast_tess_archive_manual_2.0_data_product_overview.html", "https://outerspace.stsci.edu/spaces/TESS/pages/14563420/2.0+-+Data+Product+Overview", "MAST TESS Archive Manual 2.0 Data Product Overview (PRF start_s0001 vs start_s0004)"),
    ("spoc/mast_tess_data_products.html", "https://archive.stsci.edu/missions-and-data/tess/data-products", "MAST TESS data products page"),
    ("spoc/heasarc_observing_technical.html", "https://heasarc.gsfc.nasa.gov/docs/tess/observing-technical.html", "HEASARC TESS observing-technical (PSF/PRF section)"),
    # data release notes
    ("drn/tess_sector_06_drn08_v02.pdf", "https://archive.stsci.edu/missions/tess/doc/tess_drn/tess_sector_06_drn08_v02.pdf", "TESS DRN Sector 6"),
    ("drn/tess_sector_32_drn48_v02.pdf", "https://archive.stsci.edu/missions/tess/doc/tess_drn/tess_sector_32_drn48_v02.pdf", "TESS DRN Sector 32"),
    ("drn/tess_sector_98_drn134_v01.pdf", "https://archive.stsci.edu/missions/tess/doc/tess_drn/tess_sector_98_drn134_v01.pdf", "TESS DRN Sector 98"),
    ("drn/tess_drn_index.html", "https://archive.stsci.edu/missions/tess/doc/tess_drn/", "TESS DRN directory index"),
    ("drn/heasarc_condensed_data_release_notes.html", "https://heasarc.gsfc.nasa.gov/docs/tess/data_release_notes.html", "HEASARC condensed TESS DRN"),
    ("heasarc/sector98_summary.html", "https://heasarc.gsfc.nasa.gov/docs/tess/sector98_summary.html", "HEASARC Sector 98 summary (scattered light notes)"),
    ("heasarc/sector.html", "https://heasarc.gsfc.nasa.gov/docs/tess/sector.html", "HEASARC TESS sector table"),
    # QLP
    ("qlp/mast_hlsp_qlp.html", "https://archive.stsci.edu/hlsp/qlp", "MAST QLP HLSP page"),
    ("qlp/hlsp_qlp_tess_ffi_all_tess_v1_readme.txt", "https://archive.stsci.edu/hlsps/qlp/hlsp_qlp_tess_ffi_all_tess_v1_readme.txt", "QLP HLSP readme"),
    ("qlp/hlsp_qlp_tess_ffi_all_tess_v1_data-prod-desc.pdf", "https://archive.stsci.edu/hlsps/qlp/hlsp_qlp_tess_ffi_all_tess_v1_data-prod-desc.pdf", "QLP data product description (quality flags)"),
    ("qlp/mit_qlp.html", "https://tess.mit.edu/qlp/", "MIT QLP page"),
    # TESScut / astrocut
    ("tesscut/mast_tesscut_docs.html", "https://mast.stsci.edu/tesscut/docs/getting_started.html", "TESScut API docs"),
]
res = []
for rel, url, note in DOCS:
    dest = os.path.join(D, rel)
    e = http_get(url, dest, note=note, timeout=180)
    if "error" not in e:
        # text extraction for searchability
        txt = dest + ".txt"
        if dest.endswith(".pdf"):
            subprocess.run(["pdftotext", "-layout", dest, txt], check=False)
        else:
            raw = open(dest, "rb").read().decode("utf-8", "replace")
            raw = re.sub(r"(?is)<(script|style)[^>]*>.*?</\1>", " ", raw)
            t = html.unescape(re.sub(r"<[^>]+>", " ", raw))
            t = "\n".join(l.strip() for l in t.splitlines() if l.strip())
            open(txt, "w").write(t)
        if os.path.exists(txt):
            e["text_extract"] = os.path.relpath(txt, ROOT)
            e["text_sha256"] = sha256_file(txt)
    res.append(e)
    print(rel, e.get("http_status"), e.get("bytes"), e.get("error"), flush=True)
json.dump(res, open(os.path.join(D, "docs_download_log.json"), "w"), indent=1, default=str)
