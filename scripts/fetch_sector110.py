#!/usr/bin/env python3
"""
Fetch TESS Sector N light curves for TIC 4206066 from MAST and run the frozen
pre-registered Sector 110 decision logic.

Fail-closed: if no matching product is published, exit code 2 ("AWAITING_DATA")
and nothing is evaluated. Every downloaded file is SHA-256 hashed into
results/sector110_fetch_manifest.json BEFORE evaluation, so the exact input is on record.

  python scripts/fetch_sector110.py --list                 # what MAST has (no download)
  python scripts/fetch_sector110.py --fetch                # download + hash
  python scripts/fetch_sector110.py --evaluate             # DISABLED: refuses (see below)
  python scripts/fetch_sector110.py --list --sector 6      # validate the query on a known sector

The REGISTERED analysis is the frozen code in the Zenodo deposit doi:10.5281/zenodo.23175179
(s110_final.py, run in its pinned arm64-macOS environment). The thresholds in
subprojects/tic4206066_sector110_runner.py are NOT the registered rules, so --evaluate refuses.
This script only lists/downloads/hashes. Do NOT open the FITS before the registered primary run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from typing import Any

import requests

MAST_API = "https://mast.stsci.edu/api/v0/invoke"
MAST_DOWNLOAD = "https://mast.stsci.edu/api/v0.1/Download/file"
TIC_ID = "4206066"
MAST_DIR = "data/mast"
MANIFEST = "results/sector110_fetch_manifest.json"
OUTCOME = "results/sector110_outcome.json"
EXIT_AWAITING = 2

# Per the pre-registration (Appendix G): the primary product is the first released SPOC 2-min light
# curve (mast filename tess...-s0110-...-s_lc.fits). QLP/TESS-SPOC fallbacks apply only if no SPOC
# 2-min product is public on 2027-02-28, so they are deliberately NOT selectable here.
PIPELINE_PRIORITY = ("tess",)  # MAST SPOC 2-min files: tess<date>-s<sector>-<tic>-...-s_lc.fits


def select_product(prods: list[dict[str, Any]]) -> dict[str, Any] | None:
    """Pick exactly one product by the frozen priority order; ignore all other pipelines."""
    for prefix in PIPELINE_PRIORITY:
        hits = sorted((p for p in prods if p["filename"].lower().startswith(prefix)
                       and p["filename"].lower().endswith("-s_lc.fits")), key=lambda p: p["filename"])
        if hits:
            return hits[0]
    return None


def _mast(service: str, params: dict[str, Any]) -> list[dict[str, Any]]:
    body = {"service": service, "params": params, "format": "json", "timeout": 60}
    r = requests.post(MAST_API, data={"request": json.dumps(body)}, timeout=90)
    r.raise_for_status()
    return r.json().get("data", [])


def list_products(sector: int) -> list[dict[str, Any]]:
    """Light-curve FITS products (SPOC 'lc', QLP 'llc') for this target and sector."""
    tag = f"s{sector:04d}"
    obs = _mast(
        "Mast.Caom.Filtered",
        {
            "columns": "obsid,obs_collection,target_name,sequence_number,obs_id",
            "filters": [{"paramName": "target_name", "values": [TIC_ID]}],
        },
    )
    out: list[dict[str, Any]] = []
    for o in obs:
        if tag not in str(o.get("obs_id", "")).lower() and o.get("sequence_number") != sector:
            continue
        for p in _mast("Mast.Caom.Products", {"obsid": o["obsid"]}):
            name = str(p.get("productFilename", ""))
            if name.endswith(".fits") and (name.endswith("_lc.fits") or name.endswith("_llc.fits")):
                if tag in name.lower() or o.get("sequence_number") == sector:
                    out.append({"filename": name, "uri": p["dataURI"], "size": p.get("size"),
                                "collection": o.get("obs_collection")})
    return out


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(sector: int) -> int:
    chosen = select_product(list_products(sector))
    if chosen is None:
        print(f"AWAITING_DATA: no Sector {sector} light curve from a pre-approved pipeline "
              f"{PIPELINE_PRIORITY} published on MAST for TIC {TIC_ID}.")
        return EXIT_AWAITING
    prods = [chosen]
    os.makedirs(MAST_DIR, exist_ok=True)
    entries = []
    for p in prods:
        dest = os.path.join(MAST_DIR, p["filename"])
        with requests.get(MAST_DOWNLOAD, params={"uri": p["uri"]}, stream=True, timeout=300) as r:
            r.raise_for_status()
            with open(dest, "wb") as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk)
        entries.append({**p, "path": dest, "sha256": sha256_file(dest), "bytes": os.path.getsize(dest)})
        print(f"downloaded {p['filename']} sha256={entries[-1]['sha256']}")
    os.makedirs("results", exist_ok=True)
    with open(MANIFEST, "w") as f:
        json.dump({"sector": sector, "tic": TIC_ID, "fetched_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
                   "pipeline_priority": list(PIPELINE_PRIORITY), "files": entries}, f, indent=2)
    return 0


def evaluate() -> int:
    print("REFUSED: this runner's rules are not the pre-registered ones (doi:10.5281/zenodo.23175179). "
          "Run the frozen s110_final.py in its pinned environment; results here would not count.")
    return 4


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sector", type=int, default=110)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--list", action="store_true")
    g.add_argument("--fetch", action="store_true")
    g.add_argument("--evaluate", action="store_true")
    a = ap.parse_args()
    if a.list:
        prods = list_products(a.sector)
        for p in prods:
            print(p["collection"], p["filename"], p["size"])
        if not prods:
            print(f"AWAITING_DATA: no Sector {a.sector} light curve on MAST for TIC {TIC_ID}.")
            return EXIT_AWAITING
        return 0
    return fetch(a.sector) if a.fetch else evaluate()


if __name__ == "__main__":
    sys.exit(main())
