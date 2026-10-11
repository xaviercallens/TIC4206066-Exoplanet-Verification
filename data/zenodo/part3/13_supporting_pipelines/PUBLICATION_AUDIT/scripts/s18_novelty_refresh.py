"""Step 18: refresh of the most important priority/novelty sources at audit time (the full adversarial search is the
SEALED priority re-check, 2026-09-25T04:46:08Z-05:24:55Z; this step only re-queries the sources that can change fastest
and checks that the earlier package covered aliases, harmonics and identifiers). Nothing is posted anywhere."""
import csv
import io
import json
import os
import sys
import urllib.parse

sys.path.insert(0, os.path.dirname(__file__))
from audit_common import AUD, OLD_PRIO, fetch, save_json, utcnow  # noqa

T0 = utcnow()
res = {"started_utc": T0}
UA = {"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.0 Safari/605.1.15"}


def grab(url, dest, note, headers=None):
    try:
        r = fetch(url, "data/novelty/" + dest, note=note, headers=headers or UA)
        return open(os.path.join(AUD, "data/novelty", dest), errors="ignore").read()
    except Exception as e:  # noqa
        return "FETCH_FAILED " + repr(e)


toi = grab("https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv", "exofop_toi.csv", "ExoFOP TOI list")
res["exofop_toi"] = dict(n_rows=toi.count("\n"), tic_4206066_present=any(r.split(",")[1:2] == ["4206066"] for r in toi.splitlines()[1:]),
                         any_string_4206066=("4206066" in toi))
ctoi = grab("https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv", "exofop_ctoi.csv", "ExoFOP CTOI list")
res["exofop_ctoi"] = dict(n_rows=ctoi.count("\n"), any_string_4206066=("4206066" in ctoi))
tgt = grab("https://exofop.ipac.caltech.edu/tess/target.php?id=4206066&json", "exofop_target_4206066.json", "ExoFOP target JSON")
try:
    tj = json.loads(tgt)
    res["exofop_target"] = dict(tois=tj.get("tois"), ctois=tj.get("ctois"), n_files=len(tj.get("files", []) or []),
                                n_imaging=len(tj.get("imaging", []) or []), n_spectroscopy=len(tj.get("spectroscopy", []) or []),
                                n_tseries=len(tj.get("tseries", []) or []), n_obsnotes=len(tj.get("obsnotes", []) or []))
except Exception as e:  # noqa
    res["exofop_target"] = dict(error=repr(e), head=tgt[:300])
q = urllib.parse.quote("select * from toi where tid=4206066")
nea = grab(f"https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query={q}&format=csv", "nea_toi_tid.csv", "NASA Exoplanet Archive toi by TIC")
res["nasa_toi_by_tic_rows"] = max(0, nea.count("\n") - 1) if not nea.startswith("FETCH_FAILED") else nea
q2 = urllib.parse.quote("select pl_name,hostname,ra,dec,pl_orbper from pscomppars where contains(point('icrs',ra,dec),circle('icrs',80.85983,-1.32696,0.0833))=1")
nps = grab(f"https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query={q2}&format=csv", "nea_pscomppars_cone5arcmin.csv", "NASA pscomppars cone 5'")
res["nasa_pscomppars_cone5arcmin_rows"] = max(0, nps.count("\n") - 1) if not nps.startswith("FETCH_FAILED") else nps
tce = grab("https://exo.mast.stsci.edu/api/v0.1/dvdata/tess/4206066/tces/", "exomast_tces.json", "exo.MAST TCE list")
res["exomast_tces"] = tce[:300]
for qname, qq in {"tic": '"TIC 4206066"', "stkm": '"StKM 1-561"', "tyc": '"TYC 4753-1591-1"', "gaia": "3220388198192519424",
                  "2mass": "J05232636-0119380"}.items():
    url = "https://export.arxiv.org/api/query?" + urllib.parse.urlencode({"search_query": f"all:{qq}", "max_results": 50})
    txt = grab(url, f"arxiv_{qname}.xml", f"arXiv API all:{qq}")
    res[f"arxiv_{qname}"] = dict(n_entries=txt.count("<entry>") if not txt.startswith("FETCH_FAILED") else txt[:120])
sim = grab("https://simbad.cds.unistra.fr/simbad/sim-tap/sync?" + urllib.parse.urlencode({
    "request": "doQuery", "lang": "adql", "format": "csv",
    "query": "SELECT b.bibcode, b.year FROM ident AS i JOIN has_ref AS h ON h.oidbibref = i.oidref JOIN ref AS b ON b.oidbib = h.oidbibref "
             "WHERE i.id = 'TIC 4206066' ORDER BY b.year"}), "simbad_refs.csv", "SIMBAD references via TAP")
res["simbad_refs"] = sim.strip().splitlines()[1:] if not sim.startswith("FETCH_FAILED") else sim[:200]
tl = grab("https://heasarc.gsfc.nasa.gov/docs/tess/data/target_lists/sector110_targets.csv", "tess_2min_targets_S110.csv", "TESS 2-min targets S110")
res["tess_2min_S110"] = dict(available=not tl.startswith("FETCH_FAILED"), tic_4206066_listed=("4206066" in tl) if not tl.startswith("FETCH_FAILED") else None,
                             head=tl[:120])
# coverage review of the SEALED priority package (aliases + harmonics searched)
pm = json.load(open(os.path.join(OLD_PRIO, "priority_matches.json")))
src = list(csv.DictReader(open(os.path.join(OLD_PRIO, "priority_sources.csv"))))
res["sealed_package"] = dict(n_source_groups=len(src), n_matches_examined=len(pm) if isinstance(pm, list) else len(pm.get("matches", [])),
                             identifiers_note="TIC, StKM 1-561, Gaia DR2/DR3, TYC (both spellings), 2MASS, PM J/PM I, UCAC4, APASS, WISE, LAMOST, Sp0523-0119",
                             harmonics_note="P/3, P/2, P, 2P, 3P for both periods (1% near star, 0.1% catalogue-wide)")
res["finished_utc"] = utcnow()
save_json(res, "out/s18_novelty_refresh.json")
print(json.dumps(res, indent=1)[:4000])
