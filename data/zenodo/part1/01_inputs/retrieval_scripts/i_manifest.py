"""Build DATA_MANIFEST.json: one entry per file (relpath, url, utc, sha256, bytes, key header info, matches_existing)."""
import os, sys, json, re, glob
sys.path.insert(0, os.path.dirname(__file__))
from common import *

A = "/Users/paulsmith/space <3/tess-jev/addenda"
prov = {}
for line in open(PROV):
    e = json.loads(line)
    if "error" in e:
        continue
    prov[e["relpath"]] = e  # last successful entry wins
H = json.load(open(os.path.join(ROOT, "logs", "fits_headers.json")))
C = json.load(open(os.path.join(ROOT, "logs", "comparison_with_existing.json")))
# exact TESScut URLs as issued by astroquery (from its own log line)
tc_url = {}
for line in open(os.path.join(ROOT, "logs", "b_tesscut.log")):
    m = re.match(r"Downloading URL (\S+) to .*/tesscut/(\d+x\d+)/", line)
    if m:
        u = m.group(1); sec = int(re.search(r"sector=(\d+)", u).group(1)); tc_url[(m.group(2), sec)] = u
KEYH = ["SECTOR", "CAMERA", "CCD", "TSTART", "TSTOP", "TIMESYS", "TIMEREF", "TIMEUNIT", "TELAPSE", "TELAPS", "BJDREFI", "BJDREFF", "BJDREFR",
        "INT_TIME", "EXPOSURE", "TIMEDEL", "NUM_FRM", "CROWDSAP", "FLFRCSAP", "ORIGIN", "CREATOR", "PROCVER", "DATA_REL", "DATE", "CALIB", "BESTAP",
        "FFI_TYPE", "TICVER", "VERSION", "1CRV4P", "2CRV4P", "n_rows", "columns"]
old_docs = {  # earlier snapshots of the same documents (tic4206066_PUBLICATION_AUDIT/data/refs)
    "docs/qlp/hlsp_qlp_tess_ffi_all_tess_v1_readme.txt": "qlp_readme.txt",
    "docs/qlp/hlsp_qlp_tess_ffi_all_tess_v1_data-prod-desc.pdf": "qlp_data_prod_desc.pdf",
    "docs/qlp/mast_hlsp_qlp.html": "qlp_hlsp_page.html",
    "docs/qlp/mit_qlp.html": "qlp_mit_page.html",
    "docs/heasarc/sector98_summary.html": "heasarc_sector98_summary.html",
    "docs/heasarc/sector.html": "heasarc_tess_sector_table.html",
}
entries = []
for dp, dn, fn in os.walk(ROOT):
    for f in sorted(fn):
        p = os.path.join(dp, f); rel = os.path.relpath(p, ROOT)
        if rel in ("DATA_MANIFEST.json", "README.md", "REPORT.md") or f.endswith(".part") or f == ".DS_Store":
            continue
        e = dict(relpath=rel, bytes=os.path.getsize(p), sha256=sha256_file(p))
        pe = prov.get(rel)
        kind = rel.split("/")[0]
        if pe:
            e["url"] = pe.get("url"); e["utc"] = pe.get("utc"); e["utc_start"] = pe.get("utc_start")
            if pe["sha256"] != e["sha256"]:
                e["WARNING"] = "sha256 differs from value logged at download time"
            for k in ("http_headers", "query", "nrows", "ncols", "tool", "vizier_table", "zip_sha256", "zip_bytes", "note", "mast_listed_size"):
                if pe.get(k) is not None:
                    e[k] = pe[k]
            if "stpubdata.s3" in str(pe.get("url")):
                e["requested_url"] = "https://mast.stsci.edu/api/v0.1/Download/file?uri=" + str(pe.get("note", "")).replace("MAST product ", "")
                e["url_note"] = "MAST Download API redirected to the STScI public S3 bucket (final URL recorded in 'url')"
        if kind == "tesscut" and rel.endswith(".fits"):
            size = rel.split("/")[1]; sec = int(H[rel]["SECTOR"])
            e["url"] = tc_url.get((size, sec), e.get("url"))
        if rel.endswith(".fits") and rel in H:
            e["header"] = {k: H[rel][k] for k in KEYH if k in H[rel]}
            hh = H[rel]
            e["header"]["astrocut_version"] = hh.get("PROCVER") if hh.get("CREATOR") == "astrocut" else None
        # matches_existing
        if kind == "mast" and rel.endswith(".fits"):
            e["matches_existing"] = C["mast"][os.path.basename(rel)]["old"]
        elif kind == "tesscut" and rel.endswith(".fits"):
            t = C["tesscut"][rel]
            if "21x21" in rel:
                e["matches_existing"] = {k: ("arrays identical (TIME,TIMECORR,CADENCENO,RAW_CNTS,FLUX,FLUX_ERR,FLUX_BKG,FLUX_BKG_ERR,QUALITY,POS_CORR1/2)"
                                             if all(v == "identical" for v in o["arrays"].values()) else o["arrays"]) for k, o in t["old"].items()}
            else:
                e["matches_existing"] = {"no earlier 31x31 copy": None,
                                         "central 21x21 vs fresh 21x21": "identical" if all(v == "identical" for v in t["inner21_vs_fresh21"].values()) else t["inner21_vs_fresh21"]}
        elif kind == "prf" and rel.endswith(".fits"):
            key = "/".join(rel.split("/")[-2:])
            e["matches_existing"] = {"tic4206066_PUBLICATION_AUDIT/data/prf/" + key: C["prf"].get(key)}
        elif kind == "catalogs":
            if "gaia_source_cone240" in rel:
                g = C["gaia"]["gaia_source_target"]; nb = C["gaia"]["neighbours_180as"]
                e["matches_existing"] = {"tic4206066_PUBLICATION_AUDIT/data/gaia/gaia_dr3_source.csv (target row)": f"{g['n_same']} common columns equal ({g['n_same_within_rel1e6_only']} only to float32 print precision), {g['n_diff']} different",
                                         "tic4206066_PUBLICATION_AUDIT/data/gaia/gaia_dr3_neighbours_180as.csv": f"same {nb['old_n']} source_ids within 180 arcsec" if not nb["only_fresh"] and not nb["only_old"] else nb}
            elif "astrophysical_parameters" in rel:
                g = C["gaia"]["gaiadr3_astrophysical_parameters_target"]
                e["matches_existing"] = {"tic4206066_PUBLICATION_AUDIT/data/gaia/gaia_dr3_ap.csv": f"{g['n_same']} columns equal ({g['n_same_within_rel1e6_only']} only to float32 print precision), {g['n_diff']} different"}
            elif re.search(r"nss_|vari_summary", rel):
                key = next(k for k in C["gaia"] if k in rel)
                e["matches_existing"] = {"earlier gaia/ copy": f"both empty: {C['gaia'][key]}"}
            elif "tic82_cone300as" in rel:
                t = C["tic"]
                e["matches_existing"] = {"tic4206066_PUBLICATION_AUDIT/data/catalogs/tic82_cone5arcmin.csv": f"target row: {t['target_row']['n_same']} columns equal, differs only in {list(t['target_row']['diffs'])} (query-centre dependent); cone: fresh {t['cone']['fresh_n']} vs old {t['cone']['old_n']} rows, only_old={t['cone']['only_old']} (at 300.09 arcsec, outside the 300 arcsec radius used here)"}
            else:
                e["matches_existing"] = "not compared (earlier copies are VizieR TSV exports with a different query radius)"
        elif rel in old_docs:
            op = os.path.join(A, "tic4206066_PUBLICATION_AUDIT/data/refs", old_docs[rel])
            e["matches_existing"] = {os.path.relpath(op, A): ("identical" if os.path.exists(op) and sha256_file(op) == e["sha256"] else ("DIFFERENT (web page; dynamic content expected)" if os.path.exists(op) else "absent"))}
        if kind in ("scripts", "logs") or rel.endswith((".txt", ".csv", ".json", ".ecsv")) and not pe:
            e.setdefault("generated_by", "this provisioning step (00_data/scripts)")
        # normalise service endpoints for catalogue queries
        u = str(e.get("url") or "")
        if u.startswith("vizier.cds.unistra.fr"):
            e["url"] = "https://vizier.cds.unistra.fr/viz-bin/votable" + u[len("vizier.cds.unistra.fr"):]
        if u.startswith("https://gea.esac.esa.int/tap-server/tap ("):
            e["url"] = "https://gea.esac.esa.int/tap-server/tap/sync (astroquery.gaia Gaia.launch_job, synchronous ADQL)"
        if rel.startswith("mast/mast_"):
            e["url"] = "https://mast.stsci.edu/api/v0/invoke (astroquery.mast.Observations.query_region r=10 arcsec + query_criteria(target_name=[TIC/Gaia ids]) + get_product_list)"
            e["utc"] = json.load(open(os.path.join(ROOT, "mast", "mast_query_summary.json")))["utc"]
        if rel.startswith("docs/") and rel.endswith(".txt") and os.path.exists(p[:-4]) and not pe:
            e["derived_from"] = rel[:-4]; e["generated_by"] = "text extraction (pdftotext -layout / HTML tag strip) by scripts/g_docs.py"
        entries.append(e)
man = dict(created_utc=utcnow(), root=ROOT, target=dict(tic=TIC, gaia_dr3=GAIA_DR3, ra=RA, dec=DEC),
           n_files=len(entries), files=entries)
json.dump(man, open(os.path.join(ROOT, "DATA_MANIFEST.json"), "w"), indent=1, default=str)
print("entries", len(entries), "warnings", sum("WARNING" in e for e in entries), "no_url", [e["relpath"] for e in entries if not e.get("url") and "generated_by" not in e][:40])
