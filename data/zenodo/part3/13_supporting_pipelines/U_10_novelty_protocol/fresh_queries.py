"""Audit 10 (novelty): fresh, independent spot-check queries for prior public reports of the
TIC 4206066 transit-like signals.  Written from scratch for this audit (no code imported/copied
from earlier analyses; earlier derived files are NOT used as inputs).

Every network query is logged to queries.log (UTC time, source, URL, HTTP status, bytes, sha256,
one-line result) and the raw response is stored in raw/.  Summary -> fresh_results.json.

Run:  AUDIT_PY fresh_queries.py [section ...]   (sections: exofop tce nea simbad vizier arxiv t16raven)
"""
import csv, io, json, math, os, sys, time, hashlib, datetime, urllib.parse
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw"); os.makedirs(RAW, exist_ok=True)
LOG = os.path.join(HERE, "queries.log")
RES = os.path.join(HERE, "fresh_results.json")
UA = {"User-Agent": "novelty-audit/1.0 (python-requests; independent literature/catalogue check)"}

TIC = 4206066
GAIA = "3220388198192519424"
RA0, DEC0 = 80.85982665567, -1.32696401034          # Gaia DR3 (epoch 2016.0)
# periods: manuscript values and the Sector-98 search values
PERIODS = {"A_ms": 3.1827856, "A_s98": 3.182945192, "B_ms": 11.13274, "B_s98": 11.132975672}
RELS = {"P/3": 1 / 3, "P/2": 0.5, "P": 1.0, "2P": 2.0, "3P": 3.0}
REL_TOL = 1e-3            # 0.1 %
NEAR_ARCSEC = 120.0       # 2 arcmin


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def log(source, url, status, nbytes, sha, result):
    line = "\t".join([now(), source, url, str(status), str(nbytes), sha, result.replace("\n", " ")])
    with open(LOG, "a") as f:
        f.write(line + "\n")
    print(line[:400], flush=True)


def fetch(source, url, fname, params=None, timeout=300, method="get", data=None):
    t0 = now()
    try:
        if method == "get":
            r = requests.get(url, params=params, headers=UA, timeout=timeout)
        else:
            r = requests.post(url, data=data, headers=UA, timeout=timeout)
        b, st, full = r.content, r.status_code, r.url
    except Exception as e:  # network failure is logged, never hidden
        b, st, full = repr(e).encode(), -1, url + ("?" + urllib.parse.urlencode(params) if params else "")
    sha = hashlib.sha256(b).hexdigest()
    with open(os.path.join(RAW, fname), "wb") as f:
        f.write(b)
    with open(os.path.join(RAW, "provenance.jsonl"), "a") as f:
        f.write(json.dumps(dict(request_utc=t0, source=source, url=full, status=st, bytes=len(b),
                                sha256=sha, file="raw/" + fname)) + "\n")
    return b, st, full, sha


def fetch_curl(source, url, fname, timeout=120):
    """arXiv's export API answered HTTP 406 to python-requests (tested 2026-09-25) but 200 to curl; use curl there."""
    import subprocess
    t0 = now()
    fp = os.path.join(RAW, fname)
    try:
        cp = subprocess.run(["curl", "-s", "-L", "-A", "novelty-audit/1.0", "-o", fp, "-w", "%{http_code}", "--max-time", str(timeout), url],
                            capture_output=True, text=True)
        st = int(cp.stdout.strip() or -1)
        b = open(fp, "rb").read() if os.path.exists(fp) else b""
    except Exception as e:
        b, st = repr(e).encode(), -1
        open(fp, "wb").write(b)
    sha = hashlib.sha256(b).hexdigest()
    with open(os.path.join(RAW, "provenance.jsonl"), "a") as f:
        f.write(json.dumps(dict(request_utc=t0, source=source, url=url, status=st, bytes=len(b), sha256=sha, file="raw/" + fname, via="curl")) + "\n")
    return b, st, url, sha


def sep_arcsec(ra, dec, ra0=RA0, dec0=DEC0):
    r = math.radians
    d = math.sin((r(dec) - r(dec0)) / 2) ** 2 + math.cos(r(dec0)) * math.cos(r(dec)) * math.sin((r(ra) - r(ra0)) / 2) ** 2
    return math.degrees(2 * math.asin(min(1.0, math.sqrt(d)))) * 3600.0


def sexa_to_deg(ra, dec):
    ra, dec = str(ra).strip(), str(dec).strip()
    if ":" in ra:
        h, m, s = [float(x) for x in ra.split(":")]
        rad = 15 * (h + m / 60 + s / 3600)
        sg = -1 if dec.startswith("-") else 1
        d, mm, ss = [abs(float(x)) for x in dec.split(":")]
        return rad, sg * (d + mm / 60 + ss / 3600)
    return float(ra), float(dec)


def period_matches(p):
    """return list of (signal, relation, rel_diff) for which p matches rel*P within REL_TOL"""
    out = []
    try:
        p = float(p)
    except Exception:
        return out
    if not (p > 0):
        return out
    for sname, P in PERIODS.items():
        for rn, rv in RELS.items():
            ref = P * rv
            d = abs(p - ref) / ref
            if d <= REL_TOL:
                out.append((sname, rn, round(d, 6)))
    return out


def load_results():
    return json.load(open(RES)) if os.path.exists(RES) else {}


def save_results(res):
    json.dump(res, open(RES, "w"), indent=1, default=str)


# ------------------------------------------------------------------------------------------------
def sec_exofop(res):
    out = {}
    for key, url, idcol in (("toi", "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv", "TOI"),
                            ("ctoi", "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv", "CTOI")):
        b, st, full, sha = fetch("ExoFOP " + key.upper() + " list", url, f"exofop_{key}.csv")
        rows = list(csv.DictReader(io.StringIO(b.decode("utf-8", "replace")))) if st == 200 else []
        tic_rows, near, pmatch_any, pmatch_near, nearest = [], [], 0, [], (1e99, None)
        for r in rows:
            try:
                tid = int(float(r["TIC ID"]))
            except Exception:
                tid = None
            if tid == TIC:
                tic_rows.append(r.get(idcol))
            try:
                ra, de = sexa_to_deg(r["RA"], r["Dec"])
                s = sep_arcsec(ra, de)
            except Exception:
                s = None
            per = r.get("Period (days)")
            pm = period_matches(per)
            if pm:
                pmatch_any += 1
            if s is not None:
                if s < nearest[0]:
                    nearest = (s, r.get(idcol), tid)
                if s <= NEAR_ARCSEC:
                    near.append({idcol: r.get(idcol), "tic": tid, "sep_arcsec": round(s, 2), "period": per})
                if pm and s <= 3600:
                    pmatch_near.append({idcol: r.get(idcol), "tic": tid, "sep_arcsec": round(s, 1), "period": per, "match": pm})
        # closest period-only matches (any sky position), for the record
        pm_list = []
        for r in rows:
            pm = period_matches(r.get("Period (days)"))
            if pm:
                try:
                    ra, de = sexa_to_deg(r["RA"], r["Dec"]); s = sep_arcsec(ra, de) / 3600
                except Exception:
                    s = None
                pm_list.append({idcol: r.get(idcol), "tic": r.get("TIC ID"), "period": r.get("Period (days)"), "sep_deg": None if s is None else round(s, 2), "match": pm})
        pm_list.sort(key=lambda x: (x["sep_deg"] is None, x["sep_deg"] if x["sep_deg"] is not None else 0))
        out[key] = dict(http=st, url=full, sha256=sha, n_rows=len(rows), tic_rows=tic_rows,
                        within_2arcmin=near, nearest={"sep_arcmin": round(nearest[0] / 60, 2), "id": nearest[1], "tic": nearest[2]},
                        n_period_only_matches_0p1pct_anywhere=pmatch_any, period_matches_within_1deg=pmatch_near,
                        closest_period_only_matches=pm_list[:5])
        log("ExoFOP " + key.upper() + " list", full, st, len(b), sha,
            f"rows={len(rows)}; TIC {TIC} rows={tic_rows}; within 2'={len(near)}; nearest {nearest[1]} at {nearest[0]/60:.1f}'; "
            f"period-only matches (0.1%, P/3..3P, any star)={pmatch_any}; period matches within 1 deg={len(pmatch_near)}")
        time.sleep(2)
    url = f"https://exofop.ipac.caltech.edu/tess/target.php?id={TIC}&json"
    b, st, full, sha = fetch("ExoFOP target page JSON", url, "exofop_target.json")
    try:
        j = json.loads(b)
        summ = {k: (len(v) if isinstance(v, list) else v) for k, v in j.items() if k in
                ("tois", "ctois", "planet_parameters", "imaging", "spectroscopy", "time_series", "files", "obsnotes", "tags", "stellar_companions")}
        summ["tois_content"], summ["ctois_content"] = j.get("tois"), j.get("ctois")
        names = j.get("star_names") or j.get("basic_info", {}).get("star_names")
        summ["star_names"] = names
        summ["top_keys"] = list(j.keys())
    except Exception as e:
        summ = {"parse_error": repr(e)}
    out["target_json"] = dict(http=st, url=full, sha256=sha, summary=summ)
    log("ExoFOP target page JSON", full, st, len(b), sha, "summary=" + json.dumps(summ)[:300])
    res["exofop"] = out


def sec_tce(res):
    out = {}
    url = f"https://exo.mast.stsci.edu/api/v0.1/dvdata/tess/{TIC}/tces/"
    b, st, full, sha = fetch("exo.MAST DV TCE service", url, "exomast_tces.json")
    out["exomast"] = dict(http=st, url=full, sha256=sha, response=b.decode("utf-8", "replace")[:300])
    log("exo.MAST DV TCE service", full, st, len(b), sha, "response=" + b.decode("utf-8", "replace")[:120])
    # TIC IDs within 2 arcmin (MAST TIC cone) for neighbour matching in TCE tables
    from astroquery.mast import Catalogs
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    t0 = now()
    try:
        tic_tab = Catalogs.query_region(SkyCoord(RA0, DEC0, unit="deg"), radius=2 * u.arcmin, catalog="TIC")
        neigh = {int(x): float(d) * 3600 for x, d in zip(tic_tab["ID"], tic_tab["dstArcSec"] if "dstArcSec" in tic_tab.colnames else tic_tab["distance"])}
        tic_tab.write(os.path.join(RAW, "mast_tic_cone_2arcmin.csv"), format="csv", overwrite=True)
        shab = hashlib.sha256(open(os.path.join(RAW, "mast_tic_cone_2arcmin.csv"), "rb").read()).hexdigest()
        log("MAST TIC 8.2 cone (astroquery Catalogs.query_region, catalog=TIC)", f"MAST Catalogs.Tic.Cone ra={RA0} dec={DEC0} r=2arcmin (queried {t0})", 200,
            os.path.getsize(os.path.join(RAW, "mast_tic_cone_2arcmin.csv")), shab, f"{len(neigh)} TIC objects within 2'")
    except Exception as e:
        neigh = {TIC: 0.0}
        log("MAST TIC 8.2 cone", "Catalogs.query_region", -1, 0, "-", "ERROR " + repr(e))
    out["n_tic_within_2arcmin"] = len(neigh)
    base = "https://archive.stsci.edu/missions/tess/catalogs/tce/"
    files = ["tess2018349182737-s0006-s0006_dvr-tcestats.csv",       # S6 single-sector (2-min)
             "tess2020325171311-s0032-s0032_dvr-tcestats.csv",       # S32 single-sector (2-min)
             "tess2025312234920-s0098-s0098_dvr-tcestats.csv",       # S98 single-sector (2-min)
             "tess2018206190142-s0001-s0036_dvr-tcestats.csv",       # multi-sector incl. S6, S32
             "tess2018206190142-s0001-s0039_dvr-tcestats.csv",
             "tess2018206190142-s0001-s0069_dvr-tcestats.csv",
             "tess2018206190142-s0001-s0096_dvr-tcestats.csv"]
    tabs = {}
    for fn in files:
        b, st, full, sha = fetch("MAST SPOC TCE bulk CSV", base + fn, "tce_" + fn)
        txt = b.decode("utf-8", "replace")
        lines = [l for l in txt.splitlines() if l and not l.startswith("#")]
        rows = list(csv.DictReader(lines)) if st == 200 else []
        tic_hits, neigh_hits, pm_neigh, n_pm = [], [], [], 0
        for r in rows:
            try:
                tid = int(float(r["ticid"]))
            except Exception:
                continue
            pm = period_matches(r.get("tce_period"))
            if pm:
                n_pm += 1
            if tid == TIC:
                tic_hits.append({"tceid": r.get("tceid"), "period": r.get("tce_period")})
            elif tid in neigh:
                neigh_hits.append({"tceid": r.get("tceid"), "tic": tid, "sep_arcsec": round(neigh[tid], 1), "period": r.get("tce_period"), "match": pm})
        tabs[fn] = dict(http=st, sha256=sha, n_rows=len(rows), target_rows=tic_hits, neighbour_rows_2arcmin=neigh_hits,
                        n_period_only_matches_other_stars=n_pm)
        log("MAST SPOC TCE bulk CSV", full, st, len(b), sha,
            f"rows={len(rows)}; ticid==4206066: {len(tic_hits)}; TIC within 2' (n={len(neigh)}): {len(neigh_hits)}; period-only matches (0.1%, other TICs)={n_pm}")
        time.sleep(1)
    out["bulk"] = tabs
    res["tce"] = out


def nea_tap(query, fname, res_list):
    url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
    b, st, full, sha = fetch("NASA Exoplanet Archive TAP", url, fname, params={"query": query, "format": "csv"})
    txt = b.decode("utf-8", "replace")
    rows = list(csv.DictReader(io.StringIO(txt))) if st == 200 else []
    log("NASA Exoplanet Archive TAP", full, st, len(b), sha, f"query=[{query}] -> {len(rows)} rows" + (" ; " + txt[:200] if st != 200 else ""))
    res_list.append(dict(query=query, http=st, sha256=sha, n_rows=len(rows), rows=rows[:20]))
    time.sleep(1)


def sec_nea(res):
    out = []
    dd = 2 / 60.0
    dra = dd / math.cos(math.radians(DEC0))
    box = f"ra between {RA0 - dra:.5f} and {RA0 + dra:.5f} and dec between {DEC0 - dd:.5f} and {DEC0 + dd:.5f}"
    nea_tap(f"select toi,tid,ra,dec,pl_orbper,tfopwg_disp from toi where tid={TIC}", "nea_toi_tid.csv", out)
    nea_tap(f"select toi,tid,ra,dec,pl_orbper,tfopwg_disp from toi where {box}", "nea_toi_box2arcmin.csv", out)
    nea_tap(f"select pl_name,hostname,tic_id,gaia_dr3_id,ra,dec,pl_orbper,disc_facility from ps where tic_id='TIC {TIC}'", "nea_ps_ticid.csv", out)
    nea_tap(f"select pl_name,hostname,tic_id,ra,dec,pl_orbper from ps where {box}", "nea_ps_box2arcmin.csv", out)
    nea_tap(f"select pl_name,hostname,tic_id,gaia_dr3_id,ra,dec,pl_orbper from pscomppars where tic_id='TIC {TIC}'", "nea_pscomppars_ticid.csv", out)
    nea_tap(f"select pl_name,hostname,tic_id,ra,dec,pl_orbper from pscomppars where {box}", "nea_pscomppars_box2arcmin.csv", out)
    nea_tap(f"select pl_name,hostname,gaia_dr3_id from pscomppars where gaia_dr3_id='Gaia DR3 {GAIA}'", "nea_pscomppars_gaia.csv", out)
    # positive control: a known TOI host must return rows (TOI-700 = TIC 150428135)
    nea_tap("select toi,tid,pl_orbper from toi where tid=150428135", "nea_toi_control_tic150428135.csv", out)
    res["nasa_exoplanet_archive"] = out


def simbad_tap(query, fname, out):
    url = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
    b, st, full, sha = fetch("SIMBAD TAP", url, fname, params={"request": "doQuery", "lang": "adql", "format": "csv", "query": query})
    txt = b.decode("utf-8", "replace")
    rows = list(csv.DictReader(io.StringIO(txt))) if st == 200 else []
    log("SIMBAD TAP", full, st, len(b), sha, f"query=[{query}] -> {len(rows)} rows" + ("; " + txt[:200] if st != 200 else ""))
    out.append(dict(query=query, http=st, sha256=sha, n_rows=len(rows), rows=rows[:60]))
    time.sleep(1)
    return rows


def sec_simbad(res):
    out = []
    for ident in ("TIC 4206066", "StKM 1-561", "Gaia DR3 3220388198192519424"):
        simbad_tap(f"SELECT b.main_id, b.otype, b.ra, b.dec, b.oid FROM ident AS i JOIN basic AS b ON i.oidref=b.oid WHERE i.id='{ident}'",
                   "simbad_ident_" + ident.replace(" ", "_") + ".csv", out)
    simbad_tap("SELECT i.id FROM ident AS i JOIN basic AS b ON i.oidref=b.oid WHERE b.main_id='StKM 1-561'", "simbad_all_ids.csv", out)
    simbad_tap("SELECT o.otype FROM otypes AS o JOIN basic AS b ON o.oidref=b.oid WHERE b.main_id='StKM 1-561'", "simbad_otypes.csv", out)
    simbad_tap("SELECT r.bibcode, r.\"year\", r.title FROM has_ref AS h JOIN ref AS r ON h.oidbibref=r.oidbib JOIN basic AS b ON h.oidref=b.oid "
               "WHERE b.main_id='StKM 1-561' ORDER BY 2", "simbad_refs.csv", out)
    simbad_tap(f"SELECT main_id, otype, ra, dec, DISTANCE(POINT('ICRS',ra,dec),POINT('ICRS',{RA0},{DEC0}))*3600 AS sep_arcsec FROM basic "
               f"WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{RA0},{DEC0},{30/3600:.6f}))=1 ORDER BY sep_arcsec", "simbad_cone30arcsec.csv", out)
    res["simbad"] = out


def sec_vizier(res):
    from astroquery.vizier import Vizier
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    t0 = now()
    v = Vizier(columns=["**", "_r"], row_limit=-1, timeout=600)
    try:
        tl = v.query_region(SkyCoord(RA0, DEC0, unit="deg"), radius=15 * u.arcsec)
    except Exception as e:
        log("VizieR all catalogues cone 15\"", "astroquery Vizier.query_region", -1, 0, "-", "ERROR " + repr(e))
        res["vizier"] = {"error": repr(e)}
        return
    keys = list(tl.keys())
    # store a compact dump of every table (name, description, columns, rows)
    dump = []
    kw = ("transit", "eclips", "planet", "variab", "period", "rotat", "binar", "toi", "tess", "flare", "candidate")
    flagged = []
    for k in keys:
        t = tl[k]
        desc = (t.meta.get("description") or "")
        cols = list(t.colnames)
        percols = [c for c in cols if c.lower().startswith("per") or c.lower() in ("p", "porb", "prot", "prot1", "period") or "period" in (t[c].description or "").lower()]
        entry = dict(table=k, description=desc, nrows=len(t), period_like_columns=percols)
        if percols:
            entry["period_values"] = [[str(x) for x in t[c][:5]] for c in percols]
        if any(w in (desc + " " + k).lower() for w in kw) or percols:
            flagged.append(entry)
        dump.append(entry)
    fn = os.path.join(RAW, "vizier_cone15_tables.json")
    json.dump(dump, open(fn, "w"), indent=1, default=str)
    sha = hashlib.sha256(open(fn, "rb").read()).hexdigest()
    log("VizieR all catalogues cone 15\" (astroquery Vizier.query_region, catalog=all, row_limit=-1)",
        f"https://vizier.cds.unistra.fr/viz-bin/votable?-c={RA0}%20{DEC0}&-c.rs=15 (queried {t0})", 200, os.path.getsize(fn), sha,
        f"{len(keys)} tables with rows; {len(flagged)} flagged by keyword/period column")
    res["vizier"] = dict(n_tables=len(keys), flagged=flagged, all_tables=[d["table"] for d in dump])


def sec_arxiv(res):
    base = "https://export.arxiv.org/api/query"
    qs = [('all:"TIC 4206066"', "id_tic"), ('all:"4206066"', "id_ticnum"), ('all:"StKM 1-561"', "id_stkm"),
          ('all:"3220388198192519424"', "id_gaia"), ('all:"TYC 4753-1591-1"', "id_tyc"), ('all:"J05232636-0119380"', "id_2mass"),
          ('all:"3.1828" AND all:TESS', "p_31828"), ('all:"3.183" AND all:TESS', "p_3183"), ('all:"3.18279" AND all:TESS', "p_318279"),
          ('all:"11.133" AND all:TESS', "p_11133"), ('all:"11.13" AND all:TESS', "p_1113"), ('all:"11.1327" AND all:TESS', "p_111327"),
          ('all:"5.566" AND all:TESS', "p_5566"), ('all:"22.27" AND all:TESS', "p_2227"), ('all:"22.265" AND all:TESS', "p_22265"),
          ('all:"1.591" AND all:TESS', "p_1591"), ('all:"6.366" AND all:TESS', "p_6366"), ('all:"1.061" AND all:TESS', "p_1061"),
          ('all:"9.548" AND all:TESS', "p_9548"), ('all:"3.711" AND all:TESS', "p_3711"), ('all:"33.40" AND all:TESS', "p_3340"),
          ('all:"TOI-2768"', "control_toi2768")]
    out = []
    import xml.etree.ElementTree as ET
    ns = {"a": "http://www.w3.org/2005/Atom", "o": "http://a9.com/-/spec/opensearch/1.1/"}
    for q, tag in qs:
        for attempt in range(4):
            b, st, full, sha = fetch_curl("arXiv API", base + "?search_query=" + urllib.parse.quote(q) + "&start=0&max_results=200",
                                          f"arxiv_{tag}.xml")
            if st == 200:
                break
            time.sleep(10 * (attempt + 1))
        entries, hits = [], []
        total = None
        if st == 200:
            try:
                root = ET.fromstring(b)
                tr = root.find("o:totalResults", ns)
                total = int(tr.text) if tr is not None else None
                for e in root.findall("a:entry", ns):
                    eid = e.find("a:id", ns).text
                    title = " ".join(e.find("a:title", ns).text.split())
                    summ = " ".join((e.find("a:summary", ns).text or "").split())
                    txt = (title + " " + summ)
                    target_hit = any(s in txt for s in ("4206066", "StKM 1-561", GAIA, "4753-1591", "J05232636"))
                    entries.append({"id": eid, "title": title[:150], "target_identifier_in_abstract": target_hit})
                    if target_hit:
                        hits.append(eid)
            except Exception as e:
                entries = [{"parse_error": repr(e)}]
        out.append(dict(query=q, http=st, sha256=sha, total_results=total, n_returned=len(entries), target_identifier_hits=hits, entries=entries))
        log("arXiv API (metadata: title/abstract/comments)", full, st, len(b), sha,
            f"query=[{q}] total={total} returned={len(entries)} entries-with-target-identifier={hits}")
        time.sleep(4)
    res["arxiv"] = out


def sec_t16raven(res):
    out = {}
    # T16 Planet Hunt (Roth et al. 2026) Dataverse dataset
    durl = "https://dataverse.harvard.edu/api/datasets/:persistentId/"
    b, st, full, sha = fetch("Harvard Dataverse T16 dataset metadata", durl, "t16_dataset.json", params={"persistentId": "doi:10.7910/DVN/CWEUGW"})
    log("Harvard Dataverse T16 dataset metadata", full, st, len(b), sha, "dataset metadata")
    files = []
    try:
        j = json.loads(b)
        for f in j["data"]["latestVersion"]["files"]:
            df = f["dataFile"]
            files.append((df["id"], df.get("filename"), df.get("filesize"), df.get("originalFileName")))
        out["t16_version"] = [j["data"]["latestVersion"].get("versionNumber"), j["data"]["latestVersion"].get("versionMinorNumber"),
                              j["data"]["latestVersion"].get("releaseTime")]
    except Exception as e:
        out["t16_meta_error"] = repr(e)
    out["t16_files"] = files
    t16 = []
    for fid, fnm, fsz, orig in files:
        if fnm and "candidate" in fnm.lower() and (fsz or 0) < 50e6:
            b, st, full, sha = fetch("Harvard Dataverse T16 datafile", f"https://dataverse.harvard.edu/api/access/datafile/{fid}", f"t16_{fid}.tab",
                                     params={"format": "original"})
            txt = b.decode("utf-8", "replace")
            delim = "\t" if txt.count("\t") > txt.count(",") else ","
            rows = list(csv.DictReader(io.StringIO(txt), delimiter=delim)) if st == 200 else []
            hit, near, nearest = [], [], (1e99, None)
            for r in rows:
                tid = str(r.get("TICID", r.get("TIC", ""))).strip()
                gid = str(r.get("GAIAID", "")).strip()
                if tid == str(TIC) or gid == GAIA:
                    hit.append(r)
                try:
                    s = sep_arcsec(float(r["ra"]), float(r["dec"]))
                    if s < nearest[0]:
                        nearest = (s, tid)
                    if s <= NEAR_ARCSEC:
                        near.append(tid)
                except Exception:
                    pass
            t16.append(dict(file=fnm, n_rows=len(rows), target_rows=len(hit), within_2arcmin=near, nearest_arcmin=round(nearest[0] / 60, 1)))
            log("Harvard Dataverse T16 datafile", full, st, len(b), sha,
                f"{fnm}: rows={len(rows)} target rows={len(hit)} within 2'={len(near)} nearest={nearest[0]/60:.1f}'")
    out["t16"] = t16
    # RAVEN (Zenodo 19661443)
    b, st, full, sha = fetch("Zenodo RAVEN record", "https://zenodo.org/api/records/19661443", "raven_record.json")
    log("Zenodo RAVEN record", full, st, len(b), sha, "record metadata")
    rav = []
    try:
        j = json.loads(b)
        out["raven_record"] = {"title": j.get("metadata", {}).get("title"), "publication_date": j.get("metadata", {}).get("publication_date"),
                               "version": j.get("metadata", {}).get("version")}
        for f in j.get("files", []):
            key, size = f.get("key"), f.get("size")
            link = f.get("links", {}).get("self")
            if key and key.endswith(".csv") and (size or 0) < 60e6:
                b2, st2, full2, sha2 = fetch("Zenodo RAVEN table", link, "raven_" + key)
                rows = list(csv.DictReader(io.StringIO(b2.decode("utf-8", "replace")))) if st2 == 200 else []
                hit, near, nearest = 0, [], (1e99, None)
                for r in rows:
                    if str(r.get("ticid", "")).strip() in (str(TIC), f"{TIC}.0") or str(r.get("gaia_dr3", "")).strip() == GAIA:
                        hit += 1
                    try:
                        s = sep_arcsec(float(r["ra"]), float(r["dec"]))
                        if s < nearest[0]:
                            nearest = (s, r.get("ticid"))
                        if s <= NEAR_ARCSEC:
                            near.append(r.get("ticid"))
                    except Exception:
                        pass
                rav.append(dict(file=key, n_rows=len(rows), target_rows=hit, within_2arcmin=near, nearest_deg=round(nearest[0] / 3600, 2)))
                log("Zenodo RAVEN table", full2, st2, len(b2), sha2, f"{key}: rows={len(rows)} target rows={hit} within 2'={len(near)} nearest={nearest[0]/3600:.2f} deg")
    except Exception as e:
        out["raven_error"] = repr(e)
    out["raven"] = rav
    res["t16_raven"] = out


SECTIONS = {"exofop": sec_exofop, "tce": sec_tce, "nea": sec_nea, "simbad": sec_simbad, "vizier": sec_vizier,
            "arxiv": sec_arxiv, "t16raven": sec_t16raven}

if __name__ == "__main__":
    if not os.path.exists(LOG):
        with open(LOG, "w") as f:
            f.write("# queries.log - Audit 10 (novelty) fresh spot-check queries, TIC 4206066\n"
                    "# columns: utc_logged<TAB>source<TAB>url<TAB>http_status<TAB>bytes<TAB>sha256(raw response)<TAB>result\n")
    todo = sys.argv[1:] or list(SECTIONS)
    for s in todo:
        res = load_results()
        res.setdefault("_runs", []).append({"section": s, "start_utc": now()})
        try:
            SECTIONS[s](res)
        except Exception as e:
            import traceback
            res.setdefault("_errors", {})[s] = traceback.format_exc()
            log("SECTION " + s, "-", -1, 0, "-", "ERROR " + repr(e))
        res["_runs"][-1]["end_utc"] = now()
        save_results(res)
