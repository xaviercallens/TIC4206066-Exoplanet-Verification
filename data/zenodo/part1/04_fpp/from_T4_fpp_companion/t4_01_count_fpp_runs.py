"""Count TRICERATOPS run files and recompute per-variant FPP statistics from stored run outputs (read-only)."""
import json, glob, os, csv
import numpy as np
ROOT = "/Users/paulsmith/space <3/tess-jev"
A = ROOT + "/addenda"
U = A + "/opus55_final_audit_tic4206066"
rows = []
def stats(name, src, sig, fpps, nfpps, extra=""):
    f = np.array(fpps, float); n = np.array([x for x in nfpps if x is not None], float)
    rows.append(dict(source=src, variant=name, signal=sig, n_runs=len(f), FPP_mean=f.mean(), FPP_sd=f.std(ddof=1) if len(f) > 1 else np.nan,
                     FPP_sd_ddof0=f.std(ddof=0), FPP_min=f.min(), FPP_max=f.max(),
                     NFPP_max=(n.max() if len(n) else np.nan), note=extra))
# U audit 08 runs
for f in sorted(glob.glob(U + "/08_fpp/runs/*/runs.json")):
    d = json.load(open(f)); r = d["runs"]
    a = d["args"]
    stats(d["variant"], os.path.relpath(f, ROOT), d["signal"], [x["FPP"] for x in r], [x.get("NFPP") for x in r],
          f"N={a['N']} fold={d['fold']} exptime={a['exptime_min']} Teff={a['Teff']} molusc={a['molusc']} stored_mean={d['FPP_mean']:.5f} stored_sd={d['FPP_sd']:.5f}")
# U verify runs
for f in sorted(glob.glob(U + "/08_fpp/verify/runs_v/*.json")):
    d = json.load(open(f))
    if "runs" not in d: continue
    r = d["runs"]
    if not r or "FPP" not in r[0]:
        rows.append(dict(source=os.path.relpath(f, ROOT), variant=os.path.basename(f), signal=d.get("signal"), n_runs=len(r), note="no per-run FPP key")); continue
    stats(os.path.basename(f)[:-5], os.path.relpath(f, ROOT), d.get("signal"), [x["FPP"] for x in r], [x.get("NFPP") for x in r],
          f"N={d.get('N')} mode={d.get('mode')} exptime={d.get('exptime_min')} Teff={d.get('Teff')}")
# U resolve_R5 runs
for f in sorted(glob.glob(U + "/referee/resolve_R5/runs/*.json")):
    d = json.load(open(f)); r = d["runs"]
    stats(os.path.basename(f)[:-5], os.path.relpath(f, ROOT), d.get("signal"), [x["FPP"] for x in r], [x.get("NFPP") for x in r],
          f"N={d.get('N')} mode={d.get('mode')} exptime={d.get('exptime_min')} star={d.get('star')}")
# U resolve_R-DIL runs
for f in sorted(glob.glob(U + "/referee/resolve_R-DIL/runs/*.json")):
    d = json.load(open(f))
    rows.append(dict(source=os.path.relpath(f, ROOT), variant=os.path.basename(f), signal=None, n_runs=None, note="flag-flip paired runs; keys=" + ",".join(list(d.keys())[:8])))
# PUBLICATION_AUDIT runs
for dd in sorted(glob.glob(A + "/tic4206066_PUBLICATION_AUDIT/out/fpp/runs/*")):
    nfiles = len(glob.glob(dd + "/probs_run*.csv"))
    s = json.load(open(dd + "/summary.json")) if os.path.exists(dd + "/summary.json") else {}
    fp = None
    for k in ("runs", "FPP_runs", "fpps", "FPP_list"):
        if k in s: fp = s[k]; break
    if isinstance(fp, list) and fp and isinstance(fp[0], dict):
        stats(os.path.basename(dd), os.path.relpath(dd, ROOT), os.path.basename(dd)[0], [x["FPP"] for x in fp], [x.get("NFPP") for x in fp], f"probs_run files={nfiles}")
    elif isinstance(fp, list):
        stats(os.path.basename(dd), os.path.relpath(dd, ROOT), os.path.basename(dd)[0], fp, [], f"probs_run files={nfiles}")
    else:
        rows.append(dict(source=os.path.relpath(dd, ROOT), variant=os.path.basename(dd), signal=os.path.basename(dd)[0], n_runs=nfiles,
                         note="summary keys=" + ",".join(list(s.keys())[:20]) + " | " + json.dumps({k: s[k] for k in s if 'FPP' in k.upper()})[:300]))
# pubprep own fold
for sig in "AB":
    f = A + f"/rnaas_tic4206066_pubprep/analysis/out/fpp/fpp_ownfold_{sig}.json"
    d = json.load(open(f))
    stats(f"pubprep_ownfold_{sig}", os.path.relpath(f, ROOT), sig, [x["FPP"] for x in d["runs"]], [x.get("NFPP") for x in d["runs"]], f"N={d['N']} star={d['star']}")
keys = ["source", "variant", "signal", "n_runs", "FPP_mean", "FPP_sd", "FPP_sd_ddof0", "FPP_min", "FPP_max", "NFPP_max", "note"]
out = os.path.dirname(os.path.abspath(__file__))
with open(out + "/t4_01_fpp_run_census.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, keys); w.writeheader()
    for r in rows: w.writerow({k: r.get(k) for k in keys})
for r in rows:
    m = r.get("FPP_mean")
    print(f"{r['variant'][:38]:38s} {str(r.get('signal')):4s} n={r.get('n_runs')} " + (f"{m:.4f}+-{r['FPP_sd']:.4f} [{r['FPP_min']:.4f}-{r['FPP_max']:.4f}] NFPPmax={r['NFPP_max']:.1e}" if m is not None else "") + " | " + str(r.get('note'))[:160])
