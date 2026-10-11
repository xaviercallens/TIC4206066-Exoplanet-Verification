"""T8 check 5: (i) Signal B S32 2-min per-event aperture differences with white errors of the difference light curve
(the package's diff_2.0-3.0 error uses the scatter of N = 2 events); (ii) first-pass values; (iii) paired comparison of
T9-method and per-event-method S98 values on the same 8 clean events. Writes only to OUT/_work/T8_pixel_recheck/out/."""
import json, os, sys
import numpy as np
W = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, W)
import importlib.util
spec = importlib.util.spec_from_file_location("ev", os.path.join(W, "t8_a_events_and_nulls.py"))
src = open(os.path.join(W, "t8_a_events_and_nulls.py")).read().split("res = {\"note\"")[0]
ns = {"__file__": os.path.join(W, "t8_a_events_and_nulls.py")}; exec(compile(src, "t8_a_lib", "exec"), ns)
prep, fit_event, EPH = ns["prep"], ns["fit_event"], ns["EPH"]
PK = json.load(open(os.path.join(W, "pkg_copy/out/before_global_detrend.json")))
res = {}
D = prep("S32_tpf2min"); P, T0, T14 = EPH["B"]; t = D["t"]
rows = []
for n in np.unique(np.round((t - T0) / P)):
    reason, row = fit_event(D, T0 + n * P, T14, 0.25, [EPH["A"]])
    if row: row["epoch"] = int(n); rows.append(row)
res["B_S32_tpf2min_per_event"] = rows
for k in ("diff_1.0-3.0", "diff_2.0-3.0", "diff_1.0-2.0", "diff_0.5-2.0"):
    v = np.array([r[k] for r in rows]); e = np.array([r[k + "_white"] for r in rows]); w = 1 / e**2
    res[f"B_{k}_combined_white"] = {"ppm": float((v * w).sum() / w.sum()), "err_white": float(1 / np.sqrt(w.sum())),
                                    "per_event": [[round(a, 1), round(b, 1)] for a, b in zip(v, e)]}
res["first_pass_S32_tpf2min_B"] = {r: PK["S32_tpf2min"]["apertures"][r]["B"] for r in PK["S32_tpf2min"]["apertures"]}
res["first_pass_S32_tpf2min_B_diffs"] = {k: v["B"] for k, v in PK["S32_tpf2min"]["diffs"].items()}
# (iii) paired S98 comparison on identical events
T9 = json.load(open(os.path.join(W, "out/t8_t9_decomposition.json")))["S98_clean_T9eph"]["restricted_to_clean8_events"]["per_event"]
A = json.load(open(os.path.join(W, "out/t8_a_events_and_nulls.json")))["S98_ffi"]["events"]
pe = {r["epoch"]: r["diff_1.0-2.0"] for r in A}
eps = sorted(pe)
d = np.array([T9[str(e)]["extra_ppm"] - pe[e] for e in eps])
res["S98_paired_T9method_minus_perevent_same8"] = {"epochs": eps, "per_event_diff_ppm": [round(x, 1) for x in d],
    "mean": float(d.mean()), "err_scatter": float(d.std(ddof=1) / np.sqrt(len(d)))}
json.dump(res, open(os.path.join(W, "out/t8_b_and_paired.json"), "w"), indent=1)
for k, v in res.items():
    if k != "B_S32_tpf2min_per_event": print(k, v)
for r in rows: print("B ep", r["epoch"], {k: round(v) for k, v in r["depth"].items()}, {k: round(r[k]) for k in r if k.startswith("diff")})
