"""T8 check 7: S98 per-event values of the two events removed by the v3 gap rule (epochs -1 and 8), fitted with the
gap and dump rules switched off, and the 8/9/10-event means. Also the per-aperture corrected depths (S98, S32 TPF).
Writes only to OUT/_work/T8_pixel_recheck/out/."""
import json, os
import numpy as np
W = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(W, "t8_a_events_and_nulls.py")).read().split("res = {\"note\"")[0]
src = src.replace('if np.max(np.diff(t[w])) > 0.1: return "gap_gt_2.4h_in_window", None',
                  'if GAP_RULE and np.max(np.diff(t[w])) > 0.1: return "gap_gt_2.4h_in_window", None')
ns = {"__file__": os.path.join(W, "t8_a_events_and_nulls.py"), "GAP_RULE": False}
exec(compile(src, "t8_a_lib_nogap", "exec"), ns)
D = ns["prep"]("S98_ffi"); P, T0, T14 = ns["EPH"]["A"]; B = ns["EPH"]["B"]
rows = {}
for n in (-1, 0, 1, 4, 5, 8, 9, 10, 13, 14):
    reason, row = ns["fit_event"](D, T0 + n * P, T14, 0.25, [B], dump_rule=False)
    rows[n] = row["diff_1.0-2.0"] if row else reason
v8 = np.array([rows[n] for n in (0, 1, 4, 5, 9, 10, 13, 14)])
out = {"per_event_diff_1.0-2.0_no_gap_no_dump_rule": {str(k): v for k, v in rows.items()}}
for lab, eps in (("v3_8", (0, 1, 4, 5, 9, 10, 13, 14)), ("plus_epoch_-1", (-1, 0, 1, 4, 5, 9, 10, 13, 14)),
                 ("plus_epoch_8", (0, 1, 4, 5, 8, 9, 10, 13, 14)), ("plus_both", (-1, 0, 1, 4, 5, 8, 9, 10, 13, 14))):
    v = np.array([rows[n] for n in eps if not isinstance(rows[n], str)])
    out[lab] = {"n": len(v), "mean": float(v.mean()), "err_scatter": float(v.std(ddof=1) / np.sqrt(len(v))), "median": float(np.median(v))}
PKG = json.load(open(os.path.join(W, "pkg_copy/out/after_decorrelated.json")))
out["corrected_depths_by_radius"] = {s: {r: [round(PKG["A"][s]["summary"]["corrected"][r]["depth_ppm"], 1),
                                             round(PKG["A"][s]["summary"]["corrected"][r]["err_event_scatter"], 1)] for r in ("0.5", "1.0", "2.0", "3.0")}
                                     for s in ("S98_ffi", "S32_tpf2min", "S32_ffi")}
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(W, "out/t8_s98_rule_sensitivity.json"), "w"), indent=1)
