"""T8 check 1: per-event (s03-style) aperture differences for Signal A with explicit drop reasons, and
error models for the S98 per-event statistic (event scatter, white noise, random-phase null).
Read-only on project data. Imports the package COPY's common.py (absolute ROOT) and writes only into
OUT/_work/T8_pixel_recheck/out/."""
import json, os, sys
import numpy as np
W = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(W, "pkg_copy", "scripts"))
from common import ROOT, EPH, load, radii_grid, centroid   # package copy, ROOT absolute
ODIR = os.path.join(W, "out"); os.makedirs(ODIR, exist_ok=True)
T7 = json.load(open(os.path.join(ROOT, "addenda/independent_verification_tic4206066_SEALED/out/t7_momentum_dumps.json")))
DUMPS = {"S6_ffi": T7["S06"]["dump_times"], "S32_ffi": T7["S32"]["dump_times"], "S32_tpf2min": T7["S32"]["dump_times"],
         "S98_ffi": T7["S98_tesscut"]["dump_times"]}
RAD = (0.5, 1.0, 2.0, 3.0)


def prep(name):
    t, F, bk, c, med = load(name)
    rr = radii_grid(F.shape[1:], c); cx, cy = centroid(F, c)
    ap = {r: F[:, rr <= r].sum(1) for r in RAD}
    return dict(t=t, ap=ap, cx=cx, cy=cy, bk=bk, name=name)


def fit_event(D, tc, T14, in_frac, masks, dump_rule=True, min_in=3):
    """Exactly the s03 selection and fit, for an arbitrary centre tc. masks: list of (P, T0, T14) whose
    +-0.75 T14 windows are removed from the out-of-transit set. Returns (reason, row)."""
    t = D["t"]; w = np.abs(t - tc) < 0.35
    if w.sum() < 20: return "lt20_points_in_window", None
    dt = t[w] - tc; inn = np.abs(dt) < in_frac * T14; oot = np.abs(dt) > 0.75 * T14
    for (Po, T0o, T14o) in masks:
        no = np.round((t[w] - T0o) / Po); oot &= np.abs(t[w] - (T0o + no * Po)) > 0.75 * T14o
    cad = np.median(np.diff(t)); expect = 2 * in_frac * T14 / cad
    if inn.sum() < max(min_in, 0.7 * expect): return "in_transit_coverage_lt70pct", None
    if np.max(np.diff(t[w])) > 0.1: return "gap_gt_2.4h_in_window", None
    if dump_rule and any(abs(x - tc) < 0.35 for x in DUMPS[D["name"]]): return "momentum_dump_in_window", None
    if inn.sum() < min_in or oot.sum() < 20 or (dt[oot] < 0).sum() < 5 or (dt[oot] > 0).sum() < 5: return "oot_insufficient", None
    x = D["cx"][w] - D["cx"][w][oot].mean(); y = D["cy"][w] - D["cy"][w][oot].mean(); b = D["bk"][w] - D["bk"][w][oot].mean()
    X = np.vstack([np.ones_like(dt), dt, dt**2, x, y, x*x, y*y, x*y, b]).T
    res = {}
    for r in RAD:
        fo = D["ap"][r][w] / np.median(D["ap"][r][w][oot])
        coef, *_ = np.linalg.lstsq(X[oot], fo[oot], rcond=None)
        res[r] = fo / (X @ coef)
    row = {"tc": float(tc), "n_in": int(inn.sum()), "depth": {str(r): float(1e6 * (1 - res[r][inn].mean())) for r in RAD}}
    for a, bb in ((1.0, 2.0), (0.5, 2.0), (1.0, 3.0), (2.0, 3.0)):
        d = res[a] - res[bb]
        row[f"diff_{a}-{bb}"] = float(1e6 * (-(d[inn].mean())))
        row[f"diff_{a}-{bb}_white"] = float(1e6 * np.std(d[oot], ddof=1) / np.sqrt(inn.sum()))
    return "used", row


def summarize(rows, key="diff_1.0-2.0"):
    v = np.array([r[key] for r in rows]); e = np.array([r[key + "_white"] for r in rows]); n = len(v)
    if n == 0: return {"n": 0}
    out = {"n": n, "mean": float(v.mean()), "median": float(np.median(v)),
           "err_event_scatter": float(v.std(ddof=1) / np.sqrt(n)) if n > 1 else None,
           "err_white_for_mean": float(np.sqrt((e**2).sum()) / n),
           "chi2_about_mean_white": float((((v - v.mean()) / e)**2).sum()), "dof": n - 1}
    return out


res = {"note": "s03-equivalent per-event fit re-implemented to record drop reasons and difference-light-curve white errors; "
               "values for used events must equal pkg out/after_decorrelated.json"}
A = EPH["A"]; B = EPH["B"]
for name, in_frac, min_in in (("S98_ffi", 0.25, 3), ("S32_tpf2min", 0.25, 3), ("S32_ffi", 0.25, 3), ("S6_ffi", 0.45, 2)):
    D = prep(name); P, T0, T14 = A; t = D["t"]
    book = []
    rows = []
    for n in np.unique(np.round((t - T0) / P)):
        tc = T0 + n * P
        reason, row = fit_event(D, tc, T14, in_frac, [B], min_in=min_in)
        # what would the rule-free fit give (for dropped events)?
        _, row_nr = fit_event(D, tc, T14, in_frac, [B], dump_rule=False, min_in=min_in)
        near_dump_h = float(min(abs(tc - x) for x in DUMPS[name]) * 24)
        w = np.abs(t - tc) < 0.35
        maxgap_h = float(np.max(np.diff(t[w])) * 24) if w.sum() > 1 else None
        book.append({"epoch": int(n), "tc": float(tc), "reason": reason, "nearest_dump_h": near_dump_h, "max_gap_in_window_h": maxgap_h,
                     "n_in_window": int(w.sum()),
                     "diff_1.0-2.0_if_fitted_without_dump_rule": (row_nr or {}).get("diff_1.0-2.0")})
        if row: row["epoch"] = int(n); rows.append(row)
    res[name] = {"in_frac": in_frac, "bookkeeping": book, "events": rows,
                 "summary_diff_1.0-2.0": summarize(rows), "summary_diff_2.0-3.0": summarize(rows, "diff_2.0-3.0"),
                 "summary_diff_1.0-3.0": summarize(rows, "diff_1.0-3.0")}
    # same fit with the T9-like in-transit window (full T14: |dt| < T14/2)
    if name != "S6_ffi":
        rows2 = []
        for n in np.unique(np.round((t - T0) / P)):
            reason, row = fit_event(D, T0 + n * P, T14, 0.5, [B])
            if row: row["epoch"] = int(n); rows2.append(row)
        res[name]["fullT14_window"] = {"epochs": [r["epoch"] for r in rows2], "summary_diff_1.0-2.0": summarize(rows2)}
    print(name, res[name]["summary_diff_1.0-2.0"], [(b["epoch"], b["reason"]) for b in book if b["reason"] != "used"])

# Random-phase null for the per-event statistic in S98 and S32 TPF (same rules; fake centres at T0 + n P + phi P,
# phi uniform in [0.15, 0.85] like T9; both A and B windows masked from the OOT set and fake in-transit windows
# overlapping a real A or B transit (+-0.75 T14 + 0.25 T14) rejected).
rng = np.random.default_rng(20260925)
for name in ("S98_ffi", "S32_tpf2min"):
    D = prep(name); P, T0, T14 = A; t = D["t"]
    stats, nev, singles = [], [], []
    for k in range(300):
        phi = rng.uniform(0.15, 0.85); rows = []
        for n in np.unique(np.round((t - T0) / P)):
            tc = T0 + (n + phi) * P
            bad = False
            for (Pp, T0p, T14p) in (A, B):
                nn = np.round((tc - T0p) / Pp)
                if abs(tc - (T0p + nn * Pp)) < 0.75 * T14p + 0.25 * T14: bad = True
            if bad: continue
            reason, row = fit_event(D, tc, T14, 0.25, [A, B])
            if row: rows.append(row["diff_1.0-2.0"])
        if len(rows) >= 3:
            stats.append(float(np.mean(rows))); nev.append(len(rows)); singles += rows
    stats = np.array(stats); singles = np.array(singles)
    real = res[name]["summary_diff_1.0-2.0"]
    nreal = real["n"]
    res[name]["random_phase_null"] = {
        "n_trials": int(len(stats)), "median_events_per_trial": float(np.median(nev)),
        "std_of_trial_means": float(stats.std()), "mean_of_trial_means": float(stats.mean()),
        "single_event_std": float(singles.std()), "n_single": int(len(singles)),
        "single_event_std_over_sqrtN_real": float(singles.std() / np.sqrt(nreal)),
        "frac_trials_ge_real_mean": float((stats >= real["mean"]).mean())}
    print(name, "null:", res[name]["random_phase_null"])

json.dump(res, open(os.path.join(ODIR, "t8_a_events_and_nulls.json"), "w"), indent=1)
