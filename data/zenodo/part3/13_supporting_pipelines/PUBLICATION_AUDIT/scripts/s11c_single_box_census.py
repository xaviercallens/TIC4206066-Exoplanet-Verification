"""Step 11c: census of single 1.6-h box dips in the archival light curves used by the B link test (s11, s11b).
Same preparation as s11 (detrending with only A masked, A transits removed). Box S/N at every position (step dur/4),
weights 1/sig^2 as in s11. Reported per data set: the S/N of the boxes at B's predicted archival times (linked
ephemeris), the rank of those boxes among all independent boxes, and the strongest other features with their times
(to compare with momentum dumps / sector edges). Identifies what drives the heavy tail of the S32 null distribution.
"""
import json
import os
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
import s11_nullB as s11  # noqa
from audit_common import save_json  # noqa

DUR = 1.6 / 24
EPH_CUR = json.load(open(os.path.join(lcio.AUD, "out", "ephem_current.json")))
PB, TB = EPH_CUR["B"]["P"], EPH_CUR["B"]["T0"]


def dumps(key):
    fn = lcio.PRODUCTS[key][0]
    d = fits.open(os.path.join(lcio.MAST, fn))[1].data
    q = np.asarray(d["QUALITY"]).astype(np.int64)
    t = np.asarray(d["TIME"], float)
    return t[(q & 32) != 0]                                  # bit 6 (value 32): reaction wheel desaturation (momentum dump)


def census(name, lc, dump_key=None):
    t, f, w = s11.prep(lc)
    centers = np.arange(t.min() + DUR / 2, t.max() - DUR / 2, DUR / 4)
    snr, n = np.full(len(centers), np.nan), np.zeros(len(centers), int)
    nmin = 2 if np.median(np.diff(t)) > 0.01 else 5            # 30-min S6 cadence: a 1.6-h box holds only ~3 points
    for i, c in enumerate(centers):
        m = np.abs(t - c) < DUR / 2
        if m.sum() >= nmin:
            sw = w[m].sum()
            snr[i] = (1 - np.sum(w[m] * f[m]) / sw) * np.sqrt(sw)
            n[i] = m.sum()
    ok = np.isfinite(snr)
    centers, snr = centers[ok], snr[ok]
    # B predicted times inside this data set
    ep = np.arange(np.floor((t.min() - TB) / PB), np.ceil((t.max() - TB) / PB) + 1)
    tb = TB + ep * PB
    tb = tb[(tb > t.min()) & (tb < t.max())]
    nearB = np.zeros(len(centers), bool)
    for x in tb:
        nearB |= np.abs(centers - x) < DUR
    # independent peaks: greedy non-maximum suppression (boxes > 1 duration apart), outside B
    order = np.argsort(-snr)
    peaks = []
    for i in order:
        if nearB[i]:
            continue
        if all(abs(centers[i] - p[0]) > DUR for p in peaks):
            peaks.append((float(centers[i]), float(snr[i])))
        if len(peaks) >= 8:
            break
    n_indep = int((t.max() - t.min()) / DUR)
    bvals = []
    for x in tb:
        m = np.abs(centers - x) < DUR / 2
        if m.any():
            j = np.flatnonzero(m)[np.argmax(snr[m])]
            rank = int(np.sum(np.array([p[1] for p in peaks]) > snr[j])) + 1
            bvals.append(dict(t_pred=float(x), best_box_t=float(centers[j]), snr=float(snr[j]), rank_among_non_B_peaks=rank))
    dd = dumps(dump_key) if dump_key else np.array([])
    if len(dd) == 0:
        dd = np.array([])
    for k, p in enumerate(peaks):
        near = dd[np.abs(dd - p[0]) < 0.5]
        peaks[k] = dict(t=p[0], snr=p[1], nearest_dump_h=float(np.min(np.abs(dd - p[0])) * 24) if len(dd) else None,
                        t_from_start_d=float(p[0] - t.min()), t_to_end_d=float(t.max() - p[0]))
    sd = float(np.std(snr[~nearB]))
    print(f"{name:14s} B boxes: " + ", ".join(f"{b['snr']:.1f}@{b['t_pred']:.2f}" for b in bvals) +
          f" | top non-B: " + ", ".join(f"{p['snr']:.1f}@{p['t']:.2f}" for p in peaks[:5]) + f" | sd(all) {sd:.2f}, ~{n_indep} indep boxes", flush=True)
    return dict(name=name, B_boxes=bvals, top_non_B=peaks, sd_all_positions=sd, n_independent_boxes=n_indep)


def rank_test(r):
    """distribution-free: if the two B boxes are the two strongest independent boxes of the data set, the chance that
    the top two of N exchangeable box positions form a pair consistent with the prediction is K * n_sp / C(N, 2):
    K = positions the first box may take inside the scan window (2 * 0.004 d * |n_cycles| / dur), n_sp = 3 allowed
    spacings (+-1 box) for the second box."""
    b = r["B_boxes"]
    if len(b) != 2 or any(x["rank_among_non_B_peaks"] != 1 for x in b):
        return dict(applicable=False)
    N = r["n_independent_boxes"]
    K = 2 * 0.004 * 227 / DUR
    p = K * 3 / (N * (N - 1) / 2)
    return dict(applicable=True, N_independent_boxes=N, K_window_positions=K, n_spacings=3, p_rank=p)


if __name__ == "__main__":
    out = [census("pipe S6", lcio.load_product("TESSSPOC_S6_PDC"), "TESSSPOC_S6_PDC"),
           census("pipe S32", lcio.load_product("SPOC2min_S32_PDC"), "SPOC2min_S32_PDC"),
           census("own S6", lcio.load_own("tesscut_s6", 2.5)),
           census("own S32", lcio.load_own("tesscut_s32", 2.0)),
           census("pipe S98 frozen", lcio.load_product("QLP_S98_SYS", "frozen")),
           census("pipe S98 default", lcio.load_product("QLP_S98_SYS", "default")),
           census("own S98 nostray", lcio.load_own("tesscut_s98", 2.0, drop_stray=True))]
    for r in out:
        if r["name"].endswith("S6"):
            r["rank_test"] = rank_test(r)
            print(r["name"], "rank test:", r["rank_test"])
    save_json(dict(dur_h=DUR * 24, results=out), "out/s11c_single_box_census.json")
