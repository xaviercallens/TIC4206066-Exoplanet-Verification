"""BEFORE: first-pass method (QUALITY == 0 only, NO background subtraction of the TESScut cubes, no scattered-light cut) (global running-median detrend, no pointing correction).
Per-pixel and per-aperture fractional depth for A and B; errors from a random-phase null (fixed seed).
This is the method that produced the apparent B/S32 'off-target' effect; kept for provenance."""
import json, os
import numpy as np
from common import OUT, EPH, RADII, load, phase_dist, radii_grid

rng = np.random.default_rng(1)
WIN = 0.75  # d, running-median window

def detrend(f, t, mask):
    out = np.full_like(f, np.nan)
    seg = np.concatenate([[0], np.where(np.diff(t) > 0.5)[0] + 1, [len(t)]])
    for a, b in zip(seg[:-1], seg[1:]):
        tt, ff = t[a:b], f[a:b].copy(); ff[mask[a:b]] = np.nan
        for i in range(b - a):
            w = np.abs(tt - tt[i]) < WIN / 2
            out[a + i] = np.nanmedian(ff[w]) if np.sum(np.isfinite(ff[w])) > 5 else np.nan
    return f / out

res = {}
for name in ("S6_ffi", "S32_ffi", "S32_tpf2min", "S98_ffi"):
    t, F, bk, c, med = load(name, raw=True)   # first-pass selection, kept for provenance
    dt = np.nanmedian(np.diff(t)); nb = max(1, int(round((10 / 1440) / dt)))   # bin to ~10 min
    if nb > 1:
        n = (len(t) // nb) * nb; t = t[:n].reshape(-1, nb).mean(1); F = F[:n].reshape(-1, nb, *F.shape[1:]).mean(1)
    mask = np.zeros(len(t), bool)
    for s in "AB":
        P, T0, T14 = EPH[s]; mask |= np.abs(phase_dist(t, P, T0)) < 0.75 * T14
    rr = radii_grid(F.shape[1:], c)
    lcs = {r: detrend(F[:, rr <= r].sum(1), t, mask) for r in RADII}
    def stat(lc1, lc2, s):
        P, T0, T14 = EPH[s]; core = np.abs(phase_dist(t, P, T0)) < 0.25 * T14
        f = lambda lc, w: 1 - np.nanmean(lc[w])
        d0 = f(lc1, core) - (f(lc2, core) if lc2 is not None else 0)
        nl = []
        for _ in range(400):
            w = (np.abs(phase_dist(t, P, T0 + rng.uniform(0.15, 0.85) * P)) < 0.25 * T14) & ~mask
            if w.sum() >= 0.5 * core.sum(): nl.append(f(lc1, w) - (f(lc2, w) if lc2 is not None else 0))
        return [float(d0 * 1e6), float(np.nanstd(nl) * 1e6)]
    out = {"centre_pix": c, "n_cad_10min": int(len(t)),
           "apertures": {str(r): {s: stat(lcs[r], None, s) for s in "AB"} for r in RADII},
           "diffs": {f"{a}-{b}": {s: stat(lcs[a], lcs[b], s) for s in "AB"} for a, b in ((0.5, 2.0), (1.0, 3.0), (2.0, 3.0))},
           "pixels": {}}
    for j in range(c[0] - 2, c[0] + 3):
        for i in range(c[1] - 2, c[1] + 3):
            lc = detrend(F[:, j, i], t, mask)
            out["pixels"][f"{j-c[0]:+d},{i-c[1]:+d}"] = {"flux_frac": float(med[j, i] / med[c]), **{s: stat(lc, None, s) for s in "AB"}}
    res[name] = out
    print(name, {r: (round(v["A"][0]), round(v["B"][0])) for r, v in out["apertures"].items()})
res["_note"] = ("S98_ffi result of this method is unphysical (all pixels 'brighten' by 3000-40000 ppm in A windows, growing with "
                "aperture area); cause not identified; not used for any conclusion.")
json.dump(res, open(os.path.join(OUT, "before_global_detrend.json"), "w"), indent=1)
