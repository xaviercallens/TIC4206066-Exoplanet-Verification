"""Flux-scale cross-check between products using signal A (23 transits, same star, same pixels) as a calibrator,
plus direct product-vs-product comparisons in the same sector.

For every product: detrend (biweight 0.75 d, A and B masked), subtract a 1000-ppm B template, fit a common A mid-depth
per sector with the fixed A template (k=0.0212, b=0.35, rho=2.51, LD 0.32/0.31) and per-event linear baselines in
+-0.25 d windows. Also regress each product's detrended flux on the reference product of the same sector
(slope = relative flux scale for short-timescale variations, which should be 1 if both are undiluted).
Output: out/s05_flux_scale.json
"""
from b_run import *
from s_configs import template_for, CONFIGS

tplB = template_for(CONFIGS["BASE"])
PRODUCTS = {6: ["S6_tessspoc_pdc", "S6_tessspoc_sap", "S6_qlp_ksp", "S6_qlp_sap", "S6_tglc_aper", "S6_eleanor_corr", "own6_r1.5", "own6_r2.0", "own6_r2.5", "own6_r3.0"],
            32: ["S32_spoc2m_pdc", "S32_spoc2m_sap", "S32_tessspoc10_pdc", "S32_qlp_ksp", "S32_qlp_sap", "own32_r1.5", "own32_r2.0", "own32_r2.5", "own32_r3.0"],
            98: ["S98_qlp_sysrm", "S98_qlp_det", "S98_qlp_sap", "own98_r1.5", "own98_r2.0", "own98_r2.5", "own98_r3.0"]}


def A_depth(lcd, hw=0.25):
    tA = tplA()
    y_all = lcd["fd"] + 1000e-6 * tplB.unit(phase_dt(lcd["t"], EPH_B), lcd["cad"], True)
    t = lcd["t"]
    n0 = int(np.floor((t.min() - EPH_A["T0"]) / EPH_A["P"])) - 1; n1 = int(np.ceil((t.max() - EPH_A["T0"]) / EPH_A["P"])) + 1
    blocks = []
    for n in range(n0, n1 + 1):
        tc = EPH_A["T0"] + n * EPH_A["P"]; w = np.abs(t - tc) < hw
        inn = w & (np.abs(t - tc) < EPH_A["T14"] / 2)
        if inn.sum() * lcd["cad"] / EPH_A["T14"] >= 0.5 and w.sum() >= 10:
            blocks.append((tc, w))
    N = sum(w.sum() for _, w in blocks); X = np.zeros((N, 2 * len(blocks) + 1)); y = np.zeros(N); i0 = 0
    for j, (tc, w) in enumerate(blocks):
        n = w.sum(); x = (t[w] - tc) / hw
        X[i0:i0 + n, 2 * j] = 1; X[i0:i0 + n, 2 * j + 1] = x
        X[i0:i0 + n, -1] = -tA.unit(t[w] - tc, lcd["cad"], True); y[i0:i0 + n] = y_all[w]; i0 += n
    c, cov, r = lsq(X, y, lcd["sig"])
    return float(c[-1] * 1e6), float(np.sqrt(cov[-1, -1]) * 1e6), len(blocks)


res = dict(A={}, regress={})
for sec, keys in PRODUCTS.items():
    ref = None
    for key in keys:
        lc = get_lc(key, spoc_mask=17087)
        lcd = detrend(lc)
        D, eD, n = A_depth(lcd)
        res["A"][key] = dict(sector=sec, D=D, eD=eD, n=n, sig_ppm=lcd["sig"] * 1e6, cad_min=lcd["cad"] * 1440, meta=lc.get("meta"))
        # regression on the reference product (first in list) after binning both to 30 min on common times
        if ref is None:
            ref = lcd
        else:
            tb = np.arange(max(ref["t"].min(), lcd["t"].min()), min(ref["t"].max(), lcd["t"].max()), 1 / 48)
            ib_r = np.digitize(ref["t"], tb); ib_l = np.digitize(lcd["t"], tb)
            common = np.intersect1d(np.unique(ib_r), np.unique(ib_l))
            xr = np.array([np.mean(ref["fd"][ib_r == i]) for i in common]); yl = np.array([np.mean(lcd["fd"][ib_l == i]) for i in common])
            ok = np.isfinite(xr) & np.isfinite(yl) & (np.abs(xr - 1) < 0.01) & (np.abs(yl - 1) < 0.01)
            xr, yl = xr[ok] - 1, yl[ok] - 1
            # total-least-squares-free estimate: ordinary slope of y on x and of x on y (bracket the true scale)
            s_yx = float(np.sum(xr * yl) / np.sum(xr * xr)); s_xy = float(np.sum(xr * yl) / np.sum(yl * yl))
            res["regress"][key] = dict(ref=ref["name"], n=int(ok.sum()), slope_y_on_x=s_yx, inv_slope_x_on_y=1 / s_xy,
                                       corr=float(np.corrcoef(xr, yl)[0, 1]))
        print(f"S{sec} {key:22s} A depth {D:6.0f} ± {eD:4.0f} (n={n:2d}, sig {lcd['sig']*1e6:5.0f} ppm/{lcd['cad']*1440:.1f} min)"
              + (f"  regress vs {ref['name']}: {res['regress'][key]['slope_y_on_x']:.3f}..{res['regress'][key]['inv_slope_x_on_y']:.3f} r={res['regress'][key]['corr']:.2f}" if key in res["regress"] else ""), flush=True)
jdump(res, "s05_flux_scale.json")
