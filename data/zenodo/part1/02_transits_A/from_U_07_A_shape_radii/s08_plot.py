"""Diagnostic plot: folded detrended (wotan) and local-baseline-corrected light curves with the two posterior-median
models, and the secondary-eclipse phase scan."""
import json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import *
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for tag, dv, col in (("mann_claret_jit_base", "base", "C0"), ("mann_claret_jit_raw", "raw", "C3")):
    ch = np.load(os.path.join(OUT, f"chain_{tag}.npz"))["chain"]
    T0, P, k, b, lnrho, u1, u2 = np.median(ch, axis=0)[:7]
    ars = ars_from_rho(np.exp(lnrho), P)
    allx, ally = [], []
    for s in (6, 32, 98):
        z = np.load(os.path.join(OUT, f"lc_{dv}_S{s}.npz")); t, f = z["t"], z["f"]
        dt, ep = phase_dt(t, P, T0); dB, _ = phase_dt(t, P_B, T0_B)
        m = (np.abs(dt) < 0.3) & (np.abs(dB) > T14_B_H / 24)
        t, f, dt, ep = t[m], f[m], dt[m], ep[m]
        if dv == "raw":
            mod = BatmanSet(t, float(z["cad"]), SUPERSAMPLE[s]).flux(T0, P, k, b, ars, u1, u2)
            for e in np.unique(ep):
                mm = ep == e
                c = np.polyfit(dt[mm], f[mm] - mod[mm], 1)
                f[mm] = f[mm] - np.polyval(c, dt[mm])
        allx.append(dt * 24); ally.append(f)
    x = np.concatenate(allx); y = np.concatenate(ally)
    bins = np.arange(-7.2, 7.21, 0.3)
    idx = np.digitize(x, bins)
    bx = [x[idx == i].mean() for i in range(1, len(bins)) if (idx == i).sum() > 5]
    by = [(y[idx == i].mean() - 1) * 1e6 for i in range(1, len(bins)) if (idx == i).sum() > 5]
    ax[0].plot(bx, by, "o", color=col, ms=4, label=f"{dv}: k={k:.4f}")
    tt = np.linspace(-0.3, 0.3, 2000)
    import batman
    pm = batman.TransitParams(); pm.t0 = 0; pm.per = P; pm.rp = k; pm.a = ars; pm.inc = np.degrees(np.arccos(b / ars))
    pm.ecc = 0; pm.w = 90; pm.u = [u1, u2]; pm.limb_dark = "quadratic"
    ax[0].plot(tt * 24, (batman.TransitModel(pm, tt).light_curve(pm) - 1) * 1e6, "-", color=col)
ax[0].set_xlabel("hours from mid-transit"); ax[0].set_ylabel("relative flux [ppm]"); ax[0].legend(fontsize=8)
ax[0].set_title("Signal A, all sectors (18-min bins)")
for wf in ("0.5", "1.0", "2.0"):
    d = np.loadtxt(os.path.join(OUT, f"secondary_scan_w{wf}.txt"))
    ax[1].plot(d[:, 0], d[:, 1], lw=0.8, label=f"box {wf} x T14")
ax[1].axvline(0.5, color="k", ls=":"); ax[1].set_xlabel("phase"); ax[1].set_ylabel("box flux [ppm]"); ax[1].legend(fontsize=8)
ax[1].set_title("secondary-eclipse scan (base detrending)")
plt.tight_layout(); plt.savefig(os.path.join(OUT, "fold_and_secondary_scan.png"), dpi=120)
