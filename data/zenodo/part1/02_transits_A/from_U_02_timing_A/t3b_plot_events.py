"""Diagnostic plot: every A window (detrended flux, binned), best-fit template, and the normalised chi2(dt) curve."""
import sys
import numpy as np, pandas as pd, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import t2_common as C

variant = sys.argv[1] if len(sys.argv) > 1 else "base"
df = pd.read_csv(f"{C.OUT}/t3_events_{variant}.csv")
meta = json.load(open(f"{C.OUT}/t3_meta_{variant}.json"))
P, T0 = meta["ephem_window"]["P"], meta["ephem_window"]["T0"]
tpl = C.Template(); T14 = tpl.t14
lcs = {s: C.load_lc(s) for s in (6, 32, 98)}
n = len(df); ncol = 4; nrow = int(np.ceil(n / ncol))
fig, axs = plt.subplots(nrow * 2, ncol, figsize=(16, 3.2 * nrow * 2))
for i, r in df.iterrows():
    lc = lcs[r.sector]
    _, tcs_all = C.a_epochs(lc["t"], P, T0)
    tt, yy, cov = C.prepare_window(lc, r.tc_pred, tcs_all, T14)
    x = tt - r.tc_pred
    sc = C.Scanner(x, lc["cad"], tpl)
    chi, Cc = sc.chi2(yy)
    j = int(np.argmin(chi)); s2 = chi[j] / (len(x) - 4)
    ax = axs[(i // ncol) * 2, i % ncol]
    ax.plot(x * 24, (yy - 1) * 1e6, ".", ms=2, alpha=0.4)
    nb = max(1, int(round(20 / (lc["cad"] * 1440))))
    if nb > 1:
        k = len(x) // nb
        ax.plot(x[:k * nb].reshape(k, nb).mean(1) * 24, (yy[:k * nb].reshape(k, nb).mean(1) - 1) * 1e6, "k.", ms=5)
    c = Cc[j]
    ax.plot(x * 24, (c[0] + c[1] * x - c[2] * sc.S[j] - 1) * 1e6, "r-")
    ax.axvline(0, color="g", lw=0.5); ax.set_ylim(-2500, 1500) if r.sector != 6 else ax.set_ylim(-1200, 700)
    ax.set_title(f"S{r.sector} E{r.epoch} cov{r.coverage:.2f} dt{r["dt"]:.1f} snr{r.snr:.1f} sel{int(r.sel)}", fontsize=8)
    ax2 = axs[(i // ncol) * 2 + 1, i % ncol]
    ax2.plot(sc.dts * 1440, (chi - chi[j]) / s2, "k-", lw=0.8)
    ax2.set_ylim(-0.5, 30); ax2.axvline(r["dt"], color="r", lw=0.5); ax2.axvspan(-60, 60, color="y", alpha=0.1)
    ax2.set_xlabel("dt (min)", fontsize=7)
for k in range(n, nrow * ncol):
    axs[(k // ncol) * 2, k % ncol].axis("off"); axs[(k // ncol) * 2 + 1, k % ncol].axis("off")
fig.tight_layout(); fig.savefig(f"{C.OUT}/t3_events_{variant}.png", dpi=80)
print("saved")
