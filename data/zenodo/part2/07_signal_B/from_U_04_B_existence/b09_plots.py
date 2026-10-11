"""Diagnostic figure for Audit 4."""
import numpy as np, json, glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import common_b as c
import archdata as ad
import s98data as sd

fig, ax = plt.subplots(3, 3, figsize=(15, 11))
# (a) S98 six predicted epochs (QLP DET, all cadences except 2048/dump; colour = flag class)
q = c.load_qlp(98)
t, f, qq = q["t"], q["f"], q["q"]
ns, tcs = c.pred_times(c.P_B98, c.T0_B98, t.min(), t.max())
for k, (n, tc) in enumerate(zip(ns, tcs)):
    m = np.abs(t - tc) < 0.3
    base = np.nanmedian(f[m & (np.abs(t - tc) > 0.07)])
    y = f[m] / base - 1
    cl = qq[m] == 0
    off = -k * 4e-3
    a = ax[0, 0]
    a.plot((t[m][~cl] - tc) * 24, y[~cl] + off, ".", ms=2, color="tab:orange", alpha=0.5)
    a.plot((t[m][cl] - tc) * 24, y[cl] + off, ".", ms=2, color="k", alpha=0.5)
    # 30-min bins of q==0
    tb = np.arange(-7, 7.01, 0.5)
    yy = [np.nanmean(y[cl][(np.abs((t[m][cl] - tc) * 24 - x) < 0.25)]) if np.sum(np.abs((t[m][cl] - tc) * 24 - x) < 0.25) > 3 else np.nan for x in tb]
    a.plot(tb, np.array(yy) + off, "-", color="tab:blue")
    a.text(-7, off + 1e-3, "E%d %.2f" % (n, tc), fontsize=7)
ax[0, 0].axvspan(-0.8, 0.8, color="tab:green", alpha=0.1)
ax[0, 0].set_ylim(-0.024, 0.004); ax[0, 0].set_xlabel("h from predicted B (S98 search ephemeris)")
ax[0, 0].set_title("(a) S98 QLP DET: black q=0, orange flagged", fontsize=9)
# (b) S98 BLS spectra
for var, col in [("QLP_DET_bw075", "k"), ("OWN_bhalves_bw075", "tab:red")]:
    tt, ff = sd.get(var)
    P = sd.period_grid(tt.max() - tt.min())
    r, pw, snr = sd.bls_spectrum(tt, ff, P)
    ax[0, 1].plot(P, sd.sde(pw), color=col, lw=0.6, label=var)
ax[0, 1].axvline(11.133, color="tab:green", ls=":"); ax[0, 1].legend(fontsize=7); ax[0, 1].set_xlabel("P (d)"); ax[0, 1].set_ylabel("SDE")
ax[0, 1].set_title("(b) S98-only BLS (A removed)", fontsize=9)
# (c) S98 null distributions
for fn in sorted(glob.glob("out/b03_null_*.json")):
    d = json.load(open(fn))
    ax[0, 2].hist(d["sde_max"], bins=50, histtype="step", density=True, label="%s %s (n=%d)" % (d["var"][:8], d["family"], len(d["sde_max"])))
    ax[0, 2].axvline(d["obs"]["SDE_B"], ls="--", lw=0.8)
ax[0, 2].legend(fontsize=6); ax[0, 2].set_xlabel("max SDE (null) / observed B SDE (dashed)"); ax[0, 2].set_title("(c) S98-only null", fontsize=9)
# (d-f) S6 box S/N series: pipe, own, tglc_psf
TB = [1473.7330, 1484.8657]
for a, ph in zip([ax[1, 0], ax[1, 1], ax[1, 2]], ["pipe", "own", "tglc_psf"]):
    tt, ff = ad.get(6, ph)
    g = ad.BoxGrid(tt, ff)
    x, s = g.snr_series()
    a.plot(x, s, lw=0.5, color="k")
    for tb in TB: a.axvline(tb, color="tab:green", ls="--", lw=0.8)
    for dmp in c.DUMPS[6]: a.axvline(dmp, color="tab:red", ls=":", lw=0.6)
    a.set_title("S6 1.6-h box S/N, %s (green: predicted B, red: dumps)" % ph, fontsize=9); a.set_xlabel("BTJD")
# (g) link-test S/N vs P
d = json.load(open("out/b04_link_main.json"))
for r in d:
    if r["null"] == "kept" and r["phot"] in ("pipe", "own"):
        ax[2, 0].plot(r["curve"]["P"], r["curve"]["snr"], lw=0.7, label="%s S%s" % (r["phot"], "+".join(map(str, r["sectors"]))))
ax[2, 0].axvline(c.P_B_MS, color="tab:green", ls=":"); ax[2, 0].legend(fontsize=7); ax[2, 0].set_xlabel("P (d)"); ax[2, 0].set_ylabel("combined box S/N")
ax[2, 0].set_title("(g) archival link scan (W1: +-0.004 d)", fontsize=9)
# (h) S32 box S/N series
tt, ff = ad.get(32, "pipe"); g = ad.BoxGrid(tt, ff); x, s = g.snr_series()
ax[2, 1].plot(x, s, lw=0.4, color="k")
for tb in [2175.0963, 2186.2290, 2197.3617]: ax[2, 1].axvline(tb, color="tab:green", ls="--", lw=0.8)
for dmp in c.DUMPS[32]: ax[2, 1].axvline(dmp, color="tab:red", ls=":", lw=0.6)
ax[2, 1].set_title("(h) S32 2-min PDCSAP box S/N", fontsize=9)
# (i) O-C
e = json.load(open("out/b08_ephem.json"))
E = e["ephemerides"]["S6x2 + S32x2 + S98(0,1,4)"]
ev = [x for x in e["events"] if x["used"] and not x.get("at_edge")]
for x in ev:
    oc = (x["tc"] - (E["T0"] + x["epoch"] * E["P"])) * 1440
    ax[2, 2].errorbar(x["tc"], oc, yerr=x["err_min"], fmt="o", label="S%s E%d" % (x["sector"], x["epoch"]))
ax[2, 2].axhline(0, color="k", lw=0.5); ax[2, 2].set_ylabel("O-C (min)"); ax[2, 2].legend(fontsize=6)
ax[2, 2].set_title("(i) B O-C vs S6+S32+S98 linear ephemeris", fontsize=9)
plt.tight_layout()
plt.savefig("b09_diagnostics.png", dpi=110)
print("saved")
