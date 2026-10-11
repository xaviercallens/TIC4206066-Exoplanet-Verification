"""The single RNAAS figure, built only from this directory's verified outputs (a1-a5) and the raw light curves.
(a) A per-sector folds + template; (b) B per-sector folds; (c) S6 box-S/N census with B's predicted times;
(d) PRF localisation offsets (injection-calibrated) for A and B with the combined S32+S98 result and neighbours;
(e) TRICERATOPS FPP values (this work's fold, the audit's baseline fold, the earlier fold; MOLUSC-limited).
Writes ../figure.pdf, ../figure.png and out/a6_figure_values.json."""
import csv, json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a1_timing_depths as a1
import a3_signalB_tests as a3
H = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(H, "out"); TOP = os.path.dirname(H)
AUD = os.path.abspath(os.path.join(TOP, "..", "tic4206066_PUBLICATION_AUDIT"))
S5 = json.load(open(os.path.join(OUT, "a5_summary.json"))); A1 = json.load(open(os.path.join(OUT, "a1_results.json")))
ephA = A1["signals"]["A"]["A__pipe"]["ephemeris"]["e_dt_inj_min"]; ephB = A1["signals"]["B"]["B_mann__pipe"]["ephemeris"]["e_dt_inj_min"]
EPH = {"A": dict(P=ephA["P"], T0=ephA["T0"], T14=A1["signals"]["A"]["A__pipe"]["template_T14_h"] / 24),
       "B": dict(P=ephB["P"], T0=ephB["T0"], T14=A1["signals"]["B"]["B_mann__pipe"]["template_T14_h"] / 24)}
vals = {"ephemerides_used": EPH}
plt.rcParams.update({"font.size": 7.5, "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True, "font.family": "serif", "mathtext.fontset": "dejavuserif", "pdf.fonttype": 42})
fig = plt.figure(figsize=(7.1, 5.3)); gs = fig.add_gridspec(2, 6, left=0.075, right=0.985, bottom=0.085, top=0.975, wspace=1.35, hspace=0.36)
axa, axb, axc = fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[0, 4:6]); axd, axe = fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6])
COL = {6: "#1f4e9c", 32: "#b2182b", 98: "#1b7837"}; OFF = {6: 0, 32: -1500, 98: -3000}; LAB = {6: "S6, TESS-SPOC 30 min", 32: "S32, SPOC 2 min", 98: "S98, QLP 200 s"}
lcs = {s: a1.detrend(a1.load_lc(k)) for s, k in ((6, "S6"), (32, "S32"), (98, "S98"))}
def binned(x, y, w, sig):
    edges = np.arange(-6, 6 + 1e-9, w); c = 0.5 * (edges[1:] + edges[:-1]); m, e = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (x >= lo) & (x < hi); m.append(np.mean(y[s]) if s.sum() >= 2 else np.nan); e.append(sig / np.sqrt(s.sum()) if s.sum() >= 2 else np.nan)
    return c, np.array(m), np.array(e)
for ax, sig, other in ((axa, "A", "B"), (axb, "B", "A")):
    e = EPH[sig]; tpl = a1.Template(sig, "A" if sig == "A" else "B_mann", e["P"])
    for s in (6, 32, 98):
        L = lcs[s]; t, fd = L["t"], L["fd"]; keep = ~a1.in_tr(t, other, 1.0)
        x = (((t - e["T0"]) / e["P"] + 0.5) % 1 - 0.5) * e["P"] * 24; m = keep & (np.abs(x) < 6); y = (fd[m] - 1) * 1e6; sg = 1.4826 * np.median(np.abs(y - np.median(y)))
        c, mb, eb = binned(x[m], y, 0.5 if s == 6 else 1 / 3, sg); ax.errorbar(c, mb + OFF[s], eb, fmt="o", ms=2.2, lw=0.6, color=COL[s], capsize=0)
        D = S5["B_depths"]["B_mann__pipe"]["per_sector"][str(s)][0] * 1e-6 if sig == "B" else A1["signals"]["A"]["A__pipe"]["per_sector"][str(s)]["e_D_inj"][0]
        tt = np.linspace(-0.25, 0.25, 600); ax.plot(tt * 24, -D * tpl.unit(tt, L["cad"]) * 1e6 + OFF[s], color="k", lw=0.8)
        ax.text(5.5, OFF[s] + 330, LAB[s], ha="right", va="bottom", fontsize=6.3, color=COL[s])
    ax.set_xlim(-6, 6); ax.set_ylim(-4600, 1150); ax.set_xlabel("hours from mid-transit")
    ax.text(0.03, 0.97, f"({'a' if sig == 'A' else 'b'})  P = {e['P']:.2f} d", transform=ax.transAxes, ha="left", va="top", fontweight="bold")
axa.set_ylabel("relative flux (ppm) + offset")
# (c) S6 box census
S6 = a1.load_lc("S6"); t, f, w = a3.prep_A_only(S6); cen = np.arange(t.min() + a3.DUR / 2, t.max() - a3.DUR / 2, a3.DUR / 8); snr = np.full(len(cen), np.nan)
for i, c in enumerate(cen):
    m = np.abs(t - c) < a3.DUR / 2
    if m.sum() >= 2: sw = w[m].sum(); snr[i] = (1 - np.sum(w[m] * f[m]) / sw) * np.sqrt(sw)
axc.plot(cen - 2457000 + 2457000 - 1400, snr, color="#1f4e9c", lw=0.7)
tb = [EPH["B"]["T0"] + n * EPH["B"]["P"] for n in range(-230, -220)]; tb = [x for x in tb if t.min() < x < t.max()]
for x in tb: axc.axvline(x - 1400, color="#b2182b", lw=0.8, ls="--")
for x in S6["t"][[0]]: pass
axc.set_xlabel("BTJD $-$ 1400 (d), Sector 6"); axc.set_ylabel("1.6-h box depth S/N"); axc.set_ylim(-5.5, 8.5)
axc.text(0.03, 0.97, "(c)", transform=axc.transAxes, ha="left", va="top", fontweight="bold")
axc.text(0.31, 0.97, "dashed: 11.13 d transits\npredicted from Sector 98", transform=axc.transAxes, ha="left", va="top", fontsize=6.3, color="#b2182b")
vals["S6_census"] = dict(B_pred_times=[float(x) for x in tb], B_snr=[float(snr[np.argmin(np.abs(cen - x))]) for x in tb], max_other=float(np.nanmax(np.where(np.min(np.abs(cen[:, None] - np.array(tb)[None, :]), axis=1) > a3.DUR, snr, np.nan))))
# (d) localisation
loc = S5["localisation"]; axd.set_aspect("equal"); axd.add_patch(plt.Circle((0, 0), 21, fill=False, ls=":", lw=0.6, color="k")); axd.text(0, 22.5, "21$''$ (1 pixel)", ha="center", va="bottom", fontsize=6)
mk = {"S6_ffi": "o", "S32_tpf": "s", "S98_ffi": "^"}; cc = {"S6_ffi": COL[6], "S32_tpf": COL[32], "S98_ffi": COL[98]}
for sig, fill in (("A", True), ("B", False)):
    for k, v in loc[sig]["per_set"].items():
        if k not in mk: continue
        o, s = v["offset_arcsec"], v["sigma_arcsec"]
        axd.errorbar(o[0], o[1], xerr=s[0], yerr=s[1], fmt=mk[k], ms=4, mfc=cc[k] if fill else "white", mec=cc[k], color=cc[k], lw=0.7, capsize=0)
    c = loc[sig]["combined_S32tpf_S98clean"]; o, s = c["offset_arcsec"], c["sigma_arcsec"]
    axd.add_patch(Ellipse((o[0], o[1]), 2 * s[0], 2 * s[1], fill=False, lw=1.0, ls="-" if sig == "A" else "--", color="k"))
BL = json.load(open(os.path.join(TOP, "..", "deep_followup_tic4206066", "out", "blend_check.json")))
nb = [n for n in BL["able_to_mimic"] if n["sep"] < 62]; axd.plot([n["dRA"] for n in nb], [n["dDec"] for n in nb], "x", color="k", ms=4.5, mew=0.9); axd.plot(0, 0, "+", color="k", ms=8, mew=1.0)
axd.set_xlim(62, -76); axd.set_ylim(-66, 52); axd.set_xlabel(r"$\Delta\alpha$ (arcsec, east left)"); axd.set_ylabel(r"$\Delta\delta$ (arcsec)")
hand = [Line2D([], [], marker="o", ls="none", mfc="k", mec="k", ms=4, label="3.18 d (filled)"), Line2D([], [], marker="o", ls="none", mfc="white", mec="k", ms=4, label="11.13 d (open)"),
        Line2D([], [], color="k", lw=1.0, label="S32+S98 combined, 1$\\sigma$"), Line2D([], [], marker="x", ls="none", color="k", ms=4.5, label="Gaia DR3 stars (see caption)")]
axd.legend(handles=hand, loc="lower right", frameon=False, handletextpad=0.3, borderaxespad=0.2); axd.text(0.03, 0.97, "(d)", transform=axd.transAxes, ha="left", va="top", fontweight="bold")
# (e) FPP
fp = {}
for tag in ("A_base_S32+S98", "A_molusc5_S32+S98", "B_base_S32+S98", "B_molusc5_S32+S98", "B_prevfold_repro_S32+S98", "A_repro_S32+S98", "B_tic_S32+S98", "A_tic_S32+S98"):
    j = json.load(open(os.path.join(AUD, "out", "fpp", "runs", tag, "summary.json")))["log"]; fp[tag] = (j["FPP_mean"], j["FPP_std"], j["FPP_min"], j["FPP_max"])
own = {}
for sig in ("A", "B"):
    f = os.path.join(OUT, "fpp", f"fpp_ownfold_{sig}.json")
    if os.path.exists(f): j = json.load(open(f)); own[sig] = (j["FPP_mean"], j["FPP_sd"] or 0, j["FPP_min"], j["FPP_max"])
xs = {"own": 0, "base": 1, "prev": 2, "tic": 3, "molusc": 4}; names = ["this fold", "audit fold", "earlier fold", "TIC params", "+ MOLUSC"]
for sig, col, dx, m in (("A", "k", -0.08, "o"), ("B", "#b2182b", 0.08, "D")):
    pts = {"base": fp[f"{sig}_base_S32+S98"], "prev": fp["B_prevfold_repro_S32+S98"] if sig == "B" else fp["A_repro_S32+S98"], "tic": fp[f"{sig}_tic_S32+S98"], "molusc": fp[f"{sig}_molusc5_S32+S98"]}
    if sig in own: pts["own"] = own[sig]
    for k, (mu, sd, lo, hi) in pts.items():
        axe.vlines(xs[k] + dx, lo, hi, color=col, lw=0.6); axe.errorbar(xs[k] + dx, mu, sd, fmt=m, ms=3.5, color=col, mfc=col, lw=1.2, capsize=0, label=f"{'3.18' if sig == 'A' else '11.13'} d" if k == "base" else None)
    vals[f"fpp_{sig}"] = {k: dict(mean=v[0], sd=v[1], min=v[2], max=v[3]) for k, v in pts.items()}
axe.set_yscale("log"); axe.set_ylim(1.5e-4, 2); axe.axhline(0.015, color="k", lw=0.6, ls="--"); axe.text(4.45, 0.017, "0.015", ha="right", va="bottom", fontsize=6.3)
axe.axhline(0.5, color="k", lw=0.6, ls=":"); axe.text(4.45, 0.56, "0.5", ha="right", va="bottom", fontsize=6.3)
axe.set_xticks(range(5)); axe.set_xticklabels(names, fontsize=6.5); axe.set_xlim(-0.5, 4.5); axe.set_ylabel("TRICERATOPS FPP"); axe.legend(loc="lower left", frameon=False)
axe.text(0.03, 0.97, "(e)", transform=axe.transAxes, ha="left", va="top", fontweight="bold")
fig.savefig(os.path.join(TOP, "figure.pdf")); fig.savefig(os.path.join(TOP, "figure.png"), dpi=200)
json.dump(vals, open(os.path.join(OUT, "a6_figure_values.json"), "w"), indent=1, default=float); print(json.dumps(vals["S6_census"]), "figure written")
