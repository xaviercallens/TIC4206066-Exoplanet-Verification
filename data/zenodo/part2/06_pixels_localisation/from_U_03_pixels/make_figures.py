"""Audit figures: (1) localisation offsets (panel-d analogue) and (2) difference images."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import pandas as pd
from run_real import *
from poscorr_decorr import decorrelate

S = json.load(open(os.path.join(OUT, "summary.json")))
sets, comb = S["sets"], S["combined"]


def ell(ax, mu, C, **kw):
    C = np.array(C); w, v = np.linalg.eigh(C)
    ang = np.degrees(np.arctan2(v[1, 1], v[0, 1]))
    ax.add_patch(Ellipse(mu, 2 * np.sqrt(w[1]), 2 * np.sqrt(w[0]), angle=ang, fill=False, **kw))


fig, axs = plt.subplots(1, 2, figsize=(11, 5.4))
for ax, sig in zip(axs, "AB"):
    ax.set_aspect("equal")
    ax.add_patch(plt.Circle((0, 0), 21, fill=False, ls=":", lw=0.6, color="k"))
    items = [(f"S6_tspoc_poscorr_{sig}", "o", "tab:blue", "S6 (POS_CORR-decorr.)"),
             (f"S32_tpf_{sig}", "s", "0.6", "S32 2-min, no decorr."),
             (f"S32_tpf_poscorr_{sig}", "s", "tab:red", "S32 2-min, POS_CORR-decorr."),
             (f"S98_clean29_{sig}", "^", "tab:green", "S98 clean (2048+2^29)"),
             (f"S98_clean_{sig}", "v", "yellowgreen", "S98 clean (2048 only)")]
    for k, m, c, lab in items:
        if k not in sets:
            continue
        s = sets[k]
        ax.errorbar(s["real_off"][0], s["real_off"][1], xerr=s["inj_std"][0], yerr=s["inj_std"][1], fmt=m, color=c,
                    ms=5, lw=0.8, label=lab)
    c = comb[f"{sig}:S32tpfPC+S98clean29:bias0:rob0"]
    ell(ax, c["mean"], c["cov"], color="k", lw=1.4)
    ax.plot([], [], "k-", label="this audit S32(decorr)+S98, 1$\\sigma$")
    man = {"A": ((8, -8), (5, 7)), "B": ((1, 1), (3, 4))}[sig]
    ax.add_patch(Ellipse(man[0], 2 * man[1][0], 2 * man[1][1], fill=False, ls="--", color="m", lw=1.2))
    ax.plot([], [], "m--", label="manuscript, 1$\\sigma$")
    nb = pd.read_csv(os.path.join(OUT, "neighbours_vs_fitted_position.csv"))
    cap = nb[nb["dreqA_min"] <= 1] if sig == "A" else nb[nb["dreqB_min"] <= 1]
    ax.plot(cap["dE"], cap["dN"], "kx", ms=7, mew=1.2, label="Gaia DR3 able to give depth")
    n11 = nb[nb["sep"] < 15]
    ax.plot(n11["dE"], n11["dN"], "+", color="0.4", ms=7, label="other Gaia <15\"")
    ax.plot(0, 0, "k+", ms=10)
    ax.set_xlim(95, -95); ax.set_ylim(-70, 60)
    ax.set_xlabel("offset East (arcsec)"); ax.set_ylabel("offset North (arcsec)")
    ax.set_title({"A": "3.18 d signal", "B": "11.13 d signal"}[sig])
    ax.legend(fontsize=6.5, loc="upper right", framealpha=0.6)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_localisation.png"), dpi=140)

# difference images
fig, axs = plt.subplots(2, 4, figsize=(13, 6.4))
for row, sig in enumerate("AB"):
    for col, (cfg, dec) in enumerate([("S32_tpf", False), ("S32_tpf", True), ("S98_clean29", False), ("S6_tspoc", True)]):
        name, key, var, kw = [c for c in CONFIGS if c[0] == cfg][0]
        ds = Dataset(key, var, **kw)
        if dec:
            decorrelate(ds)
        sf = ds.scene_fit()
        out, Dm, Ve, _ = analyse(ds, sig, sf, want_nb=False)
        ax = axs[row, col]
        snr = Dm / np.sqrt(Ve)
        im = ax.imshow(snr, origin="lower", cmap="RdBu_r", vmin=-5, vmax=5)
        ax.plot(sf["x_t"], sf["y_t"], "k+", ms=10)
        ax.plot(out["fit"]["x"], out["fit"]["y"], "mx", ms=8)
        cap = ds.gaia["Tmag_used"].values < 17.5
        ax.plot(ds.g_x[cap] + sf["dx"], ds.g_y[cap] + sf["dy"], "o", mfc="none", mec="0.3", ms=4)
        ax.set_xlim(-0.5, ds.cube.nx - 0.5); ax.set_ylim(-0.5, ds.cube.ny - 0.5)
        ax.set_title(f"{sig}: {cfg}{' POS_CORR-decorr.' if dec else ''}\noff=({out['off_E']:+.1f},{out['off_N']:+.1f})\" S/N={out['snr_target']:.1f}",
                     fontsize=8)
        plt.colorbar(im, ax=ax, fraction=0.046)
fig.suptitle("Difference-image S/N (in minus out; + = flux loss). +: target (OOT scene fit); x: PRF fit; o: Gaia T<17.5", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_diffimages.png"), dpi=130)
print("figures written")
