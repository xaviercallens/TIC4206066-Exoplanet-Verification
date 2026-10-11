"""Figures: A depth vs aperture before/after pointing correction; B/S32 before vs per-event; B/S32 event light curves."""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import OUT, PLOTS, EPH, load, radii_grid, centroid

bef = json.load(open(os.path.join(OUT, "before_global_detrend.json")))
aft = json.load(open(os.path.join(OUT, "after_decorrelated.json")))
bev = json.load(open(os.path.join(OUT, "b_s32_events.json")))

# 1. Signal A: depth vs aperture radius, before vs after
sets = ["S6_ffi", "S32_tpf2min", "S32_ffi", "S98_ffi"]
fig, axs = plt.subplots(1, 4, figsize=(17, 4.2), sharey=False)
for ax, s in zip(axs, sets):
    r_b = [float(r) for r in bef[s]["apertures"]]; d_b = [bef[s]["apertures"][r]["A"][0] for r in bef[s]["apertures"]]
    e_b = [bef[s]["apertures"][r]["A"][1] for r in bef[s]["apertures"]]
    key = "A" if "corrected" in aft["A"][s]["summary"] else "A_sensitivity_no_dump_rule"
    c = aft[key][s]["summary"]["corrected"]; r_a = [0.5, 1.0, 2.0, 3.0]
    d_a = [c[str(r)]["depth_ppm"] for r in r_a]; e_a = [c[str(r)]["err_event_scatter"] for r in r_a]
    if s != "S98_ffi":
        ax.errorbar(np.array(r_b) - 0.04, d_b, e_b, fmt="o-", label="before: global detrend")
    else:
        ax.text(0.5, 0.92, "before: unphysical (not shown)", transform=ax.transAxes, ha="center", fontsize=8)
    ax.errorbar(np.array(r_a) + 0.04, d_a, e_a, fmt="s-", label="after: per-event pointing correction")
    note = "" if key == "A" else " [no dump rule: all S6 events near dumps]"
    ax.set_title(f"Signal A, {s} (n={aft[key][s]['summary']['n_events']} after){note}", fontsize=9)
    ax.set_xlabel("aperture radius, px"); ax.axhline(0, lw=0.5)
axs[0].set_ylabel("transit depth, ppm"); axs[0].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "A_depth_vs_aperture_before_after.png"), dpi=100)

# 2. Signal B, S32 2-min: before (global) vs per-event local
fig, ax = plt.subplots(figsize=(6, 4.2))
r_b = [float(r) for r in bef["S32_tpf2min"]["apertures"]]
ax.errorbar(np.array(r_b) - 0.05, [bef["S32_tpf2min"]["apertures"][r]["B"][0] for r in bef["S32_tpf2min"]["apertures"]],
            [bef["S32_tpf2min"]["apertures"][r]["B"][1] for r in bef["S32_tpf2min"]["apertures"]], fmt="o-", label="before: global detrend (both events)")
for k, (ep, e) in enumerate(bev["events"].items()):
    rs = [float(r) for r in e["apertures"]]
    ax.errorbar(np.array(rs) + 0.05*(k+1), [v[0] for v in e["apertures"].values()], [v[1] for v in e["apertures"].values()], fmt="s-", label=f"per-event local, epoch {ep}")
ax.axhline(0, lw=0.5); ax.set_xlabel("aperture radius, px"); ax.set_ylabel("depth, ppm"); ax.set_title("Signal B, S32 2-min TPF")
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "B_S32_before_vs_per_event.png"), dpi=100)

# 3. B/S32 event light curves: core pixel vs r=3 aperture, and centroid drift
t, F, bk, c, med = load("S32_tpf2min", raw=True); rr = radii_grid(F.shape[1:], c); cx, cy = centroid(F, c)
P, T0, T14 = EPH["B"]
evs = list(bev["events"].items())
fig, axs = plt.subplots(len(evs), 3, figsize=(15, 3.8*len(evs)))
for k, (ep, e) in enumerate(evs):
    tc = e["t_mid_BTJD"]; w = np.abs(t - tc) < 0.35; x = (t[w] - tc) * 24
    for ax, (lab, lc) in zip(axs[k, :2], [("core pixel", F[w][:, c[0], c[1]]), ("r = 3 px aperture", F[w][:, rr <= 3].sum(1))]):
        oot = np.abs(x) > 0.75 * T14 * 24; nrm = lc / np.median(lc[oot])
        nb = 15; n = len(x) // nb * nb
        ax.plot(x, nrm, ",", color="0.6"); ax.plot(x[:n].reshape(-1, nb).mean(1), nrm[:n].reshape(-1, nb).mean(1), "o-", ms=3)
        ax.axvspan(-T14*12, T14*12, alpha=0.15); ax.axvline(-e["hours_after_nearest_A_mid"], color="r", ls=":", label="A mid-transit")
        ax.set_xlim(x.min(), x.max()); ax.set_title(f"B epoch {ep}: {lab}"); ax.set_xlabel("hours from B mid-transit"); ax.legend(fontsize=7)
    n = len(x) // 15 * 15
    axs[k, 2].plot(x[:n].reshape(-1, 15).mean(1), 1e3*(cx[w][:n].reshape(-1, 15).mean(1) - cx[w].mean()), "o-", ms=3, label="x")
    axs[k, 2].plot(x[:n].reshape(-1, 15).mean(1), 1e3*(cy[w][:n].reshape(-1, 15).mean(1) - cy[w].mean()), "s-", ms=3, label="y")
    axs[k, 2].axvspan(-T14*12, T14*12, alpha=0.15); axs[k, 2].set_title("flux-weighted centroid drift, mpix"); axs[k, 2].legend(fontsize=7)
    axs[k, 2].set_xlabel("hours from B mid-transit")
fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "B_S32_events_core_vs_wide.png"), dpi=90)
print("plots written")
