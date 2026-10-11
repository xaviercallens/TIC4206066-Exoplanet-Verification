"""Audit 2: figures. (1) injection offset distributions per sector (R144 selected, R60), (2) z distributions,
(3) O-C diagram for R144 and R60 configurations, (4) TTV periodogram and upper limits (if t6 finished)."""
import json, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import t2_common as C

ev = pd.read_csv(f"{C.OUT}/t3_events_base.csv")
t5 = json.load(open(f"{C.OUT}/t5_results.json"))
inj = {s: pd.read_csv(f"{C.OUT}/t4_inj_S{s}.csv") for s in (6, 32, 98)}
fig, axs = plt.subplots(2, 3, figsize=(15, 8))
for i, (s, d) in enumerate(inj.items()):
    sel = (d.coverage >= 0.7) & (d.snr >= 2) & (~d.edge)
    sel60 = (d.coverage >= 0.7) & (d.snr_alt >= 2) & (~d.edge_alt)
    ax = axs[0, i]
    b = np.arange(-145, 146, 5)
    ax.hist(d.dt[sel], b, histtype="step", label=f"R144 selected (SD {d.dt[sel].std():.1f}, rob {1.4826*np.median(np.abs(d.dt[sel]-d.dt[sel].median())):.1f})")
    ax.hist(d.dt_alt[sel60], b, histtype="step", label=f"R60 selected (SD {d.dt_alt[sel60].std():.1f})")
    ax.hist(d.dt[~sel & (d.coverage >= 0.7)], b, histtype="step", label="R144 rejected (S/N<2 or edge)")
    ax.set_yscale("log"); ax.set_title(f"S{s}: t_rec - t_inj (min), D_inj {d.D_inj.iloc[0]*1e6:.0f} ppm"); ax.legend(fontsize=7)
    ax = axs[1, i]
    bz = np.linspace(-8, 8, 81)
    for k, lab in (("e_white", "white"), ("e_pb", "prayer bead"), ("e_boot", "bootstrap")):
        ax.hist((d.dt / d[k])[sel], bz, histtype="step", density=True, label=lab)
    ax.hist((d.dt / d.dt[sel].std())[sel], bz, histtype="step", density=True, label="const sector SD")
    xx = np.linspace(-8, 8, 400); ax.plot(xx, np.exp(-xx ** 2 / 2) / np.sqrt(2 * np.pi), "k--", label="N(0,1)")
    ax.set_yscale("log"); ax.set_ylim(1e-4, 1); ax.set_title(f"S{s}: z = dt / sigma (selected)"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f"{C.OUT}/t7_injections.png", dpi=90)

# O-C
fig, axs = plt.subplots(1, 2, figsize=(14, 4.5))
for ax, (lab, key, tcol, errs) in zip(axs, (("R144 (untruncated search), inj SD errors", "inj_SD_selected_R144 (per-sector const)", "tc_obs", None),
                                            ("R60 (manuscript-like range)", None, "tc_obs_alt", None))):
    if key:
        e = t5["ephemeris_20"][key]; rows = ev[ev.sel]; oc = np.array(e["oc_min"]); err = np.array(e["errors_min_by_event"])
    else:
        rows = ev[(ev.coverage >= 0.7) & (ev.snr_alt >= 2) & (~ev.edge_alt)]
        e = t5["added_events_and_variants"]["search_R60_sel60"]
        sd60 = {int(s): t5["injections"][s]["selected_R60"]["sd"] for s in t5["injections"]}
        err = rows.sector.map(sd60).values
        oc = (rows[tcol].values - (e["T0"] + e["P"] * rows.epoch.values)) * 1440
    for s, m in (("6", "o"), ("32", "s"), ("98", "^")):
        k = rows.sector.values == int(s)
        ax.errorbar(rows.epoch.values[k], oc[k], err[k], fmt=m, label=f"S{s}")
    ax.axhline(0, color="k", lw=0.5); ax.set_xlabel("epoch"); ax.set_ylabel("O-C (min)")
    ax.set_title(f"{lab}: P={e['P']:.7f}+-{e['eP']:.7f}, chi2={e['chi2']:.1f}/{e['dof']}", fontsize=9); ax.legend()
fig.tight_layout(); fig.savefig(f"{C.OUT}/t7_oc.png", dpi=90)

if os.path.exists(f"{C.OUT}/t6_results.json"):
    fig, axs = plt.subplots(2, 2, figsize=(14, 8))
    for j, lab in enumerate(("R144", "R60")):
        z = np.load(f"{C.OUT}/t6_periodogram_{lab}.npz")
        ax = axs[0, j]; ax.plot(1 / z["freqs"], z["dchi"], "k-", lw=0.6); ax.set_xscale("log")
        ax.axhline(np.percentile(z["mx_perm"], 95), color="r", ls="--", label="5% FAP (permutation)")
        ax.axhline(np.percentile(z["mx_emp"], 95), color="b", ls=":", label="5% FAP (empirical noise)")
        ax.set_xlabel("TTV period (epochs)"); ax.set_ylabel("Delta chi2 (sinusoid vs linear)"); ax.legend(fontsize=8); ax.set_title(lab)
        u = pd.read_csv(f"{C.OUT}/t6_ttv_upper_limits_{lab}.csv")
        ax = axs[1, j]; ax.plot(u.P_ttv_epochs, u.A95_min, "k-", lw=0.8, label="95% UL"); ax.plot(u.P_ttv_epochs, u.A_obs_min, "g-", lw=0.5, label="fitted amplitude")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("TTV period (epochs)"); ax.set_ylabel("semi-amplitude (min)"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{C.OUT}/t7_ttv.png", dpi=90)
print("plots done")
