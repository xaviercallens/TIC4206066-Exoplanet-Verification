"""Step 19: (a) cycle-count uniqueness between sectors for both signals, from the S98-only period uncertainty and from the
per-event times; (b) audit figures: O-C diagrams, per-event depths, B null distribution, localisation summary,
per-sector folds.
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
from audit_common import AUD, save_json  # noqa

EPH = json.load(open(os.path.join(AUD, "out", "ephem_current.json")))
EV = json.load(open(os.path.join(AUD, "out", "s09_events.json")))
out = {}
# ---- (a) cycle counts
cand = {"A": json.load(open(os.path.join(AUD, "..", "..", "candidates", "tic4206066_p3.182945", "candidate.json"))),
        "B": json.load(open(os.path.join(AUD, "..", "..", "candidates", "tic4206066_p11.132976", "candidate.json")))}
for sig in ("A", "B"):
    lv = cand[sig]["leo_vetter"]
    P98, eP98 = lv["trap_per"], lv["trap_per_err"]
    T98 = lv["trap_epo"]
    res = {"S98_only_period_trapezoid": [P98, eP98], "S98_only_epoch": T98}
    for sec, tref in ((32, 2187.0), (6, 1479.0)):
        N = round((T98 - tref) / EPH[sig]["P"])
        spacing = EPH[sig]["P"] / N
        res[f"S{sec}"] = dict(cycles=int(N), alias_period_spacing_d=float(spacing), spacing_in_sigma_of_S98_period=float(spacing / eP98),
                              joint_minus_S98_period_sigma=float((EPH[sig]["P"] - P98) / eP98))
    # direct test: chi2 of the per-event times for cycle-count offsets dN = -2..2 between S98 and the archival sectors
    rows = [r for r in EV[sig]["pipe"]["rows"] if "dt_min" in r and r["coverage"] >= 0.7 and r["snr"] >= 2 and not r["timing_at_scan_edge"]]
    n = np.array([r["epoch"] for r in rows], float)
    tc = np.array([r["tc_pred"] + r["dt_min"] / 1440 for r in rows])
    et = np.array([max(r["e_dt_min"], 0.5) / 1440 for r in rows])
    arch = np.array([r["sector"] != 98 for r in rows])
    chi = {}
    for dN in (-2, -1, 0, 1, 2):
        nn = n + np.where(arch, dN, 0)
        A = np.vstack([np.ones_like(nn), nn]).T
        Wm = np.diag(1 / et ** 2)
        c = np.linalg.solve(A.T @ Wm @ A, A.T @ Wm @ tc)
        chi[str(dN)] = float(np.sum(((tc - A @ c) / et) ** 2))
    res["chi2_vs_cycle_offset"] = chi
    res["n_events_used"] = int(len(rows))
    out[sig] = res
    print(sig, json.dumps(res))
save_json(out, "out/s19_cycles.json")

# ---- (b) figures
os.makedirs(os.path.join(AUD, "plots"), exist_ok=True)
fig, axs = plt.subplots(2, 2, figsize=(15, 9))
for j, sig in enumerate(("A", "B")):
    for lab, col in (("ephem_pipe", "C0"), ("ephem_own", "C3")):
        e = EV[sig][lab]
        if not e:
            continue
        axs[0, j].errorbar(np.array(e["epochs"]) + (0.15 if lab == "ephem_own" else 0), e["oc_min"], e["err_min"], fmt="o", color=col, ms=4,
                           label=f"{lab.replace('ephem_', '')}: P={e['P']:.7f}, chi2={e['chi2']:.1f}/{e['dof']}")
    axs[0, j].axhline(0, color="k", lw=0.8)
    axs[0, j].set_xlabel("epoch")
    axs[0, j].set_ylabel("O-C (min)")
    axs[0, j].legend(fontsize=8)
    axs[0, j].set_title(f"Signal {sig}: O-C against the per-event linear ephemeris")
    for sname, mk in (("pipe", "o"), ("own", "s"), ("pipe_s98default", "^")):
        R = [r for r in EV[sig][sname]["rows"] if "depth_ppm" in r]
        axs[1, j].errorbar([r["epoch"] + {"pipe": 0, "own": 0.2, "pipe_s98default": 0.4}[sname] for r in R], [r["depth_ppm"] for r in R],
                           [r["depth_err_ppm"] for r in R], fmt=mk, ms=4, label=sname, alpha=0.8)
    axs[1, j].set_xlabel("epoch")
    axs[1, j].set_ylabel("depth (ppm)")
    axs[1, j].legend(fontsize=8)
    axs[1, j].set_title(f"Signal {sig}: every event, three photometries")
    axs[1, j].set_xscale("symlog", linthresh=20)
plt.tight_layout()
plt.savefig(os.path.join(AUD, "plots", "oc_and_event_depths.png"), dpi=100)
plt.close()

nb = json.load(open(os.path.join(AUD, "out", "s11_nullB.json")))
fig, ax = plt.subplots(1, 1, figsize=(9, 4))
for r in nb["results"][:5]:
    ax.axvline(r["obs_max_snr"], ls="--", lw=1)
    ax.text(r["obs_max_snr"], 0.9 - 0.12 * nb["results"].index(r), r["name"], fontsize=7, transform=ax.get_xaxis_transform())
ax.set_xlabel("max S/N over the dP scan")
ax.set_title("Signal B link: observed maxima (dashed); null percentiles in out/s11_nullB.json")
plt.tight_layout()
plt.savefig(os.path.join(AUD, "plots", "nullB_summary.png"), dpi=100)
plt.close()
print("figures written")
