"""Audit 8: plots from out_fpp_summary.json (audit_env)."""
import json, os
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); S = {r["variant"]: r for r in json.load(open(os.path.join(H, "out_fpp_summary.json")))}
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
for k, (sig, c) in enumerate((("A", "C0"), ("B", "C3"))):
    r = S.get(f"{sig}_base")
    if not r: continue
    mb = r["STP_mass_bins"]; labs = list(mb); ax[0].bar(np.arange(len(labs)) + 0.4 * k, [mb[l] for l in labs], 0.4, color=c, label=f"signal {sig}")
    eb = r["STP_ecc_bins"]; lab2 = list(eb); ax[1].bar(np.arange(len(lab2)) + 0.4 * k, [eb[l] for l in lab2], 0.4, color=c, label=f"signal {sig}")
ax[0].set_xticks(np.arange(len(labs)) + 0.2); ax[0].set_xticklabels(labs, rotation=45, fontsize=8); ax[0].axvline(3.7, ls="--", c="k"); ax[0].set_title("STP evidence share by companion (host) mass [Msun]\n(dashed: 0.38 Msun overluminosity limit)", fontsize=9); ax[0].legend()
ax[1].set_xticks(np.arange(len(lab2)) + 0.2); ax[1].set_xticklabels(lab2, rotation=45, fontsize=8); ax[1].set_title("STP evidence share by eccentricity", fontsize=9)
order = [v for v in S if not v.startswith("test")]
y = np.arange(len(order))
for i, v in enumerate(order):
    r = S[v]; c = "C0" if r["signal"] == "A" else "C3"
    ax[2].errorbar(r["FPP_mean"], i, xerr=r["FPP_sd"] or 0, fmt="o", c=c)
    if r.get("FPP_comp_mass_le_0.38_mean") is not None: ax[2].plot(r["FPP_comp_mass_le_0.38_mean"], i, "x", c=c)
ax[2].set_yticks(y); ax[2].set_yticklabels(order, fontsize=7); ax[2].set_xscale("log"); ax[2].axvline(0.015, ls="--", c="k")
ax[2].set_title("FPP per variant (o: TRICERATOPS; x: companions >0.38 Msun removed)", fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(H, "fig_fpp_audit08.png"), dpi=130); print("saved")
