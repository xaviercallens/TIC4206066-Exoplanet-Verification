"""Collect headline numbers (BASE), beta-scaled variant, grid ranges, per-event S98 errors; plot depths by config.
Output: out/s07_summary.json, plots/s07_depths_by_config.png"""
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
from b_lib import *
g = json.load(open(os.path.join(OUT, "s04_grid.json"))); G = {r["config"]: r for r in g}
beta = json.load(open(os.path.join(OUT, "s06_beta.json")))
injE = json.load(open(os.path.join(OUT, "s03_inj_BASE_event.json")))
b = G["BASE"]; S = {}
# beta-scaled (beta at 1.0 x T14, floor 1)
dep = {s: (b[f"D{s}"], b[f"eW{s}"] * max(beta[str(s)]["bins"]["1.0"]["beta"], 1)) for s in (6, 32, 98)}
sb = s98_vs_archival(dep)
S["BASE_beta"] = dict(dep=dep, **sb, p_common=float(stats.chi2.sf(sb["chi2_common_3sectors"], 2)))
# per-event S98 with per-event injection kappa
k98 = injE["sectors"]["98"]["kappa_event"]
S["S98_events"] = dict(epochs=b["events98"], cov=b["cov98"], D=b["Dev98"], e_white=b["eWev98"], e_inj=[e * k98 for e in b["eWev98"]], kappa_event=k98)
ns = {k: v["nsig_K"] for k, v in G.items()}
S["nsig_K_range"] = dict(min=min(ns.values()), max=max(ns.values()), median=float(np.median(list(ns.values()))),
                         argmin=min(ns, key=ns.get), argmax=max(ns, key=ns.get))
S["frac_configs_nsig_K_gt_2"] = float(np.mean([v > 2 for v in ns.values()]))
for s in (6, 32, 98):
    v = [G[k][f"D{s}"] for k in G]; S[f"D{s}_range"] = (min(v), max(v))
S["dchi2_K_range"] = (min(r["dchi2_K"] for r in g), max(r["dchi2_K"] for r in g))
jdump(S, "s07_summary.json"); print(json.dumps(S, indent=1, default=float)[:3000])
# plot
names = list(G.keys()); y = np.arange(len(names))
fig, ax = plt.subplots(1, 2, figsize=(14, 11), gridspec_kw=dict(width_ratios=[3, 1]))
for s, c, off in ((6, "C0", -0.25), (32, "C1", 0), (98, "C3", 0.25)):
    ax[0].errorbar([G[n][f"D{s}"] for n in names], y + off, xerr=[G[n][f"eK{s}"] for n in names], fmt="o", color=c, ms=3, label=f"S{s}")
ax[0].axvline(1099, color="C3", ls=":", lw=0.8); ax[0].axvline(838, color="C0", ls=":", lw=0.8); ax[0].axvline(1003, color="C1", ls=":", lw=0.8)
ax[0].set_yticks(y); ax[0].set_yticklabels(names, fontsize=7); ax[0].invert_yaxis(); ax[0].set_xlabel("B mid-transit depth (ppm), injection-scaled errors; dotted = manuscript")
ax[0].legend()
ax[1].barh(y, [G[n]["nsig_K"] for n in names], color="0.5"); ax[1].axvline(1.6, color="k", ls=":"); ax[1].set_yticks(y); ax[1].set_yticklabels([]); ax[1].invert_yaxis()
ax[1].set_xlabel("S98 - archival (sigma)")
plt.tight_layout(); plt.savefig(os.path.join(PLOTS, "s07_depths_by_config.png"), dpi=110)
