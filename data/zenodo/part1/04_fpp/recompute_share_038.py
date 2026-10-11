"""Standalone recomputation of the FPP decomposition by companion mass (runs from this package; numpy-free, stdlib only).

Input: from_U_08_fpp/out_fpp_summary.csv (written by from_U_08_fpp/f04_summarize.py from the per-run files
from_U_08_fpp/runs/*/runs.json; each run's companion-mass bookkeeping is done by the hooks in f02_trice_run.py).
For every signal-A variant without a MOLUSC file:
  share_by_FPP_drop = 1 - FPP(companions <= 0.38 Msun) / FPP
  STP_wfrac_gt_0.38 = fraction of the STP (planet on unresolved bound companion) evidence from companions > 0.38 Msun
Reproduces from_T4_fpp_companion/t4_02_checks.json keys share_gt038_A_default and share_gt038_A_range
(printed in the Note as 'about 90% from those above 0.38 Msun' and 'FPP ~0.002-0.007').
"""
import csv, os

H = os.path.dirname(os.path.abspath(__file__))
rows = [r for r in csv.DictReader(open(os.path.join(H, "from_U_08_fpp/out_fpp_summary.csv")))
        if r["signal"] == "A" and not r["molusc"]]
print(f"{'variant':12s} {'n':>3s} {'FPP':>7s} {'FPP<=0.38':>9s} {'share':>6s} {'STPw>0.38':>9s}")
sh, sw, cut = [], [], []
for r in rows:
    fpp = float(r["FPP_mean"]); c = float(r["FPP_comp_mass_le_0.38_mean"]); w = float(r["STP_wfrac_mass_gt_0.38"])
    sh.append(1 - c / fpp); sw.append(w); cut.append(c)
    print(f"{r['variant']:12s} {r['n_runs']:>3s} {fpp:7.4f} {c:9.4f} {1 - c / fpp:6.3f} {w:9.3f}")
print(f"share by FPP drop: {min(sh):.3f}-{max(sh):.3f}; STP weight > 0.38: {min(sw):.3f}-{max(sw):.3f}; "
      f"FPP after cut: {min(cut):.4f}-{max(cut):.4f}; default FPP means: "
      f"{min(float(r['FPP_mean']) for r in rows):.4f}-{max(float(r['FPP_mean']) for r in rows):.4f}")
