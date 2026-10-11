"""R4 (b): where do the >0.38 Msun companions that survive the MOLUSC-mode Gaia limits sit on the sky, and what is
left of the MOLUSC-mode STP evidence once each separation regime gets its own Gaia test?

Inputs (read-only, reused): auditor's fresh MOLUSC b26efbf survivor file (fresh5_kept.csv, sha256 63cdae36...; columns
include 'Projected Separation(AU)' at epoch 2016.0 and 'DeltaG'); verifier's fixed-companion-mass TRICERATOPS likelihood
grid for signal A (08_fpp/verify/runs_v/grid_A_v_A_S32+S98_5min.json, N=4e5 per mass) and the verifier's default runs
(same fold) for the other scenarios' evidences. The combination logic is copied from 08_fpp/verify/v06_combine.py
(MOLUSC mode as coded: survivors with a(1-e) > 10 AU, masses clamped at 0.1 Msun, padded to N = 1e6, prior 1).
Regimes (angular projected separation rho at 35.6 pc):
  rho < 0.15"      : Gaia-blended -> XP overluminosity limit (M > 0.38 excluded at 2 sigma)
  0.15" <= rho < 2": partially resolved -> C* / ipd_frac_multi_peak (speckle calibration r4_02: 0/42 target-like at
                     0.15-1.2", 0/38 at 0.2-2" for dI <= 2)
  rho >= 2"        : resolved by Gaia (no source within 11")
audit_env python (no TRICERATOPS import needed).
"""
import hashlib, json, os
import numpy as np, pandas as pd

H = os.path.dirname(os.path.abspath(__file__))
AUD = os.path.abspath(os.path.join(H, "..", ".."))
V = os.path.join(AUD, "08_fpp/verify/runs_v")
F = "/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/a08_cache/molusc_run/fresh5_kept.csv"
M, N, DIST = 0.606, 1_000_000, 1000 / 28.075
out = {"molusc_file_sha256": hashlib.sha256(open(F, "rb").read()).hexdigest()}

d = pd.read_csv(F)
aa, ee = d["semi-major axis(AU)"].to_numpy(), d["eccentricity"].to_numpy()
m_mol = np.maximum(d["mass ratio"].to_numpy(), 0.1 / M) * M
peri = aa * (1 - ee) > 10
rho = d["Projected Separation(AU)"].to_numpy() / DIST       # arcsec
dG = d["DeltaG"].to_numpy()

res_sig = {}
for sig in ("A", "B"):
    g = json.load(open(os.path.join(V, f"grid_{sig}_v_{sig}_S32+S98_5min.json")))
    rows = sorted(g["rows"], key=lambda r: r["m"])
    mg = np.array([r["m"] for r in rows])
    lS = np.array([r["lnL_STP"] for r in rows]); lP = np.array([r["lnL_PTP"] for r in rows])
    ref = max(lS[np.isfinite(lS)].max(), lP.max())
    EL = lambda ms, lg: np.exp(np.interp(ms, mg, np.where(np.isfinite(lg), lg, -1e4)) - ref)
    runs = json.load(open(os.path.join(V, f"default_{sig}_v_{sig}_S32+S98_5min.json")))["runs"]
    base = np.log(np.mean(np.exp(np.array([r["lnZ"] for r in runs]) - ref), axis=0)) + ref

    def fpp(sel, wextra=None):
        w = np.ones(sel.sum()) if wextra is None else wextra[sel]
        zS = np.log((EL(m_mol[sel], lS) * w).sum() / N) + ref
        zP = np.log((EL(m_mol[sel], lP) * w).sum() / N) + ref
        z = base.copy(); z[6] = zS; z[3] = zP
        fin = np.isfinite(z); p = np.zeros_like(z); p[fin] = np.exp(z[fin] - z[fin].max()); p /= p.sum()
        return float(1 - p[0] - p[3] - p[9])

    wS = EL(m_mol, lS)
    k = peri
    tot = wS[k].sum()
    big = m_mol > 0.38
    reg = {"blended_lt0.15": rho < 0.15, "partial_0.15-2": (rho >= 0.15) & (rho < 2.0), "resolved_ge2": rho >= 2.0}
    share = {f"M>0.38_{rn}": float(wS[k & big & rm].sum() / tot) for rn, rm in reg.items()}
    share["M<=0.38_all"] = float(wS[k & ~big].sum() / tot)
    share["M>0.38_total"] = float(wS[k & big].sum() / tot)
    r = {"STP_evidence_share_molusc_kept": share,
         "rho_quantiles_arcsec_M>0.38_kept": [float(x) for x in np.quantile(rho[k & big], [0.05, 0.25, 0.5, 0.75, 0.95])],
         "dG_range_M>0.38_kept": [float(np.nanmin(dG[k & big])), float(np.nanmax(dG[k & big]))],
         "FPP_molusc_as_coded": fpp(k),
         "FPP_molusc_minus_M>0.38_blended_only": fpp(k & ~(big & reg["blended_lt0.15"])),
         "FPP_molusc_minus_M>0.38_all_regimes": fpp(k & ~big)}
    # conservative: keep the 0.15-2" massive companions with weight = 95% UL of the speckle-calibrated miss rate / base rate
    for nm, ul in (("UL95_0.076", 0.0758), ("UL95_0.089_dI1.76", 0.0894)):
        wx = np.ones(len(m_mol)); wx[big & reg["partial_0.15-2"]] = ul / 0.3787   # likelihood ratio vs target-like base rate
        keep = k & ~(big & (reg["blended_lt0.15"] | reg["resolved_ge2"]))
        r[f"FPP_molusc_partial_regime_LR_{nm}"] = fpp(keep, wx)
    res_sig[sig] = r
out["signals"] = res_sig
out["n_kept_peri_gt10"] = int(peri.sum()); out["n_kept_peri_gt10_M_gt_0.38"] = int((peri & (m_mol > 0.38)).sum())
json.dump(out, open(os.path.join(H, "out", "r4_03_molusc_regimes.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
