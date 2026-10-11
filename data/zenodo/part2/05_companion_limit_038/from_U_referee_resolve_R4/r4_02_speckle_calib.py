"""R4 (b): empirical calibration of the Gaia DR3 'single-star' indicators against speckle-resolved companions.

Question: a companion of >=0.38 Msun at 0.2-2" (7-70 au at 35.6 pc) would have DeltaG <~ 2.2 (DeltaI <~ 1.7).
How often does such a companion leave a Gaia DR3 record that looks like TIC 4206066's
  (no other Gaia source within 11", ipd_frac_multi_peak = 0, ipd_gof_harmonic_amplitude <= 0.0222,
   RUWE <= 1.2143, corrected BP/RP flux excess C* consistent with zero, 5-parameter solution, not duplicated)?

Calibration sample: companions resolved by SOAR speckle imaging of TESS planet-candidate hosts
  Ziegler et al. 2020, AJ 159, 19 (VizieR J/AJ/159/19/table6) and Ziegler et al. 2021, AJ 162, 192 (J/AJ/162/192/table1),
  I-band separations and contrasts. Hosts are cross-matched to Gaia DR3 by position (TIC J2000 positions vs Gaia DR3
  J2000-epoch positions from VizieR I/355/gaiadr3), via the CDS xMatch service (Gaia TAP upload refused).
Comparison: the SOAR hosts with NO companion detected (to the SOAR sensitivity) -> base rate of 'target-like' records.

C* follows Riello et al. 2021 (A&A 649, A3) Table 2 and eq. 18 (checked in the arXiv source, aa39587.tex l.481-484, 860).
Run with audit_env python. Output: out/r4_02_speckle_calib.json, data/soar_gaia_match.csv
"""
import json, io, os, sys
import numpy as np, pandas as pd
from astropy.table import Table

H = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(H, "data"); O = os.path.join(H, "out")

TGT = dict(ruwe=1.2142984867095947, ipd_gofha=0.0222429521381855, ipd_fmp=0, G=10.51297378540039,
           bp_rp=1.7394800186157227, C=1.2975233793258667)


def cstar(C, x):
    x = np.asarray(x, float)
    f = np.where(x < 0.5, 1.154360 + 0.033772 * x + 0.032277 * x**2,
                 np.where(x < 4.0, 1.162004 + 0.011464 * x + 0.049255 * x**2 - 0.005879 * x**3,
                          1.057572 + 0.140537 * x))
    return np.asarray(C, float) - f


def sig_cstar(G):
    return 0.0059898 + 8.817481e-12 * np.asarray(G, float) ** 7.618399


def read_vizier_tsv(fn):
    lines = [l for l in open(fn) if not l.startswith("#") and l.strip()]
    hdr = lines[0].rstrip("\n").split("\t")
    rows = [l.rstrip("\n").split("\t") for l in lines[3:]]
    df = pd.DataFrame(rows, columns=hdr)
    return df.apply(lambda s: s.str.strip())


s2 = read_vizier_tsv(os.path.join(D, "J_AJ_162_192_table1.tsv")).rename(columns={"Con": "dm", "Sep": "sep"})
s1 = read_vizier_tsv(os.path.join(D, "J_AJ_159_19_table6.tsv")).rename(columns={"Dmag": "dm", "rho": "sep", "Year": "Obs"})
s1["src"] = "Z20"; s2["src"] = "Z21"
cols = ["src", "RAJ2000", "DEJ2000", "Comp", "Filt", "Obs", "sep", "dm"]
s = pd.concat([s1[cols], s2[cols]], ignore_index=True)
for c in ["RAJ2000", "DEJ2000", "sep", "dm", "Obs"]:
    s[c] = pd.to_numeric(s[c], errors="coerce")
s["hostkey"] = (s.RAJ2000.round(3).astype(str) + "_" + s.DEJ2000.round(3).astype(str))

# unique host positions
hosts = s.groupby("hostkey").agg(ra=("RAJ2000", "median"), dec=("DEJ2000", "median")).reset_index()
hosts["hid"] = np.arange(len(hosts))

cache = os.path.join(D, "soar_gaia_join.csv")
if not os.path.exists(cache):
    # Gaia TAP upload joins were refused (RemoteDisconnected, 2026-09-25); use the CDS xMatch service against the
    # VizieR copy of Gaia DR3 (I/355/gaiadr3), which carries RUWE, IPDgofha, IPDfmp, IPDfow, Dup, E(BP/RP), Solved
    # and J2000-epoch positions (RAJ2000/DEJ2000).
    from astroquery.xmatch import XMatch
    import astropy.units as u
    up = Table.from_pandas(hosts[["hid", "ra", "dec"]])
    r = XMatch.query(cat1=up, cat2="vizier:I/355/gaiadr3", max_distance=20 * u.arcsec, colRA1="ra", colDec1="dec").to_pandas()
    r = r.rename(columns={"Source": "source_id", "RAdeg": "ra_g", "DEdeg": "dec_g", "pmRA": "pmra", "pmDE": "pmdec",
                          "Plx": "parallax", "Gmag": "phot_g_mean_mag", "BP-RP": "bp_rp", "E(BP/RP)": "phot_bp_rp_excess_factor",
                          "RUWE": "ruwe", "IPDgofha": "ipd_gof_harmonic_amplitude", "IPDfmp": "ipd_frac_multi_peak",
                          "IPDfow": "ipd_frac_odd_win", "Dup": "duplicated_source", "Solved": "astrometric_params_solved",
                          "sepsi": "astrometric_excess_noise_sig"})
    keep = ["hid", "source_id", "ra_g", "dec_g", "RAJ2000", "DEJ2000", "pmra", "pmdec", "parallax", "phot_g_mean_mag", "bp_rp",
            "phot_bp_rp_excess_factor", "ruwe", "ipd_gof_harmonic_amplitude", "ipd_frac_multi_peak", "ipd_frac_odd_win",
            "duplicated_source", "astrometric_params_solved", "astrometric_excess_noise_sig", "angDist"]
    r[keep].to_csv(cache, index=False)
    with open(os.path.join(H, "logs", "r4_02_xmatch.log"), "a") as fh:
        fh.write(f"{pd.Timestamp.utcnow()} CDS xMatch vizier:I/355/gaiadr3, r=20 arcsec, {len(hosts)} hosts, {len(r)} rows\n")
g = pd.read_csv(cache)
g["duplicated_source"] = g["duplicated_source"].astype(str).str.lower().isin(["true", "1", "1.0"])

# identify the host: Gaia DR3 J2000-epoch position (VizieR RAJ2000/DEJ2000) closest to the TIC J2000 position
h = hosts.set_index("hid")
g["ra"] = g.ra_g; g["dec"] = g.dec_g
ra00 = g.RAJ2000.fillna(g.ra_g); de00 = g.DEJ2000.fillna(g.dec_g)
cosd = np.cos(np.radians(g.dec_g))
g["d00"] = 3600 * np.hypot((ra00 - h.loc[g.hid, "ra"].values) * cosd, de00 - h.loc[g.hid, "dec"].values)
cand = g[g.d00 < 1.5].sort_values(["hid", "phot_g_mean_mag"])
prim = cand.groupby("hid").head(1).set_index("hid")

rows = []
for hid, p in prim.iterrows():
    nb = g[(g.hid == hid) & (g.source_id != p.source_id)]
    dn = 3600 * np.hypot((nb.ra - p.ra) * np.cos(np.radians(p.dec)), nb.dec - p.dec)
    rows.append(dict(hid=hid, source_id=int(p.source_id), G=p.phot_g_mean_mag, bp_rp=p.bp_rp,
                     C=p.phot_bp_rp_excess_factor, ruwe=p.ruwe, gofha=p.ipd_gof_harmonic_amplitude,
                     fmp=p.ipd_frac_multi_peak, fow=p.ipd_frac_odd_win, dup=bool(p.duplicated_source),
                     aps=p.astrometric_params_solved,
                     nb_min_arcsec=float(dn.min()) if len(dn) else np.inf,
                     nb_n_lt3=int((dn < 3).sum()), nb_n_lt11=int((dn < 11).sum())))
P = pd.DataFrame(rows).set_index("hid")
P["Cstar"] = cstar(P.C, P.bp_rp)
P["sigC"] = sig_cstar(P.G)
P["zC"] = P.Cstar / P.sigC
hosts = hosts.join(P, on="hid")

# companions per host (I band only, unique by host+Comp, median over epochs)
c = s[s.sep.notna() & (s.Filt == "I")].groupby(["hostkey", "Comp"]).agg(sep=("sep", "median"), dm=("dm", "median"),
                                                                        n=("sep", "size")).reset_index()
det = s[s.sep.notna()].hostkey.unique()
tgt_C = float(cstar(TGT["C"], TGT["bp_rp"])); tgt_sig = float(sig_cstar(TGT["G"]))


def flags(r, c_tol_sigma=2.0):
    f = {"gaia_neighbour_lt3as": r.nb_n_lt3 > 0,
         "ipd_frac_multi_peak_gt0": r.fmp > TGT["ipd_fmp"],
         "ipd_gofha_gt_target": r.gofha > TGT["ipd_gofha"],
         "ruwe_gt_target": r.ruwe > TGT["ruwe"],
         "Cstar_gt_target_plus_2sig": r.Cstar > tgt_C + c_tol_sigma * r.sigC,
         "not_5p_or_dup": (r.aps != 31) or bool(r.dup)}
    f["any"] = any(f.values())
    return f


res = {"n_hosts_positions": int(len(hosts)), "n_hosts_matched": int(hosts.source_id.notna().sum()),
       "target": {**TGT, "Cstar": tgt_C, "sigma_Cstar_Riello_eq18": tgt_sig, "zC": tgt_C / tgt_sig}}

# comparison: hosts with no SOAR detection, G 8-13
nd = hosts[~hosts.hostkey.isin(det) & hosts.source_id.notna() & hosts.G.between(8, 13)]
fl = pd.DataFrame([flags(r) for _, r in nd.iterrows()])
res["no_companion_hosts_G8_13"] = {"n": int(len(nd)), "frac_flagged_by_indicator": fl.mean().round(3).to_dict(),
                                   "frac_target_like": float(1 - fl["any"].mean())}

# binaries: hosts whose ONLY SOAR-detected companion(s) are in the separation/contrast bin (I band)
cm = c.merge(hosts, on="hostkey")
cm = cm[cm.source_id.notna() & cm.G.between(8, 13)]
out_rows = []
for key, grp in cm.groupby("hostkey"):
    r0 = grp.iloc[0]
    f = flags(r0)
    out_rows.append(dict(hostkey=key, source_id=int(r0.source_id), G=r0.G, sep_min=grp.sep.min(),
                         sep_list=";".join(f"{a:.3f}" for a in grp.sep), dm_list=";".join(f"{a:.1f}" for a in grp.dm),
                         dm_of_closest=float(grp.sort_values("sep").dm.iloc[0]), ncomp=len(grp), ruwe=r0.ruwe,
                         gofha=r0.gofha, fmp=r0.fmp, Cstar=r0.Cstar, zC=r0.zC, nb_min=r0.nb_min_arcsec, **{f"f_{k}": v for k, v in f.items()}))
B = pd.DataFrame(out_rows)
B.to_csv(os.path.join(O, "r4_02_binaries_flags.csv"), index=False)
hosts.to_csv(os.path.join(D, "soar_gaia_match.csv"), index=False)

bins = [(0.0, 0.1), (0.1, 0.15), (0.15, 0.3), (0.3, 0.6), (0.6, 1.2), (1.2, 2.5), (0.2, 2.0), (0.15, 2.5)]
for dmmax in (1.76, 2.0, 2.5):
    for lo, hi in bins:
        sel = B[(B.ncomp == 1) & B.sep_min.between(lo, hi) & (B.dm_of_closest <= dmmax)]
        if len(sel) == 0:
            continue
        k = f"single_companion_sep{lo}-{hi}_dI_le{dmmax}"
        n = len(sel); nl = int((~sel.f_any).sum())
        # 95% upper limit on the miss rate (Clopper-Pearson / Wilson-like via beta quantile)
        from scipy.stats import beta
        ul = float(beta.ppf(0.95, nl + 1, n - nl)) if nl < n else 1.0
        res[k] = {"n": n, "n_target_like": nl, "frac_target_like": nl / n, "frac_target_like_95UL": ul,
                  "frac_flagged_by_indicator": {kk[2:]: float(sel[kk].mean()) for kk in sel.columns if kk.startswith("f_")},
                  "target_like_cases": sel[~sel.f_any][["source_id", "G", "sep_list", "dm_list", "ruwe", "gofha", "fmp", "Cstar", "zC", "nb_min"]].to_dict("records")}
# ---- robustness: looser, target-independent thresholds (the flag set is then less tuned to the target's values)
from scipy.stats import beta
VAR = {"strict_target_values": dict(ruwe=TGT["ruwe"], gofha=TGT["ipd_gofha"], fmp=0, csig=2.0),
       "loose_ruwe1.4_gofha0.05_fmp2_C3sig": dict(ruwe=1.4, gofha=0.05, fmp=2, csig=3.0),
       "Cstar_only_3sig": dict(ruwe=np.inf, gofha=np.inf, fmp=1e9, csig=3.0, nb=False, sol=False),
       "fmp_only_gt0": dict(ruwe=np.inf, gofha=np.inf, fmp=0, csig=1e9, nb=False, sol=False),
       "ruwe_only_gt1.4": dict(ruwe=1.4, gofha=np.inf, fmp=1e9, csig=1e9, nb=False, sol=False)}


def flag_any(r, v):
    f = (r.fmp > v["fmp"]) or (r.gofha > v["gofha"]) or (r.ruwe > v["ruwe"]) or (r.Cstar > tgt_C + v["csig"] * r.sigC)
    if v.get("nb", True):
        f = f or (r.nb_n_lt3 > 0)
    if v.get("sol", True):
        f = f or (r.aps != 31) or bool(r.dup)
    return bool(f)


rob = {}
Bh = cm.groupby("hostkey").agg(sep_min=("sep", "min"), n=("sep", "size"), dm=("dm", "first")).reset_index()
Bh = Bh.merge(hosts, on="hostkey")
for vn, v in VAR.items():
    base = np.array([flag_any(r, v) for _, r in nd.iterrows()])
    out = {"base_rate_target_like_no_companion_hosts": float(1 - base.mean()), "n_base": int(len(base))}
    for dmmax in (1.76, 2.0):
        for lo, hi in [(0.1, 0.2), (0.2, 2.0), (0.15, 1.2)]:
            sel = Bh[(Bh.n == 1) & Bh.sep_min.between(lo, hi) & (Bh.dm <= dmmax)]
            fl_b = np.array([flag_any(r, v) for _, r in sel.iterrows()])
            n = len(fl_b); k = int((~fl_b).sum())
            ul = float(beta.ppf(0.95, k + 1, n - k)) if k < n else 1.0
            out[f"sep{lo}-{hi}_dI_le{dmmax}"] = {"n": n, "n_target_like": k, "frac": k / n if n else None, "UL95": ul,
                                                 "LR_UL95_vs_base": ul / max(1e-9, float(1 - base.mean()))}
    rob[vn] = out
res["robustness_variants"] = rob
json.dump(res, open(os.path.join(O, "r4_02_speckle_calib.json"), "w"), indent=1, default=float)
print(json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk != "target_like_cases"}) for k, v in res.items()}, indent=1, default=float))
