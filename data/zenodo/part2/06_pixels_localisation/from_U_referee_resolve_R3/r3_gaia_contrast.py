"""R3 resolution, step 3: what Gaia DR3 can and cannot exclude between 2.2" and ~14" (objection R3-v).

Inputs (read-only): 00_data Gaia DR3 cone (240"), 03_pixels neighbour census.
Contrast sensitivity: Brandeker & Cataldi 2019, A&A 621, A86, Table 1 (Gaia DR2; separation at which a given
detection fraction is reached for a given Delta G) -- transcribed from https://arxiv.org/html/1811.05301 .
EDR3 statement (El-Badry, Rix & Heintz 2021, MNRAS 506, 2269, Sec. 4.1): contrast limit Delta G ~ 4 at 1", ~7 at 2".
MOLUSC (Wood et al. 2021) combines Ziegler et al. 2018 (<4") and Brandeker & Cataldi 2019 (>4"); its table
(MOLUSC/code/gaia_contrast.txt, read from the local MOLUSC copy) is reported for comparison.
Depth-capable: a source of flux ratio f (TESS band) can give observed depth delta only if f >= delta/delta_max.
T-G colour transformation: Stassun et al. 2019 (TIC 8), eq. 1.
"""
import os, json
import numpy as np, pandas as pd

D = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
MOLUSC = "/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/a08_cache/MOLUSC/code/gaia_contrast.txt"

# Brandeker & Cataldi 2019 Table 1: rows Delta G; columns detection fraction 1,10,50,90,99 %; entries = separation (")
BC_DG = np.array([0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
BC_FR = np.array([0.01, 0.10, 0.50, 0.90, 0.99])
BC_SEP = np.array([
    [0.46, 0.50, 0.59, 0.65, 0.66], [0.52, 0.59, 0.68, 0.76, 0.78], [0.56, 0.62, 0.73, 0.88, 0.96],
    [0.59, 0.69, 0.94, 1.1, 1.1], [0.77, 0.90, 1.1, 1.2, 1.2], [0.88, 0.97, 1.1, 1.2, 1.2], [0.93, 1.1, 1.1, 1.2, 1.3],
    [0.93, 1.1, 1.1, 1.2, 1.3], [1.0, 1.1, 1.2, 1.3, 1.4], [1.2, 1.3, 1.4, 1.5, 1.5], [1.6, 1.7, 2.0, 2.2, 2.3],
    [1.9, 2.2, 2.6, 2.9, 3.0], [2.4, 2.8, 3.2, 3.5, 3.6], [3.4, 3.9, 4.3, 4.6, 4.8], [4.8, 5.1, 5.6, 6.0, 6.1],
    [5.7, 6.2, 7.0, 7.6, 7.8], [7.3, 7.8, 8.7, 9.3, 9.6], [8.6, 9.2, 10., 11., 11.], [9.6, 10., 12., 12., 12.]])


def sep_at(dG, frac):
    """separation (") at which detection fraction `frac` is reached for contrast dG (linear interp in the table)."""
    j = int(np.where(BC_FR == frac)[0][0])
    return float(np.interp(dG, BC_DG, BC_SEP[:, j]))


def pdet(dG, sep):
    """detection probability for contrast dG at separation sep: interpolate the table in separation per row,
    then in dG. 0 below the 1% separation, 1 beyond the 99% separation."""
    def row(i):
        s = BC_SEP[i]
        if sep <= s[0]: return 0.0
        if sep >= s[-1]: return 1.0
        return float(np.interp(sep, s, BC_FR))
    if dG <= BC_DG[0]: return row(0)
    if dG >= BC_DG[-1]: return row(len(BC_DG) - 1)
    i = np.searchsorted(BC_DG, dG) - 1; w = (dG - BC_DG[i]) / (BC_DG[i + 1] - BC_DG[i])
    return (1 - w) * row(i) + w * row(i + 1)


def g_minus_t(c):  # Stassun et al. 2019 eq. 1: T = G - 0.00522555 c^3 + 0.0891337 c^2 - 0.633923 c + 0.0324473
    return 0.00522555 * c ** 3 - 0.0891337 * c ** 2 + 0.633923 * c - 0.0324473


g = pd.read_csv(D + "/00_data/catalogs/gaiadr3_gaia_source_cone240as.csv")
tgt = g[g.source_id == 3220388198192519424].iloc[0]
Gt, ct = float(tgt.phot_g_mean_mag), float(tgt.bp_rp)
Tt = Gt - g_minus_t(ct)
out = dict(target=dict(G=Gt, bp_rp=ct, T_from_G=Tt, ruwe=float(tgt.ruwe), ipd_frac_multi_peak=int(tgt.ipd_frac_multi_peak),
                       ipd_gof_harmonic_amplitude=float(tgt.ipd_gof_harmonic_amplitude), duplicated_source=bool(tgt.duplicated_source)))
oth = g[g.source_id != tgt.source_id].copy()
oth["sep"] = oth.dist_arcsec_j2016pos
near = oth[oth.sep < 15].sort_values("sep")
out["gaia_sources_within_15as"] = near[["source_id", "sep", "phot_g_mean_mag", "bp_rp"]].assign(dG=near.phot_g_mean_mag - Gt).to_dict("records")

# depth-capable contrast limits
res = []
for depth_ppm in (400, 500, 575):
    for dmax in (1.0, 0.5):
        dT = -2.5 * np.log10(depth_ppm * 1e-6 / dmax)
        for c in (0.5, 1.0, 1.74, 2.5, 3.0):
            dG = dT + g_minus_t(c) - g_minus_t(ct)
            res.append(dict(depth_ppm=depth_ppm, max_eclipse=dmax, bp_rp_blend=c, dT_max=dT, dG_max=dG,
                            sep50=sep_at(min(dG, 14), 0.50), sep90=sep_at(min(dG, 14), 0.90), sep99=sep_at(min(dG, 14), 0.99)))
out["depth_capable_limits"] = res

# expected number of depth-capable, Gaia-missed chance-aligned stars in annuli, using the local Gaia G distribution
# (sources 30-240" from the target; assume the same luminosity function behind the target)
bg = oth[(oth.sep > 30) & np.isfinite(oth.phot_g_mean_mag)]
area_bg = np.pi * (240 ** 2 - 30 ** 2)
exp = {}
for depth_ppm, dmax in ((500, 1.0), (500, 0.5), (400, 1.0)):
    dT = -2.5 * np.log10(depth_ppm * 1e-6 / dmax)
    for (r0, r1) in ((0.0, 2.2), (2.2, 4.0), (2.2, 11.0), (4.0, 11.0), (11.0, 14.0)):
        rr = np.linspace(r0, r1, 400); dr = rr[1] - rr[0]
        n_all = n_cap = n_miss = 0.0
        for _, s in bg.iterrows():
            c = s.bp_rp if np.isfinite(s.bp_rp) else 1.0
            dTs = (s.phot_g_mean_mag - g_minus_t(c)) - Tt
            dG = s.phot_g_mean_mag - Gt
            ring = 2 * np.pi * rr * dr / area_bg                 # expected number of such stars per ring
            n_all += ring.sum()
            if dTs <= dT:
                n_cap += ring.sum()
                n_miss += sum(ring[k] * (1 - pdet(dG, rr[k])) for k in range(len(rr)))
        exp[f"depth{depth_ppm}_max{dmax}_{r0}-{r1}as"] = dict(n_all_gaia_stars=n_all, n_depth_capable=n_cap,
                                                                n_depth_capable_missed_by_gaia=n_miss)
out["expected_chance_aligned"] = exp
out["n_background_sources_used"] = int(len(bg)); out["bg_area_arcsec2"] = area_bg

# MOLUSC table (Ziegler 2018 <3.1", Brandeker & Cataldi >4.2") for comparison
try:
    m = pd.read_csv(MOLUSC, sep=r"\s+")
    out["molusc_gaia_contrast_rows"] = m[(m.Sep >= 2000) & (m.Sep <= 11500)].to_dict("records")
except Exception as e:
    out["molusc_gaia_contrast_rows"] = str(e)

json.dump(out, open(os.path.join(OUT, "r3_gaia_contrast.json"), "w"), indent=1, default=float)
print("target G %.2f BP-RP %.2f T(from G) %.2f" % (Gt, ct, Tt))
print("Gaia sources <15\":", [(round(r['sep'], 1), round(r['phot_g_mean_mag'], 2), round(r['dG'], 2)) for r in out["gaia_sources_within_15as"]])
for r in res:
    if r["bp_rp_blend"] in (0.5, 1.74, 3.0):
        print("depth %d max %.1f c=%.2f dT<=%.2f dG<=%.2f sep50=%.1f sep90=%.1f sep99=%.1f" % (
            r["depth_ppm"], r["max_eclipse"], r["bp_rp_blend"], r["dT_max"], r["dG_max"], r["sep50"], r["sep90"], r["sep99"]))
for k, v in exp.items():
    print(k, {a: round(b, 4) for a, b in v.items()})
