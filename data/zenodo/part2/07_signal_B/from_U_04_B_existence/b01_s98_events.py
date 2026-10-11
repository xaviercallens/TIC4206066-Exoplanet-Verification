"""Step 1 + S98-only ephemeris.  Per-event table of predicted B events in Sector 98 (coverage, flags, background, dumps,
gaps; depth and S/N with A masked) and S98-only linear ephemeris with injection-calibrated timing errors."""
import numpy as np, json
import common_b as c

rng = np.random.default_rng(98)
q = c.load_qlp(98)
t, fdet, q_, bkg = q["t"], q["f"], q["q"], q["bkg"]
b29 = (q_ >> 29) & 1; b30 = (q_ >> 30) & 1; b2048 = (q_ >> 11) & 1; b32 = (q_ >> 5) & 1
clean = (q_ == 0) & np.isfinite(fdet)
relaxed = np.isfinite(fdet) & (b2048 == 0) & (b32 == 0) & (b30 == 0)      # keeps 2^29 (scattered-light half-orbits)
inA = c.in_A(t)
# own biweight detrending of SYS_RM_FLUX (A and B masked in the trend fit)
inB = np.zeros_like(t, bool)
for tc in c.pred_times(c.P_B98, c.T0_B98, t.min(), t.max())[1]:
    inB |= np.abs(t - tc) < c.T14_B
fsys_bw = np.full_like(t, np.nan)
for nm, m in [("clean", clean), ("relaxed", relaxed)]:
    pass
m_rel = relaxed & np.isfinite(q["sysrm"])
fsys_bw[m_rel] = c.detrend(t[m_rel], q["sysrm"][m_rel], "biweight", 0.75, mask=(inA | inB)[m_rel])

variants = {
    "QLP_DET_clean": (np.where(clean & ~inA, fdet, np.nan)),
    "QLP_DET_relaxed": (np.where(relaxed & ~inA, fdet, np.nan)),
    "SYSRM_bw075_relaxed": (np.where(relaxed & ~inA, fsys_bw, np.nan)),
    "SYSRM_bw075_clean": (np.where(clean & ~inA, fsys_bw, np.nan)),
}

# segments
gi = np.where(np.diff(t) > 0.1)[0]
seg_edges = [(t[0], t[gi[0]])] + [(t[gi[k] + 1], t[gi[k + 1]]) for k in range(len(gi) - 1)] + [(t[gi[-1] + 1], t[-1])]

def red_noise_scale(f, n_trials=2000):
    """Empirical scatter of fixed-centre depths at random epochs (avoiding A and B) divided by the median white error."""
    ok = np.isfinite(f)
    tt = t[ok]
    out_d, out_e = [], []
    tries = 0
    while len(out_d) < n_trials and tries < 20 * n_trials:
        tries += 1
        tc = rng.uniform(t.min(), t.max())
        if np.min(np.abs(tc - np.array(c.pred_times(c.P_B98, c.T0_B98, t.min(), t.max())[1]))) < 3 * c.T14_B:
            continue
        r = c.fit_event_fixed(t, f, tc)
        if r is None or r["n_in"] < 0.8 * (c.T14_B / (200 / 86400)) * 0.85:
            continue
        out_d.append(r["D"]); out_e.append(r["eD_white"])
    out_d = np.array(out_d); out_e = np.array(out_e)
    return dict(n=len(out_d), sd_emp=float(c.robust_std(out_d)), sd_std=float(np.std(out_d)), med_white=float(np.median(out_e)),
                beta=float(c.robust_std(out_d) / np.median(out_e)), mean=float(np.mean(out_d)))

ns, tcs = c.pred_times(c.P_B98, c.T0_B98, t.min(), t.max())
rows = []
rn = {k: red_noise_scale(v) for k, v in variants.items()}
print(json.dumps(rn, indent=1))
for n, tc in zip(ns, tcs):
    w = np.abs(t - tc) < c.T14_B / 2
    n_nom = c.T14_B / (200 / 86400.)
    row = dict(epoch=int(n), tc_pred=float(tc), n_cad_any=int(w.sum()), cov_any=float(w.sum() / n_nom),
               n_clean=int((w & clean).sum()), cov_clean=float((w & clean).sum() / n_nom),
               n_relaxed=int((w & relaxed).sum()),
               frac_b29=float(b29[w].mean()) if w.any() else None, frac_b30=float(b30[w].mean()) if w.any() else None,
               frac_b2048=float(b2048[w].mean()) if w.any() else None, frac_dump=float(b32[w].mean()) if w.any() else None,
               bkg_in_median=float(np.median(bkg[w])) if w.any() else None,
               bkg_clean_sector_median=float(np.median(bkg[clean])),
               dist_nearest_dump_h=float(np.min(np.abs(np.array(c.DUMPS[98]) - tc)) * 24),
               segment=[k for k, (a, b) in enumerate(seg_edges) if a - 0.01 <= tc <= b + 0.01],
               dist_seg_edge_h=float(min([min(abs(tc - a), abs(tc - b)) for (a, b) in seg_edges if a - 0.2 <= tc <= b + 0.2] or [np.nan]) * 24),
               A_overlap=bool(np.any(c.in_A(np.array([tc - c.T14_B / 2, tc, tc + c.T14_B / 2]), half=c.T14_A / 2 + c.T14_B / 2))))
    # background change in-transit relative to flanks (all cadences with finite bkg)
    fl = (np.abs(t - tc) > c.T14_B / 2) & (np.abs(t - tc) < 1.5 * c.T14_B)
    if w.any() and fl.any():
        row["bkg_in_minus_flank"] = float(np.median(bkg[w]) - np.median(bkg[fl]))
    for k, f in variants.items():
        r = c.fit_event_fixed(t, f, tc)
        if r is None:
            row[k] = None
            continue
        e = r["eD_white"] * max(rn[k]["beta"], 1.0)
        rr = c.fit_event(t, f, tc)
        row[k] = dict(D_ppm=r["D"] * 1e6, eD_white_ppm=r["eD_white"] * 1e6, eD_red_ppm=e * 1e6, snr_red=r["D"] / e, n_in=r["n_in"],
                      free_tc_dt_min=None if rr is None else rr["dt_min"], free_D_ppm=None if rr is None else rr["D"] * 1e6,
                      free_at_edge=None if rr is None else rr["at_edge"])
    rows.append(row)

for r in rows:
    print(r["epoch"], "%.3f" % r["tc_pred"], "cov_any %.2f clean %.2f" % (r["cov_any"], r["cov_clean"]), "b29", r["frac_b29"], "b2048", r["frac_b2048"], "b30", r["frac_b30"],
          "bkg", r["bkg_in_median"], "dump %.1fh" % r["dist_nearest_dump_h"], "seg", r["segment"], "edge %.1fh" % r["dist_seg_edge_h"])
    for k in variants:
        v = r[k]
        if v: print("    %-22s D=%6.0f +- %4.0f (white %4.0f) S/N %.1f n_in %d free dt %s" % (k, v["D_ppm"], v["eD_red_ppm"], v["eD_white_ppm"], v["snr_red"], v["n_in"], None if v["free_tc_dt_min"] is None else round(v["free_tc_dt_min"], 1)))

# ---- injection-calibrated timing error in clean S98 (QLP DET), trapezoid 1000 ppm, T14 1.6 h
f0 = variants["QLP_DET_clean"]
f_rel = variants["QLP_DET_relaxed"]
Bt = c.pred_times(c.P_B98, c.T0_B98, t.min(), t.max())[1]
def inject_timing(f, depth=1.0e-3, ntr=400, drop=0.0):
    dts, Ds = [], []
    k = 0
    while len(dts) < ntr and k < 50 * ntr:
        k += 1
        tc = rng.uniform(t.min(), t.max())
        if np.min(np.abs(tc - Bt)) < 5 * c.T14_B: continue
        w = np.abs(t - tc) < c.T14_B / 2
        if np.sum(np.isfinite(f[w])) < 0.9 * c.T14_B / (200 / 86400.): continue
        fi = f - depth * c.trapezoid(t - tc, c.T14_B, 0.15 * c.T14_B)
        if drop > 0:
            ww = np.abs(t - tc) < 0.5
            fi = fi.copy(); kill = ww & (rng.random(len(t)) < drop); fi[kill] = np.nan
        r = c.fit_event(t, fi, tc, scan_min=90)
        if r is None: continue
        dts.append(r["dt_min"]); Ds.append(r["D"])
    dts = np.array(dts); Ds = np.array(Ds)
    return dict(n=len(dts), sd_dt_min=float(np.std(dts)), rsd_dt_min=float(c.robust_std(dts)), p16_84=[float(x) for x in np.percentile(dts, [16, 84])],
                mean_dt=float(np.mean(dts)), sd_D_ppm=float(np.std(Ds) * 1e6), mean_D_ppm=float(np.mean(Ds) * 1e6), frac_gt30=float(np.mean(np.abs(dts) > 30)))
inj_full = inject_timing(f0)
inj_part = inject_timing(f0, drop=0.55)
inj_full_900 = inject_timing(f0, depth=0.9e-3)
print("injection timing full:", inj_full)
print("injection timing 55%-dropped (epoch-4-like):", inj_part)
print("injection timing full, 900 ppm:", inj_full_900)
# scattered-light epoch 2: inject into the relaxed data restricted to 2^29 segments
seg29 = relaxed & (b29 == 1)
f_sl = np.where(seg29 & ~inA, fdet, np.nan)
inj_sl = inject_timing(f_sl, ntr=200)
print("injection timing, 2^29 (scattered-light) data:", inj_sl)

def linfit(ep, tt, ee):
    ep = np.asarray(ep, float); tt = np.asarray(tt); ee = np.asarray(ee)
    A = np.vstack([np.ones_like(ep), ep]).T
    W = np.diag(1 / ee ** 2)
    cov = np.linalg.inv(A.T @ W @ A)
    p = cov @ A.T @ W @ tt
    res = tt - A @ p
    chi2 = float(np.sum((res / ee) ** 2))
    return dict(T0=float(p[0]), P=float(p[1]), eT0=float(np.sqrt(cov[0, 0])), eP=float(np.sqrt(cov[1, 1])), chi2=chi2, dof=len(tt) - 2,
                oc_min=(res * 1440).tolist(), epochs=ep.tolist(), times=tt.tolist(), errs_min=(ee * 1440).tolist())
# free-centre fits with the fixed-window fitter
ev = {}
for r in rows:
    n = r["epoch"]
    for k, f in [("clean", f0), ("relaxed", f_rel)]:
        rr = c.fit_event(t, f, r["tc_pred"], scan_min=90)
        if rr is not None:
            ev[(n, k)] = rr
            r["fit_" + k] = dict(dt_min=rr["dt_min"], D_ppm=rr["D"] * 1e6, eD_white_ppm=rr["eD_white"] * 1e6, n_in=rr["n_in"], at_edge=rr["at_edge"])
e_full = inj_full["sd_dt_min"] / 1440; e_part = inj_part["sd_dt_min"] / 1440; e_sl = inj_sl["sd_dt_min"] / 1440
def tt_of(n, k): return ev[(n, k)]["tc"]
sets = {
    "clean_0_1_4": ([0, 1, 4], [tt_of(0, "clean"), tt_of(1, "clean"), tt_of(4, "clean")], [e_full, e_full, e_part]),
    "relaxed_0_1_4(b30 kept)": ([0, 1, 4], [tt_of(0, "relaxed"), tt_of(1, "relaxed"), tt_of(4, "relaxed")], [e_full, e_full, e_full]),
    "clean_0_1_4 + SL_2": ([0, 1, 2, 4], [tt_of(0, "clean"), tt_of(1, "clean"), tt_of(2, "relaxed"), tt_of(4, "clean")], [e_full, e_full, e_sl, e_part]),
    "clean_0_1_only": ([0, 1], [tt_of(0, "clean"), tt_of(1, "clean")], [e_full, e_full]),
}
eph = {}
for k, (a, b, e) in sets.items():
    eph[k] = linfit(a, b, e)
    if eph[k]["dof"] == 0:
        eph[k]["P"] = float(b[1] - b[0]); eph[k]["eP"] = float(np.hypot(e[0], e[1]))
    print("S98-only ephemeris", k, "P = %.5f +- %.5f  T0 = %.4f +- %.4f chi2 %.2f dof %d O-C" % (eph[k]["P"], eph[k]["eP"], eph[k]["T0"], eph[k]["eT0"], eph[k]["chi2"], eph[k]["dof"]), np.round(eph[k]["oc_min"], 1))
json.dump(dict(rows=rows, red_noise=rn, injection_timing=dict(full=inj_full, partial55=inj_part, full900=inj_full_900, scattered_light=inj_sl),
               eph98=eph, segments=seg_edges,
               flags_note="QLP S98 QUALITY: 2048 predicted straylight (SPOC def.), 2^30 QLP internal poor quality (tess.mit.edu/qlp), 2^29 undocumented; set on 94/90/80/89% of the four a-half-orbit segments"),
          open("out/b01_s98_events.json", "w"), indent=1, default=float)
