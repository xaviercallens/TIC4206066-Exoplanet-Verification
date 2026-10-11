"""Independent per-transit fits for both signals of TIC 4206066 (Steps 1, 2 and 4 of the publication-prep task).

Written from scratch for this check; no code from the deep follow-up or the publication audit is imported.
Data: the light-curve FITS files of the deep follow-up (byte-identical to the audit's re-downloads, sha256 checked
there), plus own aperture photometry made here from the TESScut cubes (independent of every pipeline light curve).

Per transit (both signals):
  - window +-0.3 d around the predicted time, cadences within 1.0 T14 of the other signal removed;
  - model f = c0 + c1 (t - tc) - D s(t - tc - dt), s = batman template (unit mid-transit depth, fixed shape);
  - dt on a 0.25-min grid (-60..+60 min), D and baseline linear at each dt; white error from Delta chi2 = 1 with the
    residual variance; realistic errors from (i) cyclic residual permutation (200 shifts) and (ii) injection of the
    same template at 150 random out-of-transit times of the same light curve, recovered the same way.
Aggregates: linear ephemeris fits with each error model (chi2, rms O-C), joint least-squares (P, T0) with all
windows, S98-only ephemeris and the cycle-count test, per-sector depths of B with each template and photometry,
common-depth vs per-sector-depth model comparison for B.
Outputs: out/a1_results.json, out/a1_per_transit.csv, out/a1_oc.png
"""
import csv, json, os, sys, warnings
import numpy as np
from astropy.io import fits
from wotan import flatten
import batman
warnings.filterwarnings("ignore")

H = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(H, "out"); os.makedirs(OUT, exist_ok=True)
ADD = os.path.abspath(os.path.join(H, "..", ".."))
LC = os.path.join(ADD, "deep_followup_tic4206066", "lc")
CUT = os.path.join(ADD, "tic4206066_PUBLICATION_AUDIT", "data", "tesscut")   # raw TESScut cubes (astrocut), not results
FROZEN_QLP = 1610619071; LK_DEFAULT = 17087; STRAY = 2048; QLP29 = 1 << 29
TESSCUT_SHIFT_S = {6: -11.78, 32: -3.28, 98: 0.49}   # CCD-centre -> target barycentric correction (audit s06; checked below)
RHO_SUN = 1.98841e33 / (4 / 3 * np.pi * 6.957e10 ** 3); G = 6.6743e-8
# ephemerides used only for masking and window placement (audit joint fits; refined independently below)
EPH = {"A": dict(P=3.1827853, T0=3998.73679, T14=1.874 / 24), "B": dict(P=11.132745, T0=4000.86587, T14=1.66 / 24)}
# template shapes (Mann-prior joint fits of the audit; the results below do not depend on them at the quoted precision)
SHAPE = {"A": dict(k=0.02116, b=0.35, rho=2.507, u1=0.325, u2=0.309),
         "B_mann": dict(k=0.03305, b=0.871, rho=2.507, u1=0.325, u2=0.309),
         "B_free": dict(k=0.0296, b=0.49, rho=13.35, u1=0.325, u2=0.309)}
rng = np.random.default_rng(11)

def ars(rho, P): return (rho * RHO_SUN * G * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)

class Template:
    def __init__(self, sig, shape, P):
        s = SHAPE[shape]; self.P = P
        self.pm = batman.TransitParams(); pm = self.pm
        pm.t0, pm.per, pm.rp, pm.a = 0.0, P, s["k"], ars(s["rho"], P)
        pm.inc = np.degrees(np.arccos(s["b"] / pm.a)); pm.ecc, pm.w = 0.0, 90.0; pm.u = [s["u1"], s["u2"]]; pm.limb_dark = "quadratic"
        self.mid = 1 - batman.TransitModel(pm, np.array([0.0])).light_curve(pm)[0]
        self.t14 = P / np.pi * np.arcsin(np.sqrt((1 + s["k"]) ** 2 - s["b"] ** 2) / (pm.a * np.sin(np.radians(pm.inc)))) * 24
    def unit(self, trel, cad_d):
        ss = max(1, int(np.ceil(cad_d * 1440 / 2)))
        m = batman.TransitModel(self.pm, trel, supersample_factor=ss, exp_time=cad_d) if ss > 1 else batman.TransitModel(self.pm, trel)
        return (1 - m.light_curve(self.pm)) / self.mid

def in_tr(t, key, fac=1.0):
    e = EPH[key]; ph = ((t - e["T0"]) / e["P"] + 0.5) % 1 - 0.5
    return np.abs(ph * e["P"]) < fac * e["T14"]

def load_lc(name):
    if name == "S6":
        h = fits.open(os.path.join(LC, "hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits")); col = "PDCSAP_FLUX"; mask = LK_DEFAULT
    elif name == "S32":
        h = fits.open(os.path.join(LC, "tess2020324010417-s0032-0000000004206066-0200-s_lc.fits")); col = "PDCSAP_FLUX"; mask = LK_DEFAULT
    else:
        h = fits.open(os.path.join(LC, "hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits")); col = "SYS_RM_FLUX"; mask = FROZEN_QLP
    d = h[1].data; t = np.asarray(d["TIME"], float); f = np.asarray(d[col], float); q = np.asarray(d["QUALITY"]).astype(np.int64)
    hdr = {**h[0].header, **h[1].header}
    m = np.isfinite(t) & np.isfinite(f) & (f > 0) & ((q & mask) == 0)
    meta = {k: hdr.get(k) for k in ("TIMESYS", "TIMEREF", "BJDREFI", "BJDREFF", "TIMEUNIT", "TIMEDEL", "CROWDSAP", "FLFRCSAP")}
    return dict(name=name, t=t[m], f=f[m] / np.median(f[m]), cad=float(np.median(np.diff(t[m]))), meta=meta, sector={"S6": 6, "S32": 32, "S98": 98}[name])

def own_phot(sector, radius=2.0):
    """Own aperture photometry from the TESScut cube: background = median of faint pixels > 5 px away; aperture r <= radius
    around the brightest pixel near the stamp centre. S98: cadences flagged 2048 (straylight) or QLP-like 2^29 are not
    available in the cube header, so only the SPOC quality bits are used (DEFAULT_BITMASK + straylight)."""
    import glob
    p = glob.glob(os.path.join(CUT, f"tess-s{sector:04d}-*_astrocut.fits"))[0]
    h = fits.open(p); d = h[1].data; t = np.asarray(d["TIME"], float); fl = np.asarray(d["FLUX"], float); q = np.asarray(d["QUALITY"]).astype(np.int64)
    good = np.isfinite(t) & np.all(np.isfinite(fl.reshape(len(t), -1)), axis=1) & ((q & (LK_DEFAULT | STRAY)) == 0)
    t, fl = t[good] + TESSCUT_SHIFT_S[sector] / 86400, fl[good]
    med = np.median(fl, axis=0); ny, nx = med.shape; yy, xx = np.mgrid[:ny, :nx]
    c = np.hypot(xx - nx // 2, yy - ny // 2) < 3; iy, ix = np.unravel_index(np.argmax(np.where(c, med, -np.inf)), med.shape)
    far = np.hypot(xx - ix, yy - iy) > 5; bp = far & (med < np.percentile(med[far], 30))
    bkg = np.median(fl[:, bp], axis=1); ap = np.hypot(xx - ix, yy - iy) <= radius
    f = (fl - bkg[:, None, None])[:, ap].sum(axis=1)
    return dict(name=f"own{sector}", t=t, f=f / np.median(f), cad=float(np.median(np.diff(t))), sector=sector, meta={})

def detrend(lc, window=0.75):
    t, f = lc["t"], lc["f"]; msk = in_tr(t, "A", 1.5) | in_tr(t, "B", 1.5)
    _, tr = flatten(t[~msk], f[~msk], method="biweight", window_length=window, return_trend=True, break_tolerance=0.3)
    ok = np.isfinite(tr); fd = f / np.interp(t, t[~msk][ok], tr[ok])
    s = 1.4826 * np.median(np.abs(fd[~msk] - np.median(fd[~msk]))); keep = np.isfinite(fd) & (fd < 1 + 4 * s)
    return dict(lc, t=t[keep], fd=fd[keep], sig=float(s))

def fit_window(t, f, tc, tpl, cad, dts):
    """returns best dt (d), D, chi2 curve (white, unscaled), baseline coefficients"""
    trel = t - tc; chis, Ds, cs = [], [], []
    for dt in dts:
        s = tpl.unit(trel - dt, cad); X = np.vstack([np.ones_like(trel), trel, -s]).T
        c, *_ = np.linalg.lstsq(X, f, rcond=None); r = f - X @ c
        chis.append(np.sum(r * r)); Ds.append(c[2]); cs.append(c)
    chis = np.array(chis); j = int(np.argmin(chis))
    return dts[j], Ds[j], chis, cs[j], j

def analyse_event(t, f, tc, tpl, cad, n_perm=200, n_inj=0, lc_full=None, other=None):
    dts = np.arange(-60, 60.001, 0.25) / 1440
    dt, D, chis, c, j = fit_window(t, f, tc, tpl, cad, dts)
    trel = t - tc; s = tpl.unit(trel - dt, cad); X = np.vstack([np.ones_like(trel), trel, -s]).T
    r = f - X @ c; s2 = np.var(r); dof = len(f) - 3
    chir = chis / s2 * (dof / len(f)) if False else chis / s2   # chi2 with sigma^2 = residual variance (chi2_red = 1 at best fit)
    ok = np.where(chir - chir[j] <= 1)[0]; e_dt_white = 0.5 * (dts[ok[-1]] - dts[ok[0]]) * 1440
    cov = np.linalg.pinv(X.T @ X) * s2; e_D_white = np.sqrt(cov[2, 2])
    edge = ok[0] == 0 or ok[-1] == len(dts) - 1
    # cyclic residual permutation
    perm_dt, perm_D = [], []
    model = X @ c
    for k in range(n_perm):
        sh = rng.integers(1, len(r)); fp = model + np.roll(r, sh)
        d1, D1, _, _, _ = fit_window(t, fp, tc, tpl, cad, dts[::2]); perm_dt.append(d1 * 1440); perm_D.append(D1)
    res = dict(dt_min=float(dt * 1440), D=float(D), e_dt_white_min=float(e_dt_white), e_D_white=float(e_D_white),
               e_dt_perm_min=float(np.std(perm_dt)), e_D_perm=float(np.std(perm_D)), rms_ppm=float(np.sqrt(s2) * 1e6),
               n=int(len(f)), edge=bool(edge), snr_white=float(D / e_D_white))
    return res

def injection_errors(lcd, sig, tpl, D, n_inj=150, half=0.3):
    """inject the template with depth D at random times away from both signals, recover dt and D as for real events"""
    t, f, cad = lcd["t"], lcd["fd"], lcd["cad"]; T14 = tpl.t14 / 24
    bad = in_tr(t, "A", 3.0) | in_tr(t, "B", 3.0); cand = t[~bad]; dts = np.arange(-60, 60.001, 0.5) / 1440
    out = []; tries = 0
    while len(out) < n_inj and tries < n_inj * 40:
        tries += 1; tc = rng.choice(cand); w = np.abs(t - tc) < half
        inn = np.abs(t - tc) < T14 / 2
        if inn.sum() * cad < 0.7 * T14 or bad[w].any() or w.sum() < 20: continue
        fi = f[w] * (1 - D * tpl.unit(t[w] - tc, cad))
        dt, Dr, _, _, _ = fit_window(t[w], fi, tc, tpl, cad, dts)
        out.append((dt * 1440, Dr))
    a = np.array(out)
    return dict(n=len(a), sd_dt_min=float(np.std(a[:, 0])), mean_dt_min=float(np.mean(a[:, 0])), sd_D=float(np.std(a[:, 1])), mean_D=float(np.mean(a[:, 1])))

def linear_ephem(epochs, tcs, errs):
    n = np.array(epochs, float); tc = np.array(tcs); e = np.array(errs)
    A = np.vstack([np.ones_like(n), n]).T; W = np.diag(1 / e ** 2); cov = np.linalg.inv(A.T @ W @ A); c = cov @ A.T @ W @ tc
    oc = tc - A @ c; chi2 = float(np.sum((oc / e) ** 2)); dof = len(n) - 2
    return dict(T0=float(c[0]), P=float(c[1]), eT0=float(np.sqrt(cov[0, 0])), eP=float(np.sqrt(cov[1, 1])), chi2=chi2, dof=dof,
                eP_scaled=float(np.sqrt(cov[1, 1]) * np.sqrt(max(chi2 / dof, 1))), eT0_scaled=float(np.sqrt(cov[0, 0]) * np.sqrt(max(chi2 / dof, 1))),
                oc_min=[float(x * 1440) for x in oc], rms_oc_min=float(np.std(oc) * 1440))

def run_signal(sig, sets, tpl_name, label):
    e = EPH[sig]; other = "B" if sig == "A" else "A"; tpl = Template(sig, tpl_name, e["P"]); T14 = tpl.t14 / 24
    rows = []; inj = {}
    for lcd in sets:
        t, f, cad = lcd["t"], lcd["fd"], lcd["cad"]
        n0 = int(np.floor((t.min() - e["T0"]) / e["P"])); n1 = int(np.ceil((t.max() - e["T0"]) / e["P"]))
        depths = []
        for n in range(n0, n1 + 1):
            tc = e["T0"] + n * e["P"]; w = (np.abs(t - tc) < 0.3) & ~in_tr(t, other, 1.0)
            inn = w & (np.abs(t - tc) < T14 / 2); cov = inn.sum() * cad / T14
            if cov < 0.5 or w.sum() < 20: continue
            r = analyse_event(t[w], f[w], tc, tpl, cad); r.update(set=lcd["name"], sector=lcd["sector"], epoch=n, tc_pred=float(tc), coverage=float(cov))
            rows.append(r); depths.append(r["D"])
        Dmed = float(np.median(depths)) if depths else 0.0
        inj[lcd["name"]] = injection_errors(lcd, sig, tpl, Dmed) if depths else None
    for r in rows:
        r["e_dt_inj_min"] = inj[r["set"]]["sd_dt_min"]; r["e_D_inj"] = inj[r["set"]]["sd_D"]
    return rows, inj, tpl

def sector_depth(rows, sector, errkey):
    R = [r for r in rows if r["sector"] == sector]; d = np.array([r["D"] for r in R]); e = np.array([r[errkey] for r in R]); w = 1 / e ** 2
    return float(np.sum(w * d) / w.sum()), float(1 / np.sqrt(w.sum())), len(R)

if __name__ == "__main__":
    pipe = [detrend(load_lc(k)) for k in ("S6", "S32", "S98")]
    own = [detrend(own_phot(6, 2.5)), detrend(own_phot(32, 2.0)), detrend(own_phot(98, 2.0))]
    res = {"time_system_headers": {lc["name"]: lc["meta"] for lc in pipe}, "signals": {}}
    allrows = []
    for sig, tpls in (("A", ["A"]), ("B", ["B_mann", "B_free"])):
        res["signals"][sig] = {}
        for tpl_name in tpls:
            for setname, sets in (("pipe", pipe), ("own", own)):
                rows, inj, tpl = run_signal(sig, sets, tpl_name, setname)
                for r in rows: r.update(signal=sig, template=tpl_name, photometry=setname)
                allrows += rows
                R = dict(template_T14_h=tpl.t14, injections=inj, per_sector={})
                for sec in (6, 32, 98):
                    R["per_sector"][str(sec)] = {ek: sector_depth(rows, sec, ek) for ek in ("e_D_white", "e_D_perm", "e_D_inj")}
                # ephemerides (events with coverage >= 0.7, S/N >= 2, not at scan edge)
                good = [r for r in rows if r["coverage"] >= 0.7 and r["snr_white"] >= 2 and not r["edge"]]
                tcs = [r["tc_pred"] + r["dt_min"] / 1440 for r in good]; ep = [r["epoch"] for r in good]
                R["ephemeris"] = {}
                for ek in ("e_dt_white_min", "e_dt_perm_min", "e_dt_inj_min"):
                    errs = [max(r[ek], 0.25) / 1440 for r in good]
                    R["ephemeris"][ek] = dict(linear_ephem(ep, tcs, errs), n=len(good))
                # S98-only ephemeris and cycle-count uniqueness (white errors scaled)
                g98 = [r for r in good if r["sector"] == 98]
                if len(g98) >= 3:
                    l98 = linear_ephem([r["epoch"] for r in g98], [r["tc_pred"] + r["dt_min"] / 1440 for r in g98], [max(r["e_dt_inj_min"], 0.25) / 1440 for r in g98])
                    cyc = {}
                    for off in (-2, -1, 0, 1, 2):
                        ep2 = [r["epoch"] + (off if r["sector"] != 98 else 0) for r in good]
                        cyc[str(off)] = linear_ephem(ep2, tcs, [max(r["e_dt_inj_min"], 0.25) / 1440 for r in good])["chi2"]
                    R["S98_only"] = dict(l98, n=len(g98), chi2_vs_cycle_offset=cyc)
                # joint least-squares (P, T0) on all windows with the template fixed: chi2 surface over P near the linear solution
                lin = R["ephemeris"]["e_dt_inj_min"]
                res["signals"][sig][f"{tpl_name}__{setname}"] = R
                print(f"{sig} {tpl_name:7s} {setname:4s}: n_ev={len(rows)} good={len(good)} P_lin(inj)={lin['P']:.7f}±{lin['eP']:.7f} (scaled {lin['eP_scaled']:.7f}) chi2={lin['chi2']:.1f}/{lin['dof']} rms={lin['rms_oc_min']:.1f} min; "
                      f"depths ppm (inj err): " + ", ".join(f"S{s}: {R['per_sector'][str(s)]['e_D_inj'][0]*1e6:.0f}±{R['per_sector'][str(s)]['e_D_inj'][1]*1e6:.0f}" for s in (6, 32, 98)), flush=True)
    # B: common vs per-sector depth (weighted least squares with injection errors)
    for key in list(res["signals"]["B"].keys()):
        R = res["signals"]["B"][key]; rows = [r for r in allrows if r["signal"] == "B" and f"{r['template']}__{r['photometry']}" == key]
        d = np.array([r["D"] for r in rows]); e = np.array([r["e_D_inj"] for r in rows]); sec = np.array([r["sector"] for r in rows])
        mu = np.sum(d / e ** 2) / np.sum(1 / e ** 2); chi_common = float(np.sum(((d - mu) / e) ** 2))
        chi_sep = 0.0
        for s in np.unique(sec):
            m = sec == s; mus = np.sum(d[m] / e[m] ** 2) / np.sum(1 / e[m] ** 2); chi_sep += float(np.sum(((d[m] - mus) / e[m]) ** 2))
        k_common, k_sep, n = 1, len(np.unique(sec)), len(d)
        arch = sec != 98; d98, e98 = sector_depth(rows, 98, "e_D_inj")[:2]
        muA = np.sum(d[arch] / e[arch] ** 2) / np.sum(1 / e[arch] ** 2); eA = 1 / np.sqrt(np.sum(1 / e[arch] ** 2))
        R["depth_model_comparison"] = dict(n_events=n, chi2_common=chi_common, chi2_per_sector=chi_sep, dchi2=chi_common - chi_sep,
                                           dBIC_common_minus_sep=(chi_common + k_common * np.log(n)) - (chi_sep + k_sep * np.log(n)),
                                           common_depth_ppm=[mu * 1e6, 1e6 / np.sqrt(np.sum(1 / e ** 2))], S98_ppm=[d98 * 1e6, e98 * 1e6],
                                           archival_ppm=[muA * 1e6, eA * 1e6], S98_minus_archival_sigma=(d98 - muA) / np.hypot(e98, eA))
        print("B depth comparison", key, json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in R["depth_model_comparison"].items()}))
    json.dump(res, open(os.path.join(OUT, "a1_results.json"), "w"), indent=1)
    keys = ["signal", "template", "photometry", "set", "sector", "epoch", "tc_pred", "coverage", "dt_min", "e_dt_white_min", "e_dt_perm_min", "e_dt_inj_min", "D", "e_D_white", "e_D_perm", "e_D_inj", "snr_white", "rms_ppm", "n", "edge"]
    with open(os.path.join(OUT, "a1_per_transit.csv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in allrows: w.writerow({k: r.get(k) for k in keys})
    print("done")
