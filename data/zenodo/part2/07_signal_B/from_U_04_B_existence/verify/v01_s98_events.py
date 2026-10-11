"""V01: Sector 98 predicted B events (manuscript ephemeris): coverage, QLP flag fractions, per-event box depth with an
empirical (random-epoch) error; DET_FLUX vs SYS_RM_FLUX comparison after identical detrending."""
import json, numpy as np
import vcommon as v

rng = np.random.default_rng(11)
d, t, q = v.load("qlp98")
det = np.asarray(d["DET_FLUX"], float); sys_ = np.asarray(d["SYS_RM_FLUX"], float)
bits = {b: ((q >> b) & 1).astype(bool) for b in (5, 11, 29, 30)}
A = v.inA(t)
ok0 = (q == 0) & np.isfinite(det) & np.isfinite(sys_)
okr = np.isfinite(det) & ~bits[11] & ~bits[5] & np.isfinite(sys_)       # relaxed: keep 2^29 and 2^30
det0 = np.where(ok0, det / np.nanmedian(det[ok0]), np.nan)
n = np.arange(-2, 6); tcs = v.T0_MS + n * v.P_MS
tcs_in = [(k, tc) for k, tc in zip(n, tcs) if t.min() - 0.1 < tc < t.max() + 0.1]
mB = np.zeros_like(t, bool)
for _, tc in tcs_in:
    mB |= np.abs(t - tc) < v.DUR
# own detrending of SYS_RM (A and B masked) on clean and relaxed sets
sysd0 = np.full_like(t, np.nan); sysd0[ok0] = v.biweight(t[ok0], sys_[ok0] / np.nanmedian(sys_[ok0]), (A | mB)[ok0])
sysdr = np.full_like(t, np.nan); sysdr[okr] = v.biweight(t[okr], sys_[okr] / np.nanmedian(sys_[okr]), (A | mB)[okr])
detdr = np.full_like(t, np.nan); detdr[okr] = det[okr] / np.nanmedian(det[okr])


def box(tt, ff, tc, half_in=v.DUR / 2, win=3 * v.DUR):
    """depth = linear baseline (fit to |dt| in (half_in+0.01, win)) minus in-box mean (|dt| < half_in)."""
    m = np.isfinite(ff) & (np.abs(tt - tc) < win) & ~v.inA(tt)
    x, y = tt[m] - tc, ff[m]
    inn = np.abs(x) < half_in; out = np.abs(x) > half_in + 0.01
    if inn.sum() < 3 or out.sum() < 10:
        return None
    c = np.polyfit(x[out], y[out], 1)
    base = np.polyval(c, x[inn])
    return float(np.mean(base - y[inn])), int(inn.sum())


def emp_err(tt, ff, n_in_target, ntr=1500):
    """std of box depths at random epochs (no B, no A overlap) with n_in within 15% of the event's n_in."""
    vals = []
    k = 0
    while len(vals) < ntr and k < 40 * ntr:
        k += 1
        tc = rng.uniform(t.min(), t.max())
        if np.min(np.abs(tc - tcs)) < 3 * v.DUR:
            continue
        r = box(tt, ff, tc)
        if r is None or abs(r[1] - n_in_target) > 0.15 * n_in_target:
            continue
        vals.append(r[0])
    vals = np.array(vals)
    return float(v.rstd(vals)), float(np.std(vals)), len(vals)


rows = []
for k, tc in tcs_in:
    w = np.abs(t - tc) < v.DUR / 2
    row = dict(epoch=int(k), tc=float(tc), n_cad=int(w.sum()), n_q0=int((w & (q == 0)).sum()),
               frac_b29=float(bits[29][w].mean()), frac_b30=float(bits[30][w].mean()), frac_b2048=float(bits[11][w].mean()),
               frac_b32=float(bits[5][w].mean()))
    for name, ff in [("DET_q0", det0), ("SYSRM_bw_q0", sysd0), ("SYSRM_bw_relaxed", sysdr), ("DET_relaxed", detdr)]:
        r = box(t, ff, tc)
        if r is None:
            row[name] = None
            continue
        e_r, e_s, ne = emp_err(t, ff, r[1])
        row[name] = dict(D_ppm=r[0] * 1e6, n_in=r[1], err_emp_ppm=e_r * 1e6, err_emp_std_ppm=e_s * 1e6, n_err=ne, snr=r[0] / e_r)
    rows.append(row)
    print(json.dumps(row), flush=True)

# DET vs SYS_RM (own biweight) comparison on clean data, out of A and B
m = ok0 & ~A & ~mB & np.isfinite(sysd0)
diff = det0[m] - sysd0[m]
cmp_ = dict(rstd_det_ppm=v.rstd(det0[m]) * 1e6, rstd_sysbw_ppm=v.rstd(sysd0[m]) * 1e6, rstd_diff_ppm=v.rstd(diff) * 1e6,
            corr=float(np.corrcoef(det0[m], sysd0[m])[0, 1]),
            DETMET=str(d.columns.names))
# where does the difference concentrate? hourly binned residual difference: largest excursions
tb = t[m]; bins = np.arange(tb.min(), tb.max() + 1 / 24, 1 / 24)
idx = np.digitize(tb, bins)
bd = np.array([np.mean(diff[idx == i]) if np.sum(idx == i) > 5 else np.nan for i in range(1, len(bins))])
cmp_["hourly_diff_rstd_ppm"] = v.rstd(bd) * 1e6
cmp_["hourly_diff_max_abs_ppm"] = float(np.nanmax(np.abs(bd)) * 1e6)
top = np.argsort(-np.nan_to_num(np.abs(bd)))[:8]
cmp_["hourly_diff_top"] = [(float(bins[i]), float(bd[i] * 1e6)) for i in top]
print(json.dumps(cmp_, default=str)[:1500])
json.dump(dict(rows=rows, det_vs_sysrm=cmp_, manifest=v.check_manifest()), open(v.HERE + "/out/v01_s98_events.json", "w"), indent=1, default=str)
