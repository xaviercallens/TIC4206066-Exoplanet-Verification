"""Executable regression tests for every bug found in the earlier TIC 4206066 work and inside this audit.
Run: audit_env/bin/python tests/test_regressions.py   (PRF-specific tests: tests/test_prf_regressions.py)
Each test prints PASS/FAIL with the quantity checked; the file exits non-zero on any failure.
"""
import inspect
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
AUD = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(AUD, "scripts"))
import lcio  # noqa
import transit_tools as tt  # noqa

EPH = json.load(open(os.path.join(AUD, "out", "ephem_current.json")))
results = []


def check(name, ok, detail):
    results.append(dict(test=name, passed=bool(ok), detail=detail))
    print(("PASS " if ok else "FAIL ") + name + " :: " + detail)


# T1  deep_followup lc_analysis.py fitted signal A without masking signal B (two S32 B transits lie 2.0-2.4 h after
#     A transits inside the +-0.35 d windows). The audit windows must exclude B's in-transit cadences.
lcs = [lcio.detrend(lcio.load_product("SPOC2min_S32_PDC"), eph=EPH)]
W = tt.windows(lcs, "A", EPH, half=0.30)
t = np.concatenate([w["t"] for w in W])
n_b_inside = int(lcio.in_transit(t, "B", 0.5, EPH).sum())
check("T1_A_windows_exclude_B_transits", n_b_inside == 0, f"B in-transit cadences inside A windows: {n_b_inside}")

# T2  detrending must mask BOTH signals (otherwise the trend eats part of the other transit)
src = inspect.getsource(lcio.detrend)
check("T2_detrend_masks_both_signals", 'in_transit(t, "A"' in src and 'in_transit(t, "B"' in src, "lcio.detrend masks A and B")

# T3  np.isclose default rtol silently changed the fast null-test scan; fast scan must equal the brute force exactly
import s11_nullB as nb  # noqa
rng = np.random.default_rng(3)
tt_ = np.sort(rng.uniform(1468, 1490, 3000))
ff_ = 1 + rng.normal(0, 5e-4, tt_.size)
ww_ = np.full(tt_.size, 4e6)
a = nb.scan_stat(tt_, ff_, ww_, nb.T98, 1.6 / 24)
b = nb.scan_stat_brute(tt_, ff_, ww_, nb.T98, 1.6 / 24)
fin = np.isfinite(a) & np.isfinite(b)
check("T3_fast_null_scan_equals_brute", np.array_equal(np.isfinite(a), np.isfinite(b)) and np.max(np.abs(a[fin] - b[fin])) < 1e-9,
      f"max |fast-brute| = {np.max(np.abs(a[fin] - b[fin])):.2e}; 'isclose' absent from scan_stat: {'isclose(' not in inspect.getsource(nb.scan_stat)}")

# T4  Mann+2015 [Fe/H] radius relation: published Table 1 printed f = -0.04458 (press error); erratum/arXiv: +0.04458
s04 = open(os.path.join(AUD, "scripts", "s04_stellar.py")).read()
check("T4_mann2015_feh_sign", "(1 + 0.04458 * feh)" in s04, "s04 uses +0.04458")

# T5  TESScut times are barycentred to the CCD centre; product times must agree after the target-specific shift
tim = json.load(open(os.path.join(AUD, "out", "s06_timing.json")))
worst = max(abs(v["residual_s"]["median"]) for s in ("S6", "S32", "S98") for k, v in tim["compare"][s]["products"].items()
            if k.startswith(("TESS-SPOC", "QLP")))
check("T5_times_consistent_after_ccd_centre_shift", worst < 6.0, f"worst median residual SPOC/QLP vs shifted TESScut = {worst:.2f} s")

# T6  QLP barycentric bugs (S74-79, S99-104): S98 residual must be stable across the sector (no drift)
q = tim["compare"]["S98"]["products"]["QLP_S98"]["residual_by_quarter_s"]
check("T6_QLP_S98_no_timing_drift", np.ptp(q) < 1.0, f"QLP S98 residual by quarter {np.round(q, 2).tolist()} s")

# T7  Pecaut & Mamajek table parsing: Msun is column 30 (0-based); M0V must read 0.57, K8V 0.62
pm = {}
for line in open(os.path.join(AUD, "data/refs/EEM_dwarf_UBVIJHK_colors_Teff.txt")):
    p = line.split()
    if p and p[0] in ("M0V", "K8V") and len(p) > 30 and not line.startswith("#"):
        pm.setdefault(p[0], float(p[30]))
check("T7_PM_table_mass_column", pm.get("M0V") == 0.57 and pm.get("K8V") == 0.62, f"parsed {pm}")

# T8  bolometric zero point: m_bol,sun = -26.8167 (Mann+2015) -> M_bol,sun = m + 5 log10(206265) + 5 = 4.7554
st = json.load(open(os.path.join(AUD, "out", "s04_stellar.json")))
check("T8_Mbol_sun", abs(st["Mbol_sun_used"] - 4.7554) < 1e-3, f"M_bol,sun used = {st['Mbol_sun_used']:.4f}")

# T9  Gaussian prior on rho while sampling ln(rho): log-density must be log N(rho) + ln(rho) (Jacobian)
s08 = open(os.path.join(AUD, "scripts", "s08_globalfit.py")).read()
x = np.linspace(np.log(1.5), np.log(3.8), 20001)
dens = np.exp(-0.5 * ((np.exp(x) - 2.507) / 0.23) ** 2 + x)
dens /= np.trapezoid(dens, x)
mean_rho = np.trapezoid(np.exp(x) * dens, x)
check("T9_lnrho_prior_jacobian", "+ np.log(rho)" in s08 and abs(mean_rho - 2.507) < 0.01, f"implied mean rho {mean_rho:.4f} (target 2.507)")

# T10 S6: every A transit lies 0.4-3.5 h from a momentum dump -> S6 excluded from the headline localisation
d6 = json.load(open(os.path.join(AUD, "out", "s12d_s6_dumps.json")))["A"]
loc = json.load(open(os.path.join(AUD, "out", "s12c_clean_localisation.json")))
check("T10_S6_dump_coincidence_flagged", max(abs(h) for h in d6["transit_minus_dump_h"]) < 4 and "combined_without_S6" in loc,
      f"transit-dump offsets (h) {np.round(d6['transit_minus_dump_h'], 2).tolist()}")

# T11 TRICERATOPS exptime must equal the bin width of a binned fold (6 min) in the audit variants
s13 = open(os.path.join(AUD, "scripts", "s13_fpp_run.py")).read()
check("T11_fpp_exptime_equals_bin", "EXPT = 6 / 1440" in s13, "base/loc/molusc/tic variants use 6-min exptime")

# T12 S98 localisation headline uses clean intervals only (QLP bits 29/30 + straylight clear) for data AND injections
check("T12_S98_clean_localisation", loc["clean_S98"]["A"]["clean"] is True and loc["clean_S98"]["B"]["clean"] is True, "clean flag set")

# T13 localisation significances must come from injection-recovery in real pixels, not from bootstrap over epochs
s12c = open(os.path.join(AUD, "scripts", "s12c_clean_localisation.py")).read()
check("T13_injection_based_errors", "inject_on_target" in s12c and "bootstrap" not in s12c.split("def localise")[1][:3000],
      "localise() uses real-pixel injections")

# T14 (audit bug C7) epoch-only null trials at the fixed S98 period overlap: the +-0.004 d scan spans +-0.9 d in S6, so
#     the trials contain only a handful of independent windows and their upper percentiles saturate at one feature.
#     The null that is quoted must randomise the period as well (s11b) and yield many distinct windows.
s6 = nb.prep(lcio.load_product("TESSSPOC_S6_PDC"))
rng14 = np.random.default_rng(14)


def argmax_cfg(t, f, w, T0, P0, grid):
    old = nb.DPS
    nb.DPS = grid
    try:
        sc = nb.scan_stat(t, f, w, T0, 1.6 / 24, P0)
    finally:
        nb.DPS = old
    if not np.isfinite(sc).any():
        return None, None
    i = int(np.nanargmax(np.where(np.isfinite(sc), sc, -np.inf)))
    Pm = P0 + grid[i]
    inn = np.abs((((t - T0) / Pm + 0.5) % 1 - 0.5) * Pm) < 0.8 / 24
    return tuple(int(round((T0 + c * Pm) / (0.8 / 24))) for c in np.unique(np.round((t - T0) / Pm)[inn])), float(sc[i])


ep_cfg, ep_val = set(), []
for u in rng14.uniform(0.1, 0.9, 1500):
    c, v = argmax_cfg(*s6, nb.T98 + u * nb.P98, nb.P98, nb.DPS)
    if c is not None:
        ep_cfg.add(c), ep_val.append(v)
rp_cfg = set()
for _ in range(1500):
    P0 = rng14.uniform(9.0, 13.5)
    c, _v = argmax_cfg(*s6, nb.T98 + rng14.uniform() * P0, P0, nb.DPS * P0 / nb.P98)
    if c is not None:
        rp_cfg.add(c)
nb11b = json.load(open(os.path.join(AUD, "out", "s11b_nullB_periods.json")))
sat = abs(np.percentile(ep_val, 99) - np.max(ep_val)) < 1e-6 * np.max(ep_val)      # upper tail saturated at one feature
check("T14_null_trials_independent", len(ep_cfg) < 40 and sat and len(rp_cfg) > 300
      and all(r["n_distinct_argmax_configurations"] > 1000 for r in nb11b["results"]),
      f"S6 epoch-only nulls: {len(ep_cfg)} distinct windows in 1500 trials (p99 = max: {sat}); period+epoch nulls: {len(rp_cfg)}; "
      f"s11b distinct windows min {min(r['n_distinct_argmax_configurations'] for r in nb11b['results'])}")

# T15 (audit note C8) s10 draws its random boxes at phases 0.1-0.9 while its phase scan covers 0.05-0.95; the empirical
#     sigma re-derived over 0.05-0.95 (s10c) must agree within 5 % so that no quoted significance depends on it
eb = json.load(open(os.path.join(AUD, "out", "s10_eb_tests.json")))
pf = json.load(open(os.path.join(AUD, "out", "s10c_phase_features.json")))
r1 = pf["B_own_phase0.925"]["emp_sigma_ppm"] / eb["B_own"]["primary"]["err_emp_robust_ppm"]
r2 = pf["A_own_phase0.27"]["emp_sigma_ppm"] / eb["A_own"]["primary"]["err_emp_robust_ppm"]
check("T15_random_box_phase_range", abs(r1 - 1) < 0.05 and abs(r2 - 1) < 0.05, f"sigma ratio (0.05-0.95 vs 0.1-0.9): B own {r1:.3f}, A own {r2:.3f}")

json.dump(results, open(os.path.join(HERE, "test_regressions_result.json"), "w"), indent=1)
bad = [r for r in results if not r["passed"]]
print(f"{len(results) - len(bad)}/{len(results)} passed")
sys.exit(1 if bad else 0)
