"""Audit 9, step 3: Gaia DR3 RV constraints -> K and M sin i limits vs period.
Field definitions: Gaia DR3 documentation, gaia_source datamodel (downloaded 2026-09-25, sha256 in REPORT):
  radial_velocity_error = sqrt(sigma_med^2 + 0.11^2), sigma_med = sqrt(pi/2) sigma(V_j)/sqrt(rv_nb_transits)
  (rv_method_used = 1; Katz et al. 2023, A&A 674, A5);  rv_amplitude_robust = max - min of epoch RVs after
  3xIQR outlier removal; rv_chisq_pvalue = P(chi2 > observed | constant RV), k = N-1;
  rv_renormalised_gof = 4.5 sqrt(nu) (RUWE^(2/3) + 2/(9 nu) - 1), RUWE = UWE/mode(UWE) in the (grvs_mag, template Teff) bin;
  suggested variability selection: rv_chisq_pvalue < 0.01 & rv_renormalised_gof > 4 & rv_nb_transits >= 10.
Limit method (own): Monte-Carlo of Gaia-like epochs (33 transits in 15 visibility periods over rv_time_duration =
922.9 d; pairs of field-of-view transits 106.5 min apart, 6 h spin) with white epoch noise equal to the typical
epoch scatter of single stars in this bin; circular orbit, random phase. K_95(P) = K for which 95% of simulations give
a sample standard deviation larger than the observed sigma(V_j) (and, separately, a larger robust amplitude than
rv_amplitude_robust). M sin i from the mass function with M1 = 0.606 Msun."""
import json
from pathlib import Path
import numpy as np
from astropy.table import Table
from scipy.stats import chi2 as chi2d
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
g = Table.read(HERE / "data/gaia_dr3_source.ecsv")[0]
eps, N, amp = float(g["radial_velocity_error"]), int(g["rv_nb_transits"]), float(g["rv_amplitude_robust"])
pval, rgof, dur, nvis = float(g["rv_chisq_pvalue"]), float(g["rv_renormalised_gof"]), float(g["rv_time_duration"]), int(g["rv_visibility_periods_used"])
nu = N - 1
sigV = np.sqrt((eps ** 2 - 0.11 ** 2) * 2 * N / np.pi)
chi2_obs = chi2d.isf(pval, nu)                 # chi2 statistic implied by the p-value
uwe_p = np.sqrt(chi2_obs / nu)                 # unit-weight error implied by the p-value
ruwe_rv = (rgof / (4.5 * np.sqrt(nu)) + 1 - 2 / (9 * nu)) ** 1.5
sig_eps = sigV / uwe_p                         # implied mean formal epoch uncertainty
sig_typ = sigV / ruwe_rv                       # typical epoch scatter of stars in this (grvs, Teff) bin
out = dict(radial_velocity_error=eps, rv_nb_transits=N, rv_amplitude_robust=amp, rv_chisq_pvalue=pval,
           rv_renormalised_gof=rgof, rv_time_duration_d=dur, rv_visibility_periods_used=nvis,
           sigma_Vj_kms=sigV, chi2_from_pvalue=chi2_obs, UWE_from_pvalue=uwe_p, RUWE_rv_from_rgof=ruwe_rv,
           implied_formal_epoch_sigma_kms=sig_eps, typical_epoch_scatter_kms=sig_typ,
           excess_scatter_kms=float(np.sqrt(max(sigV ** 2 - sig_typ ** 2, 0))),
           katz_variability_criterion=bool(pval < 0.01 and rgof > 4 and N >= 10))
out["K_if_excess_is_orbital_kms"] = out["excess_scatter_kms"] * np.sqrt(2)
out["K_conservative_all_scatter_kms"] = float(sigV * np.sqrt(2))

rng = np.random.default_rng(42)
def epochs(nsim):
    ts = []
    for _ in range(nsim):
        c = np.sort(np.concatenate([[0, dur], rng.uniform(0, dur, nvis - 2)]))
        per = np.full(nvis, N // nvis); per[rng.choice(nvis, N - per.sum(), replace=False)] += 1
        t = []
        for ci, npv in zip(c, per):
            base = ci + rng.uniform(0, 3)
            for j in range(npv):                       # alternate FoV (106.5 min) and spin (6 h) steps
                t.append(base + (j // 2) * 0.25 + (j % 2) * 106.5 / 1440)
        ts.append(np.array(t)[:N])
    return np.array(ts)

NS = 4000
T = epochs(NS)
noise = rng.normal(0, sig_typ, (NS, N))
def stats_for(K, P):
    ph = rng.uniform(0, 2 * np.pi, NS)[:, None]
    v = K * np.sin(2 * np.pi * T / P + ph) + noise
    sd = v.std(axis=1, ddof=1)
    q1, q3 = np.percentile(v, [25, 75], axis=1); iqr = q3 - q1
    vm = np.where((v < (q1 - 3 * iqr)[:, None]) | (v > (q3 + 3 * iqr)[:, None]), np.nan, v)
    am = np.nanmax(vm, axis=1) - np.nanmin(vm, axis=1)
    return sd, am

# calibration check: K = 0 distribution vs observed
sd0, am0 = stats_for(0.0, 10.0)
out["K0_sim"] = dict(sd_median=float(np.median(sd0)), frac_sd_ge_obs=float(np.mean(sd0 >= sigV)),
                     amp_median=float(np.median(am0)), frac_amp_ge_obs=float(np.mean(am0 >= amp)))
def klim(P, which, conf=0.95, obs=None):
    def f(K):
        sd, am = stats_for(K, P)
        s = sd if which == "sd" else am
        return np.mean(s > obs) - conf
    lo, hi = 0.0, 5.0
    while f(hi) < 0 and hi < 2000: hi *= 2
    return brentq(f, lo, hi, xtol=0.01) if f(lo) < 0 else 0.0

G = 6.67430e-11; MSUN = 1.98841e30
def msini(K_kms, P_d, M1=0.606):
    # solve K = (2 pi G / P)^(1/3) m sin i / (M1 + m)^(2/3) (circular, sin i = 1 for m)
    fn = lambda m: (2 * np.pi * G / (P_d * 86400)) ** (1 / 3) * m * MSUN / ((M1 + m) * MSUN) ** (2 / 3) - K_kms * 1e3
    return brentq(fn, 1e-6, 50)
def K_of(m, P_d, M1=0.606, sini=1.0):
    return (2 * np.pi * G / (P_d * 86400)) ** (1 / 3) * m * sini * MSUN / ((M1 + m) * MSUN) ** (2 / 3) / 1e3

rows = []
for P in [1.0, 3.1827856, 6.365571, 11.13274, 22.265, 30, 100, 300, 1000, 3000, 10000]:
    Ksd = klim(P, "sd", obs=sigV); Kam = klim(P, "amp", obs=amp)
    rows.append(dict(P_d=P, K95_sd_kms=Ksd, K95_amp_kms=Kam, Msini95_sd_Msun=msini(Ksd, P), Msini95_amp_Msun=msini(Kam, P),
                     Msini_if_K3p2_Msun=msini(3.2, P), K_for_0p08Msun_kms=K_of(0.08, P), K_for_0p38Msun_kms=K_of(0.38, P)))
    print({k: round(v, 4) for k, v in rows[-1].items()}, flush=True)
out["limits_vs_period"] = rows
out["manuscript_check"] = dict(Msini_K3p2_P3p18=msini(3.2, 3.1827856), Msini_K3p2_P11p13=msini(3.2, 11.13274),
                               K_0p08_P3p18=K_of(0.08, 3.1827856), K_0p08_P6p37=K_of(0.08, 6.365571))
print(json.dumps({k: v for k, v in out.items() if k != "limits_vs_period"}, indent=1, default=float))
json.dump(out, open(HERE / "out/rv_limits.json", "w"), indent=1, default=float)
