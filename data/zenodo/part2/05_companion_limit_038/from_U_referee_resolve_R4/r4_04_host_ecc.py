"""R4 (d),(e): circular-orbit host-mass threshold for the 3.18 d signal, and the eccentricity needed by lighter hosts.

(d) Undiluted threshold: Mann+2015 R(M_K) x Mann+2019 M(M_K) (coefficients as verified by audit 06 from the arXiv
    sources) solved for rho_host = rho_max, with rho_max = 3.5 (manuscript, S6+S32+S98), 3.60 (verifier 09, diluted-host
    fits, dchi2 = 4) and 3.64 (referee pass1, S32+S98 only, 97.7%). Diluted fits (audit 09 and its verifier) are quoted
    from their JSON files.
(e) Photo-eccentric requirement (Kipping 2010; Dawson & Johnson 2012): the transit speed relative to circular is
    g = (1 + e sin w)/sqrt(1 - e^2); T ~ 1/g and a circular fit returns rho_circ = rho_host * g^3. A host denser than
    rho_max therefore needs g <= (rho_max / rho_host)^(1/3) < 1 (transit near apoapse, w ~ 270 deg).
    e_min (w = 270 deg): g = sqrt((1 - e)/(1 + e)) -> e_min = (1 - g^2)/(1 + g^2).
    Prior fraction: e ~ Beta(0.867, 3.03) (TRICERATOPS 1.1.0 sample_ecc, planet=True; Kipping 2013), w ~ U(0, 360 deg),
    unweighted and weighted by the eccentric transit-probability factor (1 + e sin w)/(1 - e^2)... for a transit at true
    anomaly f = 90 - w the geometric factor is (1 + e cos f)/(1 - e^2) = (1 + e sin w)/(1 - e^2).
Tidal e-damping time (Goldreich & Soter 1966 form, tau_e = (4/63) Q'_p (m_p/M_*) (a/R_p)^5 / n) for context.
audit_env python.
"""
import json, os
import numpy as np
from scipy.optimize import brentq

H = os.path.dirname(os.path.abspath(__file__))
AUD = os.path.abspath(os.path.join(H, "..", ".."))
res = {}
R15 = lambda mk: 1.9515 - 0.3520 * mk + 0.01680 * mk**2
a6 = np.array([-0.642, -0.208, -8.43e-4, 7.87e-3, 1.42e-4, -2.13e-4])
M19 = lambda mk: 10 ** sum(a6[i] * (mk - 7.5) ** i for i in range(6))
rho = lambda mk: M19(mk) / R15(mk) ** 3
thr = {}
for rmax in (3.5, 3.6, 3.64):
    mk = brentq(lambda x: rho(x) - rmax, 4.6, 7.5)
    thr[f"rho_max_{rmax}"] = {"MK": mk, "M": M19(mk), "R": R15(mk)}
res["undiluted_threshold_Mann"] = thr

v = json.load(open(os.path.join(AUD, "09_companion/verify/v3_fit_A.json")))
res["verifier09_diluted_fits_A"] = {"rho_max_dchi2_4": v["target_rho_max_dchi2_4"],
                                    "hosts": [{k: h[k] for k in ("M", "R", "rho", "dchi2_circ", "dchi2_free", "rho_app")} for h in v["hosts"]]}
vb = json.load(open(os.path.join(AUD, "09_companion/verify/v3_fit_B.json")))
res["verifier09_diluted_fits_B"] = {"hosts": [{k: h.get(k) for k in ("M", "dchi2_circ", "dchi2_free")} for h in vb["hosts"]]}
# dchi2 = 4 crossing of the circular fit for A (linear interpolation over the grid)
hs = sorted(v["hosts"], key=lambda h: h["M"])
Ms = np.array([h["M"] for h in hs]); dc = np.array([h["dchi2_circ"] for h in hs])
cross = None
for i in range(len(Ms) - 1):
    if dc[i] >= 4 > dc[i + 1]:
        cross = float(Ms[i] + (4 - dc[i]) * (Ms[i + 1] - Ms[i]) / (dc[i + 1] - dc[i]))
res["verifier09_circular_dchi2_4_crossing_A"] = cross

# eccentricity requirement and prior fractions
rng = np.random.default_rng(42)
n = 4_000_000
e = rng.beta(0.867, 3.03, n)
w = rng.uniform(0, 2 * np.pi, n)
g = (1 + e * np.sin(w)) / np.sqrt(1 - e**2)
ptr = (1 + e * np.sin(w)) / (1 - e**2)
ecc = {}
for Mh in (0.14, 0.17, 0.20, 0.25, 0.30, 0.35, 0.38):
    mk = brentq(lambda x: M19(x) - Mh, 4.6, 11.0)
    rh = rho(mk)
    row = {"MK": mk, "rho_host": rh}
    for rmax in (3.5, 3.64):
        gmax = (rmax / rh) ** (1 / 3)
        emin = (1 - gmax**2) / (1 + gmax**2)
        ok = g <= gmax
        row[f"rho_max_{rmax}"] = {"g_max": gmax, "e_min_w270": emin, "prior_frac_e_gt_emin": float((e > emin).mean()),
                                  "prior_frac_g_le_gmax": float(ok.mean()),
                                  "prior_frac_g_le_gmax_transit_weighted": float((ptr * ok).sum() / ptr.sum())}
    ecc[f"{Mh:.2f}"] = row
res["eccentricity_requirement"] = ecc
res["note_mann_range"] = "Mann+2019 valid for M_K 4.5-10.5 (0.075-0.70 Msun); 0.14 Msun is inside."

# tidal damping context: host 0.3 Msun, P = 3.18 d, planet 2.0-2.7 R_earth, 4-8 M_earth
G, Msun, Me, Re, AU, yr = 6.674e-11, 1.989e30, 5.972e24, 6.371e6, 1.496e11, 3.156e7
P = 3.1827856 * 86400
tid = {}
for Mh in (0.2, 0.3):
    a = (G * Mh * Msun * P**2 / (4 * np.pi**2)) ** (1 / 3)
    nn = 2 * np.pi / P
    for Rp, mp in ((2.0, 5.0), (2.7, 8.0)):
        for Qp in (1e2, 1e3, 1e4, 1e5):
            tau = (4 / 63) * Qp * ((mp * Me) / (Mh * Msun)) * (a / (Rp * Re)) ** 5 / nn
            tid[f"M{Mh}_Rp{Rp}_mp{mp}_Q'{Qp:.0e}"] = tau / yr
res["tidal_tau_e_yr"] = tid
json.dump(res, open(os.path.join(H, "out", "r4_04_host_ecc.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if not k.startswith("verifier09_diluted")}, indent=1))
