"""T5 signal B minimal reruns (read-only inputs; outputs only in this folder).
1. Linear ephemeris fits with full covariance from U b08 per-event times (4-event manuscript set, 7-event set, S98-only).
2. S110 predictions (BJD_TDB and geocentric UTC) with 1-sigma from the covariance (formal and conservative).
3. P_B/P_A ratio; S32 B-vs-A contact timing.
4. p-value -> sigma conversions, Fisher combinations, look-elsewhere for S98-only FAP.
5. File-time chronology of the post hoc sequence."""
import json, os, numpy as np
from scipy import stats
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
ROOT = "/Users/paulsmith/space <3/tess-jev"
U = ROOT + "/addenda/opus55_final_audit_tic4206066"
OUT = ROOT + "/addenda/tic4206066_FINAL_RECONCILIATION/_work/T5_signal_B"
res = {}
b08 = json.load(open(U + "/04_B_existence/out/b08_ephem.json"))
ev = {(e["sector"], e["epoch"]): e for e in b08["events"] if e.get("used")}
def linfit(keys):
    n = np.array([k[1] for k in keys], float); t = np.array([ev[k]["tc"] for k in keys]); s = np.array([ev[k]["err_min"] for k in keys]) / 1440
    A = np.vstack([np.ones_like(n), n]).T; W = np.diag(1 / s**2)
    C = np.linalg.inv(A.T @ W @ A); p = C @ A.T @ W @ t
    chi2 = float(np.sum(((t - A @ p) / s) ** 2)); dof = len(n) - 2
    scale = max(1.0, np.sqrt(chi2 / dof)) if dof > 0 else 1.0
    return dict(T0_btjd=p[0], P=p[1], eT0=np.sqrt(C[0, 0]), eP=np.sqrt(C[1, 1]), cov=C.tolist(), corr=C[0, 1] / np.sqrt(C[0, 0] * C[1, 1]),
                chi2=chi2, dof=dof, scale=scale, n_events=len(n), events=[list(k) for k in keys])
fits = {"manuscript_4ev": linfit([("6", -227), ("6", -226), ("98", 0), ("98", 1)]),
        "seven_ev": linfit([("6", -227), ("6", -226), ("32", -164), ("32", -162), ("98", 0), ("98", 1), ("98", 4)]),
        "S98_only_3ev": linfit([("98", 0), ("98", 1), ("98", 4)]),
        "S98_only_2ev_E0E1": None}
k2 = [("98", 0), ("98", 1)]
t0, t1 = ev[k2[0]]["tc"], ev[k2[1]]["tc"]; s0, s1 = ev[k2[0]]["err_min"] / 1440, ev[k2[1]]["err_min"] / 1440
fits["S98_only_2ev_E0E1"] = dict(P=t1 - t0, eP=float(np.hypot(s0, s1)), T0_btjd=t0, eT0=s0)
res["fits"] = fits
# final (adopted) ephemeris: U/pubprep conservative values
PB, ePB, T0B, eT0B = 11.13274, 0.00003, 4000.865, 0.003
res["adopted"] = dict(P=PB, eP=ePB, T0_btjd=T0B, eT0=eT0B, T0_bjd=T0B + 2457000)
# S110 predictions
ra, dec = 80.85983, -1.32696   # from U/01_ephemeris_A/s13_s110.py L11
target = SkyCoord(ra * u.deg, dec * u.deg); geo = EarthLocation.from_geocentric(0 * u.m, 0 * u.m, 0 * u.m)
def to_utc(bjd):
    t = Time(bjd, format="jd", scale="tdb", location=geo)
    for _ in range(3):
        ltt = t.light_travel_time(target, kind="barycentric"); t = Time(bjd, format="jd", scale="tdb", location=geo) - ltt
    return t.utc.iso
win = (Time("2026-10-31 00:00:00", scale="utc").tdb.jd - 2457000, Time("2026-11-27 00:00:00", scale="utc").tdb.jd - 2457000)
rows = []
f4 = fits["manuscript_4ev"]; C4 = np.array(f4["cov"]) 
fac = ePB / f4["eP"]   # scale formal covariance so that sigma_P equals the adopted 3e-5
fS = fits["S98_only_3ev"]; CS = np.array(fS["cov"]) * fS["scale"] ** 2
for n in range(29, 36):
    tb = T0B + n * PB
    if not (win[0] - 1 <= tb <= win[1] + 1): continue
    s_uncorr = np.hypot(eT0B, n * ePB) * 1440
    s_cov = np.sqrt(C4[0, 0] + n**2 * C4[1, 1] + 2 * n * C4[0, 1]) * 1440
    t4 = f4["T0_btjd"] + n * f4["P"]
    tS = fS["T0_btjd"] + n * fS["P"]; sS = np.sqrt(CS[0, 0] + n**2 * CS[1, 1] + 2 * n * CS[0, 1]) * 1440
    rows.append(dict(epoch=n, btjd=tb, bjd_tdb=tb + 2457000, utc_geocentre=to_utc(tb + 2457000), inside_Oct31_Nov26=bool(win[0] <= tb <= win[1] - 1),
                     sigma_min_adopted_uncorrelated=s_uncorr, sigma_min_4ev_formal_cov=s_cov, sigma_min_4ev_cov_scaled_to_eP_3e5=s_cov * fac,
                     t_4ev_fit_minus_adopted_min=(t4 - tb) * 1440, S98only_pred_minus_adopted_min=(tS - tb) * 1440, S98only_sigma_min=sS))
res["S110"] = dict(window_utc="2026-10-31 00:00 to 2026-11-26 24:00 (HEASARC sector110_summary.html, accessed 2026-09-25T15:30Z)", rows=rows)
# ratio with final A ephemeris (U revision, P_A = 3.182785 +- 0.000006)
PA, ePA = 3.182785, 0.000006
r = PB / PA; er = r * np.hypot(ePB / PB, ePA / PA)
res["ratio"] = dict(PB_over_PA=r, err=er, offset_from_7_2_percent=(r / 3.5 - 1) * 100,
                    audit_values=dict(PB=11.132745380039134, PA_ms=3.1827856, r=11.132745380039134 / 3.1827856))
# S32 contact timing (A: P, T0 from revision; T14_A = 1.87 h; B T14 1.62-1.66 h)
T0A = 2460998.736 - 2457000; T14A = 1.87; T14B = 1.66
s32 = []
for nB in (-164, -162):
    tB = T0B + nB * PB; nA = round((tB - T0A) / PA); tA = T0A + nA * PA
    dt = (tB - tA) * 24
    s32.append(dict(epochB=nB, tB_btjd=tB, epochA=nA, tA_btjd=tA, B_minus_A_mid_h=dt,
                    gap_Aend_to_Bstart_h=dt - T14A / 2 - T14B / 2, fitted_tB_b08=ev[("32", nB)]["tc"],
                    fitted_gap_h=(ev[("32", nB)]["tc"] - tA) * 24 - T14A / 2 - T14B / 2))
res["S32_B_vs_A"] = s32
# significance arithmetic
z = lambda p: float(stats.norm.isf(p))
ps = {"link_S6S32_pubprep_19of10000": 1.9998e-3, "link_S6S32_audit_24of20000": 1.25e-3, "link_trials2_pubprep": 4e-3,
      "link_window_pm0.011_manuscript_null_upper": 6e-3, "link_clean_null_upper": 1e-4, "link_clean_null_lower": 1e-5,
      "S6_rank_pubprep": 1.5521e-3, "S6_rank_U": 2.0e-3, "S6_fullprior_null_verify": 1.14e-3, "S32_only": 0.032}
res["sigma_one_sided"] = {k: dict(p=v, z=z(v)) for k, v in ps.items()}
def fisher(*p):
    X = -2 * np.sum(np.log(p)); return float(stats.chi2.sf(X, 2 * len(p)))
res["fisher"] = {"S98_0.35_x_link_2e-3": fisher(0.35, 2e-3), "S98_1.6e-3_x_link_2e-3": fisher(1.6e-3, 2e-3)}
res["fisher"] = {k: dict(p=v, z=z(v)) for k, v in res["fisher"].items()}
N = 126246
res["S98_only_lookelsewhere"] = {str(p): dict(per_star=p, expected_false_alarms_in_search=p * N) for p in (2e-4, 1.6e-3, 0.01, 0.35)}
# chronology
chron = {}
for rel in ["addenda/pool_multisector/multisector_results_pass1.json", "addenda/pool_multisector/multisector_results.json", "addenda/pool_multisector/SUMMARY.md",
            "addenda/deep_followup_tic4206066/out", "addenda/deep_followup_tic4206066/signal_11d.py", "addenda/deep_followup_tic4206066/out/signal_11d.json",
            "addenda/deep_followup_tic4206066/signal_11d_null.py", "addenda/deep_followup_tic4206066/out/signal_11d_null.json",
            "addenda/deep_followup_tic4206066/deep_followup_report.md"]:
    st = os.stat(ROOT + "/" + rel)
    chron[rel] = dict(birth_local=Time(st.st_birthtime, format="unix").to_datetime().isoformat() + "Z", mtime_utc=Time(st.st_mtime, format="unix").to_datetime().isoformat() + "Z")
res["chronology_utc"] = chron
json.dump(res, open(OUT + "/t5_rerun.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
