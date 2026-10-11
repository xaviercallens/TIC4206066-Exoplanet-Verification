"""T1_A_timing rerun: weighted linear ephemerides of signal A from STORED per-transit times.

Read-only inputs (no source file is written; outputs go to this folder only):
  pubprep   addenda/rnaas_tic4206066_pubprep/analysis/out/a1_per_transit.csv  (+ a1_results.json injection SDs)
  auditS09  addenda/tic4206066_PUBLICATION_AUDIT/out/s09_events.json          (red-noise-inflated Delta-chi2=1 errors)
  U02aud    addenda/opus55_final_audit_tic4206066/02_timing_A/out/t3_events_base.csv (+ t6/t5 injection SDs)
  U02ver    addenda/opus55_final_audit_tic4206066/02_timing_A/verify/v2_real_events.csv (+ v4_results.json SDs)
  URS6      addenda/opus55_final_audit_tic4206066/referee/resolve_R-S6/out/s1_timing_R60.json / _R144.json
For each pipeline and event selection: WLS P, T0 at epoch 0 (BTJD ~3998.737), errors, chi2/dof, rms O-C,
T0 at the covariance-minimising epoch, and the effect of dropping Sector 6.
"""
import json, os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/Users/paulsmith/space <3/tess-jev"
OUT = os.path.dirname(os.path.abspath(__file__))
PUB = f"{ROOT}/addenda/rnaas_tic4206066_pubprep/analysis/out"
AUD = f"{ROOT}/addenda/tic4206066_PUBLICATION_AUDIT/out"
U = f"{ROOT}/addenda/opus55_final_audit_tic4206066"


def wls(ep, t, e):
    ep = np.asarray(ep, float); t = np.asarray(t, float); e = np.asarray(e, float)
    A = np.vstack([np.ones_like(ep), ep]).T
    W = 1 / e ** 2
    C = np.linalg.inv(A.T @ (A * W[:, None]))
    c = C @ (A.T @ (W * t))
    r = t - A @ c
    chi2 = float(np.sum((r / e) ** 2))
    dof = len(t) - 2
    # covariance-minimising (decorrelated) epoch
    n_dec = -C[0, 1] / C[1, 1]
    n_int = int(np.round(n_dec))
    T_dec = c[0] + c[1] * n_int
    eT_dec = float(np.sqrt(C[0, 0] + 2 * n_int * C[0, 1] + n_int ** 2 * C[1, 1]))
    return dict(n=int(len(t)), P=float(c[1]), eP=float(np.sqrt(C[1, 1])), T0_BTJD=float(c[0]),
                T0_BJD=float(c[0] + 2457000), eT0=float(np.sqrt(C[0, 0])),
                corr_T0_P=float(C[0, 1] / np.sqrt(C[0, 0] * C[1, 1])),
                chi2=chi2, dof=dof, p_upper=float(stats.chi2.sf(chi2, dof)),
                rms_oc_min=float(np.sqrt(np.mean((r * 1440) ** 2))),
                max_abs_z=float(np.max(np.abs(r / e))),
                decorr_epoch=n_int, T_decorr_BJD=float(T_dec + 2457000), eT_decorr=eT_dec,
                oc_min={f"{int(a)}": round(float(b), 2) for a, b in zip(ep, r * 1440)})


res = {}

# ---------------- pubprep a1 (pipeline photometry) ----------------
pt = pd.read_csv(f"{PUB}/a1_per_transit.csv")
pa = pt[(pt.signal == "A") & (pt.photometry == "pipe")].copy()
inj = json.load(open(f"{PUB}/a1_results.json"))["signals"]["A"]["A__pipe"]["injections"]
sd_pub = {6: inj["S6"]["sd_dt_min"], 32: inj["S32"]["sd_dt_min"], 98: inj["S98"]["sd_dt_min"]}
pa["t"] = pa.tc_pred + pa.dt_min / 1440
sel = (pa.coverage >= 0.7) & (pa.snr_white >= 2) & (~pa.edge)
p20 = pa[sel]
e = p20.sector.map(sd_pub) / 1440
res["pubprep_20_inj"] = wls(p20.epoch, p20.t, e)
res["pubprep_20_inj"]["events"] = [f"S{s}:{n}" for s, n in zip(p20.sector, p20.epoch)]
res["pubprep_20_white"] = wls(p20.epoch, p20.t, np.maximum(p20.e_dt_white_min, 0.25) / 1440)
two = p20[~(((p20.sector == 98) & (p20.epoch == -1)) | ((p20.sector == 32) & (p20.epoch == -573)))]
res["pubprep_18_twosided_inj"] = wls(two.epoch, two.t, two.sector.map(sd_pub) / 1440)
res["pubprep_18_twosided_inj"]["note"] = "pubprep pipe has no S6 E-795 row, so the two-sided rule gives 5/5/8 = 18"
noS6 = p20[p20.sector != 6]
res["pubprep_20_noS6_inj"] = wls(noS6.epoch, noS6.t, noS6.sector.map(sd_pub) / 1440)

# ---------------- publication audit s09 (red-noise-inflated Delta-chi2=1 errors) ----------------
s9 = json.load(open(f"{AUD}/s09_events.json"))["A"]["pipe"]["rows"]
r9 = pd.DataFrame([r for r in s9 if r.get("dt_min") is not None and r["coverage"] >= 0.7 and r["snr"] >= 2
                   and not r["timing_at_scan_edge"]])
r9["t"] = r9.tc_pred + r9.dt_min / 1440
res["audit_s09_21_rednoise"] = wls(r9.epoch, r9.t, np.maximum(r9.e_dt_min, 0.5) / 1440)
res["audit_s09_21_rednoise"]["events"] = [f"S{s}:{n}" for s, n in zip(r9.sector, r9.epoch)]
res["audit_s09_21_rednoise"]["median_e_dt_min"] = {int(s): float(g.e_dt_min.median()) for s, g in r9.groupby("sector")}
# same times, injection-calibrated per-sector SDs (pubprep and U verifier values)
res["audit_s09_21_pubprepINJ"] = wls(r9.epoch, r9.t, r9.sector.map(sd_pub) / 1440)
r9b = r9[~(((r9.sector == 98) & (r9.epoch == -1)) | ((r9.sector == 32) & (r9.epoch == -573)))]
res["audit_s09_19_twosided_pubprepINJ"] = wls(r9b.epoch, r9b.t, r9b.sector.map(sd_pub) / 1440)
res["audit_s09_19_twosided_rednoise"] = wls(r9b.epoch, r9b.t, np.maximum(r9b.e_dt_min, 0.5) / 1440)

# ---------------- U 02 auditor (t3_events_base; dt = R144, dt_alt = R60) ----------------
t3 = pd.read_csv(f"{U}/02_timing_A/out/t3_events_base.csv")
t6 = json.load(open(f"{U}/02_timing_A/out/t6_results.json"))["ttv"]
sd_aud = {R: {int(k): v for k, v in t6[R]["sigmas_min"].items()} for R in ("R60", "R144")}
aud20 = t3[t3.sel]
for R, col in (("R60", "tc_obs_alt"), ("R144", "tc_obs")):
    res[f"U02aud_20_{R}"] = wls(aud20.epoch, aud20[col], aud20.sector.map(sd_aud[R]) / 1440)
    a19 = aud20[~((aud20.sector == 32) & (aud20.epoch == -573))]
    res[f"U02aud_19_twosided_{R}"] = wls(a19.epoch, a19[col], a19.sector.map(sd_aud[R]) / 1440)
    a19n = a19[a19.sector != 6]
    res[f"U02aud_19_twosided_noS6_{R}"] = wls(a19n.epoch, a19n[col], a19n.sector.map(sd_aud[R]) / 1440)
res["U02aud_19_twosided_R60"]["events"] = [f"S{s}:{n}" for s, n in zip(a19.sector, a19.epoch)]

# ---------------- U 02 verifier (v2_real_events; dt60 / dt144) ----------------
v2 = pd.read_csv(f"{U}/02_timing_A/verify/v2_real_events.csv").set_index(["sector", "epoch"])
v4 = json.load(open(f"{U}/02_timing_A/verify/v4_results.json"))["injections"]
set19 = [(6, -795), (6, -794), (6, -793), (6, -791), (6, -790), (6, -789), (32, -572), (32, -571), (32, -570),
         (32, -569), (32, -568), (98, 0), (98, 1), (98, 4), (98, 5), (98, 9), (98, 10), (98, 13), (98, 14)]
pkg20 = [k for k in set19 if k != (6, -795)] + [(32, -573), (98, -1)]
for R in (60, 144):
    sd = {s: v4[f"S{s}_R{R}"]["sd"] for s in (6, 32, 98)}
    for name, S in (("19_twosided", set19), ("pkg20", pkg20), ("19_twosided_noS6", [k for k in set19 if k[0] != 6])):
        d = v2.loc[S].reset_index()
        res[f"U02ver_{name}_R{R}"] = wls(d.epoch, d.tc + d[f"dt{R}"] / 1440, d.sector.map(sd) / 1440)

# ---------------- U referee resolve_R-S6 pipeline ----------------
for R in ("R60", "R144"):
    j = json.load(open(f"{U}/referee/resolve_R-S6/out/s1_timing_{R}.json"))
    ev = pd.DataFrame(j["events"])
    sdr = {int(k): v for k, v in j["sigma_used_min"].items()}
    ev = ev[ev.sel]
    ev["t"] = ev.tc_pred + ev.dt
    res[f"URS6_all_sel_{R}"] = wls(ev.epoch, ev.t, ev.sector.map(sdr) / 1440)
    res[f"URS6_all_sel_{R}"]["events"] = [f"S{s}:{n}" for s, n in zip(ev.sector, ev.epoch)]
    e19 = ev[~(((ev.sector == 98) & (ev.epoch == -1)) | ((ev.sector == 32) & (ev.epoch == -573)))]
    res[f"URS6_19_twosided_{R}"] = wls(e19.epoch, e19.t, e19.sector.map(sdr) / 1440)
    e19n = e19[e19.sector != 6]
    res[f"URS6_19_twosided_noS6_{R}"] = wls(e19n.epoch, e19n.t, e19n.sector.map(sdr) / 1440)
    res[f"URS6_19_twosided_{R}"]["events"] = [f"S{s}:{n}" for s, n in zip(e19.sector, e19.epoch)]

# ---------------- summary table ----------------
rows = []
for k, v in res.items():
    rows.append(dict(key=k, n=v["n"], P=round(v["P"], 8), eP=f"{v['eP']:.2e}", T0_BJD=round(v["T0_BJD"], 5),
                     eT0=round(v["eT0"], 5), chi2=round(v["chi2"], 2), dof=v["dof"], rms=round(v["rms_oc_min"], 2),
                     decorr_n=v["decorr_epoch"], T_dec=round(v["T_decorr_BJD"], 5), eT_dec=round(v["eT_decorr"], 5)))
tab = pd.DataFrame(rows)
print(tab.to_string(index=False))
tab.to_csv(f"{OUT}/r1_ephemerides_table.csv", index=False)

# ---------------- S6-drop effect, spread of 19-set across pipelines ----------------
extra = {}
for base, nos6 in (("U02aud_19_twosided_R60", "U02aud_19_twosided_noS6_R60"),
                   ("U02aud_19_twosided_R144", "U02aud_19_twosided_noS6_R144"),
                   ("U02ver_19_twosided_R60", "U02ver_19_twosided_noS6_R60"),
                   ("U02ver_19_twosided_R144", "U02ver_19_twosided_noS6_R144"),
                   ("URS6_19_twosided_R60", "URS6_19_twosided_noS6_R60"),
                   ("pubprep_20_inj", "pubprep_20_noS6_inj")):
    a, b = res[base], res[nos6]
    extra[f"dropS6:{base}"] = dict(dP=b["P"] - a["P"], dP_over_eP=(b["P"] - a["P"]) / a["eP"],
                                   dT0_min=(b["T0_BTJD"] - a["T0_BTJD"]) * 1440,
                                   dT0_over_eT0=(b["T0_BTJD"] - a["T0_BTJD"]) / a["eT0"],
                                   eP_ratio=b["eP"] / a["eP"], chi2_noS6=f"{b['chi2']:.2f}/{b['dof']}")
k19 = ["U02aud_19_twosided_R60", "U02aud_19_twosided_R144", "U02ver_19_twosided_R60", "U02ver_19_twosided_R144",
       "URS6_19_twosided_R60", "audit_s09_19_twosided_pubprepINJ"]
extra["19set_cross_pipeline"] = dict(P=[res[k]["P"] for k in k19], T0=[res[k]["T0_BJD"] for k in k19],
                                     P_range=[min(res[k]["P"] for k in k19), max(res[k]["P"] for k in k19)],
                                     T0_range=[min(res[k]["T0_BJD"] for k in k19), max(res[k]["T0_BJD"] for k in k19)],
                                     chi2=[f"{res[k]['chi2']:.2f}/{res[k]['dof']}" for k in k19], keys=k19)
# pixel-recheck ephemeris (audit ephem_current.json) vs final 19-set ephemeris: max time offset over the data
Pp, Tp = 3.182785293657263, 3998.7367930333035
f = res["U02ver_19_twosided_R60"]
off = {s: (Tp + Pp * n - (f["T0_BTJD"] + f["P"] * n)) * 1440 for s, n in (("S6", -792), ("S32", -570), ("S98", 7))}
extra["pixel_recheck_ephem_minus_final_min"] = off
# manuscript pubprep ephemeris vs final at S110 (epoch ~ (4-Nov-2026) ) : predicted offset
extra["pubprep_minus_final_at_E0_min"] = (2460998.7375446941755 - 2457000 - f["T0_BTJD"]) * 1440
print(json.dumps(extra, indent=1))
json.dump(dict(results=res, extra=extra), open(f"{OUT}/r1_ephemerides.json", "w"), indent=1)
