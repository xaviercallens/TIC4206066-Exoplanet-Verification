"""Verifier T1_A_timing: independent re-derivation (own WLS code, np.polyfit cross-check) of the
signal-A linear ephemerides from STORED per-transit times. Read-only inputs; writes only next to itself."""
import json, os, numpy as np, pandas as pd
from scipy import stats
ROOT = "/Users/paulsmith/space <3/tess-jev"
HERE = os.path.dirname(os.path.abspath(__file__))
U = f"{ROOT}/addenda/opus55_final_audit_tic4206066"
PUB = f"{ROOT}/addenda/rnaas_tic4206066_pubprep/analysis/out"
AUD = f"{ROOT}/addenda/tic4206066_PUBLICATION_AUDIT/out"

def fit(ep, t, e):
    ep = np.asarray(ep, float); t = np.asarray(t, float); e = np.asarray(e, float)
    # np.polyfit with weights 1/sigma, cov='unscaled' gives the formal WLS covariance
    c, C = np.polyfit(ep, t, 1, w=1/e, cov="unscaled")
    P, T0 = c; eP, eT0 = np.sqrt(np.diag(C)); r = t - (T0 + P*ep)
    chi2 = float(np.sum((r/e)**2)); dof = len(t)-2
    return dict(n=len(t), P=P, eP=eP, T0_BJD=T0+2457000, eT0=eT0, chi2=chi2, dof=dof,
                p_sf=float(stats.chi2.sf(chi2, dof)), rms_min=float(np.sqrt(np.mean((r*1440)**2))),
                corr=float(C[0,1]/np.sqrt(C[0,0]*C[1,1])))

out = {}
v2 = pd.read_csv(f"{U}/02_timing_A/verify/v2_real_events.csv").set_index(["sector","epoch"])
inj = json.load(open(f"{U}/02_timing_A/verify/v4_results.json"))["injections"]
s19 = [(6,-795),(6,-794),(6,-793),(6,-791),(6,-790),(6,-789),(32,-572),(32,-571),(32,-570),(32,-569),(32,-568),
       (98,0),(98,1),(98,4),(98,5),(98,9),(98,10),(98,13),(98,14)]
for R in (60,144):
    sd = {s: inj[f"S{s}_R{R}"]["sd"] for s in (6,32,98)}
    sd1 = {s: inj[f"S{s}_R{R}"]["sd_onesided"] for s in (6,32,98)}
    def run(S, onesided=()):
        d = v2.loc[S].reset_index()
        e = np.array([sd1[s] if (s,n) in onesided else sd[s] for s,n in zip(d.sector,d.epoch)])/1440
        return fit(d.epoch, d.tc + d[f"dt{R}"]/1440, e)
    out[f"ver19_R{R}"] = run(s19)
    # principled alternative: keep S98 E-1 but give it the one-sided-injection SD (U v4 sd_onesided)
    out[f"ver19+Em1_onesidedSD_R{R}"] = run(s19+[(98,-1)], onesided=[(98,-1)])
    out[f"ver19+Em1+Em573_onesidedSD_R{R}"] = run(s19+[(98,-1),(32,-573)], onesided=[(98,-1),(32,-573)])
    out[f"ver19_noS6_R{R}"] = run([k for k in s19 if k[0]!=6])
    # also drop S98 E1 (post-egress baseline 0.18 T14, comparable to E-1's 0.12 T14)
    out[f"ver19_minusE1_R{R}"] = run([k for k in s19 if k!=(98,1)])
# audit s09 pipe (21) with its own red-noise errors, and with injection SDs
s9 = json.load(open(f"{AUD}/s09_events.json"))["A"]
eph = s9["ephem_pipe"]
rows = {(r["sector"], r["epoch"]): r for r in s9["pipe"]["rows"]}
ep = np.array(eph["epochs"]); sec = np.array([rows[(6 if d.startswith("TESSSPOC") else 32 if "S32" in d else 98, n)]["sector"] for d,n in zip(eph["datasets"],ep)])
t = np.array([rows[(s,n)]["tc_pred"] + rows[(s,n)]["dt_min"]/1440 for s,n in zip(sec,ep)])
et = np.array([max(rows[(s,n)]["e_dt_min"],0.5) for s,n in zip(sec,ep)])/1440
out["audit21_own_errors"] = fit(ep, t, et)
a1 = json.load(open(f"{PUB}/a1_results.json"))["signals"]["A"]["A__pipe"]["injections"]
sdp = {6:a1["S6"]["sd_dt_min"],32:a1["S32"]["sd_dt_min"],98:a1["S98"]["sd_dt_min"]}
sdv = {s: inj[f"S{s}_R60"]["sd"] for s in (6,32,98)}
out["audit21_pubprepINJ"] = fit(ep, t, np.array([sdp[s] for s in sec])/1440)
out["audit21_UverINJ_R60"] = fit(ep, t, np.array([sdv[s] for s in sec])/1440)
m = ~(((sec==98)&(ep==-1))|((sec==32)&(ep==-573)))
out["audit19_UverINJ_R60"] = fit(ep[m], t[m], np.array([sdv[s] for s in sec[m]])/1440)
# pubprep 20
pt = pd.read_csv(f"{PUB}/a1_per_transit.csv"); pa = pt[(pt.signal=="A")&(pt.photometry=="pipe")]
pa = pa[(pa.coverage>=0.7)&(pa.snr_white>=2)&(~pa.edge)]
out["pubprep20"] = fit(pa.epoch, pa.tc_pred+pa.dt_min/1440, pa.sector.map(sdp)/1440)
out["pubprep20_events"] = [f"S{s}:{n}" for s,n in zip(pa.sector,pa.epoch)]
for k,v in out.items():
    if isinstance(v, dict):
        print(f"{k:32s} n={v['n']} P={v['P']:.7f}+-{v['eP']:.2e} T0={v['T0_BJD']:.5f}+-{v['eT0']:.5f} chi2={v['chi2']:.2f}/{v['dof']} p={v['p_sf']:.2f} rms={v['rms_min']:.2f} corr={v['corr']:.2f}")
    else: print(k, v)
json.dump(out, open(f"{HERE}/v1_independent_ephem.json","w"), indent=1, default=float)
