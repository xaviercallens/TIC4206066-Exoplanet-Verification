import numpy as np, pandas as pd, json
from scipy import stats
import v1_common as V
ev=pd.read_csv("v2_real_events.csv")
out={}
sd={}
for s in (6,32,98):
    d=pd.read_csv(f"v3_inj_S{s}.csv")
    for R in (60,144):
        sel=(d["cov"]>=0.7)&(d[f"dchi{R}"]>=4)&(~d[f"edge{R}"])
        x=d.loc[sel,f"dt{R}"]
        z=x/d.loc[sel,f"ew{R}"]
        sd[(s,R)]=float(x.std())
        out[f"S{s}_R{R}"]=dict(n=int(sel.sum()),sd=float(x.std()),rob=float(1.4826*np.median(np.abs(x-np.median(x)))),mean=float(x.mean()),
            frac_gt60=float((np.abs(x)>60).mean()),kurt=float(stats.kurtosis(x)),z_white_sd=float(z.std()),
            sd_sub150_p2_5=float(np.percentile([np.std(np.random.default_rng(i).choice(x.values,150,replace=False),ddof=1) for i in range(2000)],2.5)))
        # both-sides baseline restricted injections
        sel2=sel&(d.post>=0.3)&(d.pre>=0.3)
        out[f"S{s}_R{R}"]["sd_twosided"]=float(d.loc[sel2,f"dt{R}"].std())
        sel3=sel&(d.post<0.2)|sel&(d.pre<0.2)
        out[f"S{s}_R{R}"]["sd_onesided"]=float(d.loc[sel3,f"dt{R}"].std()); out[f"S{s}_R{R}"]["n_onesided"]=int(sel3.sum())
for k,v in out.items(): print(k,{a:round(b,3) for a,b in v.items()})

def wls(tobs,ep,e):
    A=np.vstack([np.ones_like(ep),ep]).T; Wt=1/e**2
    C=np.linalg.inv(A.T@(A*Wt[:,None])); c=C@(A.T@(Wt*tobs)); r=tobs-A@c
    return c[1],np.sqrt(C[1,1]),c[0],np.sqrt(C[0,0]),float(np.sum((r/e)**2)),r*1440

key=lambda r:(int(r.sector),int(r.epoch))
aud=[(6,-795),(6,-794),(6,-793),(6,-791),(6,-790),(6,-789),(32,-573),(32,-572),(32,-571),(32,-570),(32,-569),(32,-568),(98,0),(98,1),(98,4),(98,5),(98,9),(98,10),(98,13),(98,14)]
pkg=[k for k in aud if k!=(6,-795)]+[(98,-1)]
two=[k for k in aud+[(98,-1)] if not(k in [(32,-573),(98,-1)])]
sets=dict(auditor20=aud,package20=pkg,drop_onesided_573_and_m1=two,aud_minus573=[k for k in aud if k!=(32,-573)])
res={}
for name,S in sets.items():
    e2=ev.set_index(["sector","epoch"]).loc[S].reset_index()
    for R in (60,144):
        tobs=e2.tc+e2[f"dt{R}"]/1440; ep=e2.epoch.values.astype(float)
        e=np.array([sd[(s,R)] for s in e2.sector])/1440
        P,eP,T0,eT0,chi,oc=wls(tobs.values,ep,e)
        ew=np.maximum(e2[f"ew{R}"].values,0.25)/1440
        Pw,ePw,_,_,chiw,_=wls(tobs.values,ep,ew)
        dof=len(S)-2
        res[f"{name}_R{R}"]=dict(n=len(S),P=P,eP=eP,T0=T0+2457000,eT0=eT0,chi2=chi,dof=dof,p_up=float(stats.chi2.sf(chi,dof)),
            chi2_white=chiw,P_white=Pw,eP_white=ePw,rms_oc=float(np.std(oc)),maxz=float(np.max(np.abs(oc/(e*1440)))))
        # MC: draw errors from injection pools of each sector
        rng=np.random.default_rng(1)
        pools={s:pd.read_csv(f"v3_inj_S{s}.csv") for s in (6,32,98)}
        pl={}
        for s in (6,32,98):
            d=pools[s]; sel=(d["cov"]>=0.7)&(d[f"dchi{R}"]>=4)&(~d[f"edge{R}"]); pl[s]=d.loc[sel,f"dt{R}"].values/1440
        Ps=[];chis=[]
        for i in range(5000):
            err=np.array([rng.choice(pl[s]) for s in e2.sector]); tt=T0+P*ep+err
            p_,_,_,_,c_,_=wls(tt,ep,e); Ps.append(p_); chis.append(c_)
        res[f"{name}_R{R}"].update(mc_sdP=float(np.std(Ps)),mc_chi2_frac_below=float(np.mean(np.array(chis)<chi)),mc_chi2_median=float(np.median(chis)))
for k,v in res.items(): print(k,{a:(round(b,8) if isinstance(b,float) else b) for a,b in v.items()})
json.dump(dict(injections=out,ephem=res),open("v4_results.json","w"),indent=1)
