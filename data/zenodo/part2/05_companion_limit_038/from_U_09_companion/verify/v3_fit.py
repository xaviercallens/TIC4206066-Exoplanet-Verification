"""Independent re-fit (own code): per-transit local quadratic detrending (not wotan), S6 SPOC-FFI PDCSAP (30 min),
S32 SPOC PDCSAP (2 min), S98 QLP DET_FLUX (not SYS_RM_FLUX), QUALITY==0. Circular transit on a diluted companion host
(TRICERATOPS flux_relation dilution, Mann15 radius). Free k,b; rho fixed (circular) or free (apparent density).
Usage: v3_fit.py A|B"""
import sys, json, numpy as np, batman
from astropy.io import fits
from scipy.optimize import minimize
SIG=sys.argv[1]
D="/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/00_data/mast/"
F={6:("hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits","PDCSAP_FLUX",15),
   32:("tess2020324010417-s0032-0000000004206066-0200-s_lc.fits","PDCSAP_FLUX",1),
   98:("hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits","DET_FLUX",5)}
EPH={"A":(3.1827856,3998.7375,1.87/24),"B":(11.13274,4000.865,1.66/24)}
P,T0,T14=EPH[SIG]; Po,T0o,T14o=EPH["B" if SIG=="A" else "A"]
G=6.67430e-8; RS=1.98841e33/(4/3*np.pi*6.957e10**3)
def aR(rho): return (G*rho*RS*(P*86400)**2/(3*np.pi))**(1/3)
data={}
for s,(fn,col,ss) in F.items():
    h=fits.open(D+fn); d=h[1].data; t=np.array(d["TIME"],float); f=np.array(d[col],float); q=np.array(d["QUALITY"])
    ok=np.isfinite(t)&np.isfinite(f)&(q==0); t,f=t[ok],f[ok]
    ph=(t-T0+0.5*P)%P-0.5*P; pho=(t-T0o+0.5*Po)%Po-0.5*Po
    n=np.round((t-T0)/P); keep_t=[];keep_f=[]
    for e in np.unique(n[np.abs(ph)<0.35]):
        w=(n==e)&(np.abs(ph)<0.35)&(np.abs(pho)>T14o)
        oot=w&(np.abs(ph)>0.75*T14)
        if oot.sum()<10 or (w&(np.abs(ph)<0.4*T14)).sum()<2: continue
        c=np.polyfit(ph[oot],f[oot],2); keep_t.append(ph[w]); keep_f.append(f[w]/np.polyval(c,ph[w]))
    tt=np.concatenate(keep_t); ff=np.concatenate(keep_f)
    oo=np.abs(tt)>0.75*T14; sig=1.4826*np.median(np.abs(ff[oo]-np.median(ff[oo])))
    pm=batman.TransitParams(); pm.t0=0;pm.per=P;pm.rp=.02;pm.a=13;pm.inc=89;pm.ecc=0;pm.w=90;pm.limb_dark="quadratic";pm.u=[.4,.25]
    m=batman.TransitModel(pm,tt,supersample_factor=ss,exp_time=float(h[1].header["TIMEDEL"]))
    data[s]=(tt,ff,sig,pm,m); print(s,len(tt),"sigma ppm",round(sig*1e6))
def chi2(k,a,b,u,fh):
    if not(0<=b<1+k) or k<=0 or k>1 or b/a>=1: return 1e30
    c=0
    for tt,ff,sig,pm,m in data.values():
        pm.rp=k;pm.a=a;pm.inc=np.degrees(np.arccos(b/a));pm.u=list(u)
        c+=np.sum(((ff-(1-fh*(1-m.light_curve(pm))))/sig)**2)
    return c
def fit(rho,u,fh):
    best=(1e31,None); k0=np.sqrt((5e-4 if SIG=="A" else 1e-3)/fh)
    for b0 in (0.05,0.4,0.8):
        x0=[np.log(k0),b0]+([] if rho else [np.log(3.)])
        fun=lambda x: chi2(np.exp(x[0]),aR(rho if rho else np.exp(x[2])),x[1],u,fh)
        r=minimize(fun,x0,method="Nelder-Mead",options=dict(xatol=1e-5,fatol=1e-4,maxiter=3000))
        if r.fun<best[0]: best=(r.fun,r.x)
    return best
H=json.load(open("v2_hosts.json")); u_t=(0.409,0.253)
res={"signal":SIG}
c0,x0=fit(None,u_t,1.0); res["target_free"]=dict(chi2=c0,k=float(np.exp(x0[0])),b=float(x0[1]),rho=float(np.exp(x0[2])))
# rho upper limit, target host
from scipy.optimize import brentq
prof=lambda lr: fit(np.exp(lr),u_t,1.0)[0]-c0-4
try: res["target_rho_max_dchi2_4"]=float(np.exp(brentq(prof,x0[2],x0[2]+1.5,xtol=3e-3)))
except Exception as ex: res["target_rho_max_dchi2_4"]=str(ex)
print(res,flush=True)
rows=[]
Ms=[0.10,0.14,0.17,0.2,0.26,0.3,0.37,0.38,0.44,0.48,0.5,0.52,0.53,0.54,0.56,0.606] if SIG=="A" else [0.10,0.14,0.17,0.26,0.37,0.5,0.606]
for M in Ms:
    h=min(H,key=lambda r:abs(r["M"]-M)); R=h["R_mann"]; rho=M/R**3; fh=h["fluxfrac"]; u=(h["u1"],h["u2"])
    cc,xc=fit(rho,u,fh); ce,xe=fit(None,u,fh)
    row=dict(M=M,R=R,rho=rho,fh=fh,dchi2_circ=cc-c0,b_circ=float(xc[1]),k_circ=float(np.exp(xc[0])),dchi2_free=ce-c0,rho_app=float(np.exp(xe[2])))
    rows.append(row); print({k:round(v,3) for k,v in row.items()},flush=True)
res["hosts"]=rows; json.dump(res,open(f"v3_fit_{SIG}.json","w"),indent=1)
