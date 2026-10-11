# Independent coarse global fit: trapezoid template, per-window linear baseline, grid over (T0, P).
import numpy as np, json
from astropy.io import fits
D="/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/00_data/mast/"
def load(fn,col,cad):
    d=fits.open(D+fn)[1].data; t=np.array(d['TIME'],float); f=np.array(d[col],float); q=np.array(d['QUALITY'])
    m=np.isfinite(t)&np.isfinite(f)&(q==0); t,f=t[m],f[m]; f=f/np.nanmedian(f); return t,f,cad
S={6:load("hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits","PDCSAP_FLUX",1800/86400),
   32:load("tess2020324010417-s0032-0000000004206066-0200-s_lc.fits","PDCSAP_FLUX",120/86400),
   98:load("hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits","SYS_RM_FLUX",200/86400)}
P0,T00=3.1827856,3998.7375; PB,TB=11.13274,4000.865
T14=1.87/24; Tin=0.12/24; HW=0.35
def trap(dt):
    a=np.abs(dt); y=np.clip((T14/2-a)/Tin,0,1); return y
wins=[]
for s,(t,f,cad) in S.items():
    # B mask
    ph=((t-TB+PB/2)%PB)-PB/2; keep=np.abs(ph)>1.0/24
    t,f=t[keep],f[keep]
    n0=int(np.floor((t.min()-T00)/P0)); n1=int(np.ceil((t.max()-T00)/P0))
    sd=np.std(np.diff(f))/np.sqrt(2)
    for n in range(n0,n1+1):
        tc=T00+n*P0; w=np.abs(t-tc)<HW
        inn=np.abs(t-tc)<T14/2
        if inn.sum()*cad<0.5*T14 or w.sum()<10: continue
        # need baseline both sides
        if (t[w]<tc-T14/2).sum()<2 or (t[w]>tc+T14/2).sum()<2: continue
        wins.append(dict(s=s,n=n,t=t[w],f=f[w],cad=cad,sig=sd))
print("windows",len(wins),{s:sum(1 for w in wins if w['s']==s) for s in S})
ss=np.linspace(-0.5,0.5,15)
def model_win(w,tc):
    if w['cad']>0.01:
        dt=(w['t'][:,None]+ss[None,:]*w['cad'])-tc; return trap(dt).mean(1)
    return trap(w['t']-tc)
def chi2(T0,P,return_depth=False):
    # linear in (depth, a_i, b_i): solve jointly -> profile depth by 1D scan-free approach: solve each window with depth common via stacking
    # Build normal equations
    rows=[];
    tot=0; A_list=[]; 
    # two-stage: for given depth, baselines are linear; solve depth analytically via projection
    num=0; den=0; res0=0
    for w in wins:
        tc=T0+w['n']*P; m=model_win(w,tc); X=np.vstack([np.ones_like(w['t']),w['t']-tc]).T
        # project out baseline
        Q,_=np.linalg.qr(X); pf=w['f']-Q@(Q.T@w['f']); pm=m-Q@(Q.T@m)
        wt=1/w['sig']**2; num+=wt*(pf@pm); den+=wt*(pm@pm); res0+=wt*(pf@pf)
    dep=-num/den
    c=res0-num**2/den
    return (c,dep) if return_depth else c
dT=np.arange(-12,12.01,0.5)/1440; dP=np.arange(-20e-6,20.01e-6,1e-6)
G=np.array([[chi2(T00+a,P0+b) for b in dP] for a in dT])
i,j=np.unravel_index(G.argmin(),G.shape)
c,dep=chi2(T00+dT[i],P0+dP[j],True)
# reduced chi2 rescale
npts=sum(len(w['t']) for w in wins); red=c/(npts-2*len(wins)-3)
Gs=(G-G.min())/red
# marginal 1-sigma from dchi2=1 profile
profT=Gs.min(1); profP=Gs.min(0)
def ci(x,p):
    ok=x[p<=1]; return float(ok.min()),float(ok.max())
out=dict(n_windows=len(wins),per_sector={s:sum(1 for w in wins if w['s']==s) for s in S},
 T0_best=float(T00+dT[i]),P_best=float(P0+dP[j]),depth_ppm=float(dep*1e6),red_chi2=float(red),
 T0_1sig=[T00+v for v in ci(dT,profT)],P_1sig=[P0+v for v in ci(dP,profP)],
 dT0_min=float(dT[i]*1440),dP=float(dP[j]))
print(json.dumps(out,indent=1))
json.dump(out,open('v03_grid.json','w'),indent=1)
# per-sector time shift at manuscript ephemeris
sh={}
for s in S:
    ws=[w for w in wins if w['s']==s]; best=[]
    for a in np.arange(-30,30.1,0.5)/1440:
        num=den=res0=0
        for w in ws:
            tc=T00+w['n']*P0+a; m=model_win(w,tc); X=np.vstack([np.ones_like(w['t']),w['t']-tc]).T
            Q,_=np.linalg.qr(X); pf=w['f']-Q@(Q.T@w['f']); pm=m-Q@(Q.T@m)
            num+=pf@pm; den+=pm@pm; res0+=pf@pf
        best.append((res0-num**2/den)/ws[0]['sig']**2)
    best=np.array(best); aa=np.arange(-30,30.1,0.5); k=best.argmin(); ok=aa[best-best.min()<=1]
    sh[s]=dict(shift_min=float(aa[k]),lo=float(ok.min()),hi=float(ok.max()))
print(sh); json.dump(sh,open('v03_sector_shifts.json','w'),indent=1)
