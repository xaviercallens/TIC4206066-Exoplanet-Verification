import numpy as np, pandas as pd, sys
import v1_common as V
s=int(sys.argv[1]); N=int(sys.argv[2]); rng=np.random.default_rng(7000+s)
t,f,cad,sh=V.load(s); tpl=V.model_fn(cad); Dinj=5.2e-4
nA=np.arange(-900,40); tA=V.T0A+nA*V.PA
rows=[]; tries=0
while len(rows)<N and tries<50*N:
    tries+=1; ti=rng.uniform(t.min(),t.max())
    if np.min(np.abs(ti-tA))<V.W+V.T14A: continue
    tt,yy=V.window(t,f,ti)
    if len(tt)<10: continue
    pre,post,cov=V.baseline(tt,ti); covf=cov*cad/V.T14A
    if covf<0.5: continue
    yy=yy*(1-Dinj*tpl(tt-ti))
    dts,chis,Ds,cn=V.scan(tt,yy,ti,tpl)
    s2=np.min(chis)/(len(tt)-5); r=dict(t_inj=ti,cov=covf,pre=pre*cad/V.T14A,post=post*cad/V.T14A,n=len(tt))
    for R in (60,144):
        dt,Dd,j,edge=V.best(dts,chis,Ds,R)
        m=np.abs(dts)<=R; ok=chis[m]/s2<=chis[j]/s2+1; ew=(dts[m][ok].max()-dts[m][ok].min())/2
        r.update({f"dt{R}":dt,f"D{R}":Dd,f"edge{R}":edge,f"ew{R}":ew,f"dchi{R}":(cn-chis[j])/s2})
    rows.append(r)
df=pd.DataFrame(rows); df.to_csv(f"v3_inj_S{s}.csv",index=False); print(s,len(df),tries)
