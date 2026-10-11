import numpy as np, pandas as pd, sys
import v1_common as V
rows=[]
for s in (6,32,98):
    t,f,cad,sh = V.load(s); tpl=V.model_fn(cad)
    n0=int(np.floor((t.min()-V.T0A)/V.PA)); n1=int(np.ceil((t.max()-V.T0A)/V.PA))
    for n in range(n0,n1+1):
        tc=V.T0A+n*V.PA
        tt,yy=V.window(t,f,tc)
        if len(tt)<10: continue
        pre,post,cov=V.baseline(tt,tc); covf=cov*cad/V.T14A
        if covf<0.5: continue
        dts,chis,Ds,cn=V.scan(tt,yy,tc,tpl)
        s2=np.min(chis)/(len(tt)-5)
        r={}
        for R in (60,144):
            dt,Dd,j,edge=V.best(dts,chis,Ds,R)
            # white err: delta chi2=1 half width
            m=np.abs(dts)<=R; ok=chis[m]/s2<=chis[j]/s2+1; ew=(dts[m][ok].max()-dts[m][ok].min())/2
            r.update({f"dt{R}":dt,f"D{R}":Dd,f"edge{R}":edge,f"ew{R}":ew,f"dchi_null{R}":(cn-chis[j])/s2})
        # also unconstrained-sign fit at R144
        dtsU,chU,DU,_=V.scan(tt,yy,tc,tpl,dpos=False); jU=np.argmin(chU)
        rows.append(dict(sector=s,epoch=n,tc=tc,cov=covf,n=len(tt),pre=pre*cad/V.T14A,post=post*cad/V.T14A,dtU=dtsU[jU],DU=DU[jU],**r))
df=pd.DataFrame(rows); pd.set_option("display.width",250)
df.to_csv("v2_real_events.csv",index=False)
print(df.round(3).to_string())
