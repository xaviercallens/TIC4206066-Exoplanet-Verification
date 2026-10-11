import numpy as np, vcommon as v
d,t,q=v.load("qlp98")
det=np.asarray(d["DET_FLUX"],float); s=np.asarray(d["SYS_RM_FLUX"],float)
ok=(q==0)&np.isfinite(det)
tt=t[ok]; dd=det[ok]; ss=s[ok]
gi=np.where(np.diff(tt)>0.5)[0]; ed=np.r_[0,gi+1,len(tt)]
for a,b in zip(ed[:-1],ed[1:]):
    print("seg %.2f-%.2f n=%d med DET %.5f rstd DET %.0f ppm | rstd SYS %.0f ppm; p2p DET %.0f"%(tt[a],tt[b-1],b-a,np.median(dd[a:b]),v.rstd(dd[a:b])*1e6,v.rstd(ss[a:b]/np.median(ss[a:b]))*1e6, v.rstd(np.diff(dd[a:b]))/np.sqrt(2)*1e6))
