import numpy as np, json
from astropy.io import fits
D="/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/00_data/mast/"
f=fits.open(D+"hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits")
h=f[1].header; d=f[1].data
print([ (k,h.get(k)) for k in ['TIMESYS','TIMEREF','BJDREFI','BJDREFF','TIMEUNIT','TIMEPIXR','TSTART','TSTOP']])
print(d.columns.names)
q=d['QUALITY'].astype(np.int64); t=d['TIME']
print("N",len(q),"q==0",(q==0).sum())
out={}
for b in range(32):
    n=((q>>b)&1).sum()
    if n: 
        m=((q>>b)&1)==1
        # intervals
        tt=t[m]; br=np.where(np.diff(tt)>0.1)[0]
        segs=[(round(float(tt[0]),2),)] if len(tt) else []
        st=np.r_[0,br+1]; en=np.r_[br,len(tt)-1]
        segs=[(round(float(tt[a]),2),round(float(tt[b_]),2)) for a,b_ in zip(st,en)]
        out[f"2^{b}"]=dict(n=int(n),segs=segs[:10])
        print(f"2^{b}",n,segs[:8])
json.dump(out,open("v01_qlp_bits.json","w"),indent=1)
# header comments on quality
for k in h: 
    if 'QUAL' in k or 'FLAG' in k.upper(): print(k,h[k])
print(repr(f[0].header)[:3000])
