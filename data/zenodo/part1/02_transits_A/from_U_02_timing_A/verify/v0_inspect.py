import numpy as np
from astropy.io import fits
D="/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/00_data/mast/"
h=fits.open(D+"hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits")
d=h[1].data; print(h[1].columns.names)
t=d['TIME']; q=d['QUALITY'].astype(np.int64); f=d['SYS_RM_FLUX']
tc=3995.5540; T14=1.887/24
for lo,hi in [(tc-0.4,tc+0.4)]:
    m=(t>lo)&(t<hi)
    bits=np.unique(q[m]); print("bits in window",bits)
    # summarize by contiguous runs of quality
    tt=t[m]; qq=q[m]; ff=f[m]
    run_start=0
    for i in range(1,len(tt)+1):
        if i==len(tt) or qq[i]!=qq[run_start] or (tt[i]-tt[i-1])>0.01:
            print(f"  {tt[run_start]:.4f}-{tt[i-1]:.4f} (dt from tc {(tt[run_start]-tc)*24:+.2f}h..{(tt[i-1]-tc)*24:+.2f}h) q={qq[run_start]} n={i-run_start} finite={np.isfinite(ff[run_start:i]).sum()}")
            run_start=i
print("egress at +%.2f h"%(T14/2*24))
# B times
PB,T0B=11.13274,4000.865
for s,tcA in [(32,2175.0008)]:
    pass
nB=np.round((2175.0008-T0B)/PB); print("B near S32 E-573:",T0B+nB*PB, (T0B+nB*PB-2175.0008)*24*60,"min")
