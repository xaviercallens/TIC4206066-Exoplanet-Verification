import numpy as np, json
from astropy.io import fits
D="/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/00_data/tesscut/21x21/"
out={}
for s in ["0006-1-1","0032-1-2","0098-1-2"]:
    with fits.open(D+f"tess-s{s}_80.859830_-1.326960_21x21_astrocut.fits",memmap=True) as h:
        d=h[1].data
        p1=np.array(d["POS_CORR1"],float);p2=np.array(d["POS_CORR2"],float)
        out[s]=dict(n=len(p1),pc1_absmax=float(np.nanmax(np.abs(p1))),pc2_absmax=float(np.nanmax(np.abs(p2))))
json.dump(out,open("/Users/paulsmith/space <3/tess-jev/addenda/tic4206066_FINAL_RECONCILIATION/_work/T3_pixels_localisation/verify/v_poscorr_tesscut.json","w"),indent=1)
print(out)
