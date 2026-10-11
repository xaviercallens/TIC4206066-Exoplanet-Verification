"""Host table (TRICE_PY): TRICERATOPS stellar_relations radii/Teff, flux_relation TESS fraction, Claret LD lookup
(same lookup rule as lnZ_STP), plus Mann 2019 mass->M_Ks and Mann 2015 radius (own implementation)."""
import numpy as np, json, pandas as pd, triceratops, os
from triceratops.funcs import stellar_relations, flux_relation
Ms, Rs, Teff = 0.606, 0.623, 4012
A19=np.array([-0.642,-0.208,-8.43e-4,7.87e-3,1.42e-4,-2.13e-4])
mk=np.linspace(4,11,7001); mm=10**sum(a*(mk-7.5)**i for i,a in enumerate(A19))
ld=pd.read_csv(os.path.join(os.path.dirname(triceratops.__file__),"data/ldc_tess.csv")); ld=ld[ld.Z==0.0]
rows=[]
for M in np.round(np.arange(0.10,0.6061,0.02),3).tolist()+[0.17,0.37,0.38,0.53,0.606]:
    R,T=stellar_relations(np.array([M]),np.array([Rs]),np.array([Teff])); R=float(R[0]); T=float(T[0])
    mks=float(np.interp(M,mm[::-1],mk[::-1])); Rm=1.9515-0.3520*mks+0.01680*mks**2
    fr=float((flux_relation(np.array([M]))/(flux_relation(np.array([M]))+flux_relation(np.array([Ms]))))[0])
    logg=np.log10(6.674e-8*M*1.989e33/(Rm*6.957e10)**2)
    Tr=min(max(round(T/250)*250,3500),10000); gr=min(max(round(logg/0.5)*0.5,3.5),5.0)
    r=ld[(ld.Teff==Tr)&(ld.logg==gr)].iloc[0]
    rows.append(dict(M=M,R_trice=R,Teff_trice=T,MKs=mks,R_mann=Rm,fluxfrac=fr,u1=float(r.aLSM),u2=float(r.bLSM)))
json.dump(rows,open("v2_hosts.json","w"),indent=1)
for r in rows: print({k:round(v,3) for k,v in r.items()})
