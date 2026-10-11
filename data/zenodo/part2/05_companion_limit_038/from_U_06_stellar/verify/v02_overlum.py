import numpy as np, pandas as pd, json
D="/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/06_stellar_dl/overluminous_catalog.csv"
c=pd.read_csv(D)
t=c[c.source_id==3220388198192519424].iloc[0]; print(t.to_dict())
O=c.pabsG_dr-c.absG_dr; Ot=t.pabsG_dr-t.absG_dr
res={'O_target':float(Ot)}
for dM,dc in [(0.1,0.02),(0.2,0.03),(0.3,0.05),(0.5,0.1)]:
    s=(abs(c.pabsG_dr-t.pabsG_dr)<dM)&(abs(c.pg_rp_dr-t.pg_rp_dr)<dc)
    x=O[s]; xn=x[~c.ol[s]]
    mad=1.4826*np.median(abs(x-np.median(x)))
    p=np.percentile(x,[16,50,84])
    res[f'win_{dM}_{dc}']=dict(n=int(s.sum()),med=float(np.median(x)),mad=float(mad),p16_84_half=float((p[2]-p[0])/2),std_normolum=float(np.std(xn)),frac_above=float((x>Ot).mean()))
# companion mass limit: Pecaut-Mamajek parsed table from auditor used ONLY as table of (M_G? ) -> read raw EEM file myself
import re
rows=[]
for line in open("/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/refs/EEM_dwarf_UBVIJHK_colors_Teff.txt"):
    if line.startswith('#SpT'): hdr=line[1:].split(); started=True; continue
    if 'hdr' in dir() and re.match(r'^[OBAFGKMLTY]\d',line):
        f=line.split()
        if len(f)==len(hdr): rows.append(f)
pm=pd.DataFrame(rows,columns=hdr)
print(hdr)
pm=pm[['SpT','M_G','Msun']].replace('...',np.nan).dropna()
pm['M_G']=pm.M_G.astype(float); pm['Msun']=pm.Msun.astype(float)
pm=pm[(pm.M_G>6)&(pm.M_G<16)].sort_values('M_G')
MG1=float(t.absG_dr)
lims={}
for sig in [0.052,0.068,0.10]:
    for k in [2,3]:
        Olim=Ot+k*sig
        # total flux brighter than predicted single by Olim: F_tot/F1 = 10^(0.4 Olim) assuming primary at predicted mag
        f=10**(0.4*Olim)-1  # companion/primary flux ratio
        MG2=float(t.pabsG_dr)-2.5*np.log10(f)
        m=np.interp(MG2,pm.M_G.values,pm.Msun.values)
        lims[f'{k}sig_{sig}']=dict(Olim=Olim,f=f,MG2=MG2,mass=float(m))
res['limits']=lims
print(json.dumps(res,indent=1))
json.dump(res,open('v02_overlum.json','w'),indent=1)
