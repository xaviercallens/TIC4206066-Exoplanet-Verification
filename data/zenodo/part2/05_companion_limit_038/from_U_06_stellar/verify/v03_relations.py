import numpy as np, json
from astropy.table import Table
g=Table.read("../data/gaia_dr3_source.ecsv")[0]
plx=float(g['parallax']); eplx=float(g['parallax_error'])
print('plx',plx,eplx,'ruwe',float(g['ruwe']),'gmag',float(g['phot_g_mean_mag']),'bp_rp',float(g['bp_rp']))
try:
    from zero_point import zpt; zpt.load_tables()
    zp=float(zpt.get_zpt(np.array([float(g['phot_g_mean_mag'])]),np.array([float(g['nu_eff_used_in_astrometry'])]),np.array([np.nan_to_num(float(g['pseudocolour']))]),np.array([float(g['ecl_lat'])]),np.array([int(g['astrometric_params_solved'])]))[0])
except Exception as e: zp=None; print('zp fail',e)
print('zp',zp)
K,eK=7.694,0.029
rng=np.random.default_rng(7); N=1_000_000
p=rng.normal(plx-(zp or 0),eplx,N); k=rng.normal(K,eK,N)
MK=k+5*np.log10(p/1000)+5
R=(1.9515-0.3520*MK+0.01680*MK**2)*(1+rng.normal(0,0.0289,N))
a=[-0.642,-0.208,-8.43e-4,7.87e-3,1.42e-4,-2.13e-4]
x=MK-7.5; M=10**sum(a[i]*x**i for i in range(6))*(1+rng.normal(0,0.024,N))
rho=M/R**3
pc=lambda v:[float(q) for q in np.percentile(v,[16,50,84])]
out=dict(zp=zp,d_raw=1000/plx,d_zp=1000/(plx-(zp or 0)),MK=pc(MK),R=pc(R),M=pc(M),rho=pc(rho),MG=float(g['phot_g_mean_mag'])+5*np.log10((plx-(zp or 0))/1000)+5)
# host threshold rho=3.5
mk=np.linspace(4.5,7,25001); x=mk-7.5
Rm=1.9515-0.3520*mk+0.01680*mk**2; Mm=10**sum(a[i]*x**i for i in range(6))
r=Mm/Rm**3; j=np.argmin(abs(r-3.5)); out['thr']=dict(MK=float(mk[j]),M=float(Mm[j]),R=float(Rm[j]))
print(json.dumps(out,indent=1)); json.dump(out,open('v03_relations.json','w'),indent=1)
