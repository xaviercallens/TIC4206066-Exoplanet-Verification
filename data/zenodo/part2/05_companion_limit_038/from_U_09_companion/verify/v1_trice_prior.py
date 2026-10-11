"""Verify: TRICERATOPS 1.1.0 companion-mass prior and eccentricity prior fractions (runs in TRICE_PY)."""
import numpy as np, json
from scipy.stats import beta
from triceratops.priors import sample_q_companion, sample_ecc
from triceratops.funcs import stellar_relations, flux_relation
np.random.seed(1)
N=400000; Ms=0.606
q=sample_q_companion(np.random.rand(N),Ms); m=q*Ms
out=dict(m_min=float(m.min()),m_max=float(m.max()),frac_gt_038=float(np.mean(m>0.38)),frac_gt_052=float(np.mean(m>0.52)),
         frac_017_038=float(np.mean((m>0.17)&(m<0.38))))
# ecc prior: fraction with e>=emin AND (1+e sin w)/sqrt(1-e^2) large enough to reach rho_app<=rho_max (any w)
e=beta.rvs(0.867,3.03,size=N,random_state=2); w=np.random.uniform(0,2*np.pi,N)
g=(1+e*np.sin(w))/np.sqrt(1-e**2)   # rho_app = rho_true * g^3
for emin in [0.19,0.3,0.5]:
    out[f"P(e>={emin})"]=float(np.mean(e>=emin))
# fraction with g <= x where x = (1-emin)/(1+emin) inverse: e_min at w=270 -> g=(1-e)/sqrt(1-e^2)=sqrt((1-e)/(1+e))
for emin in [0.19,0.3,0.5]:
    gmax=np.sqrt((1-emin)/(1+emin))
    out[f"P(g<={gmax:.3f}) [needs slowdown equiv e_min={emin}]"]=float(np.mean(g<=gmax))
# TRICERATOPS flux ratio for a few masses vs Ms
for mc in [0.17,0.3,0.38,0.52]:
    fr=flux_relation(np.array([mc]))/(flux_relation(np.array([mc]))+flux_relation(np.array([Ms])))
    out[f"trice_fluxfrac_{mc}"]=float(fr[0])
print(json.dumps(out,indent=1)); json.dump(out,open("v1_trice_prior.json","w"),indent=1)
