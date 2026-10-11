import numpy as np, json, pandas as pd
from scipy import stats
U="/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/"
out={}
def comb(sets):
    w=np.array([[1/s[1][0]**2,1/s[1][1]**2] for s in sets]); x=np.array([s[0] for s in sets])
    m=(w*x).sum(0)/w.sum(0); e=1/np.sqrt(w.sum(0)); chi2=((m/e)**2).sum()
    return dict(mean=m.round(3).tolist(),sig=e.round(3).tolist(),zE=round(m[0]/e[0],3),chi2=round(chi2,3),p=round(float(stats.chi2.sf(chi2,2)),4))
S32PC=([5.2055,0.3928],[3.0811,2.7826]); S98U=([5.6489,-0.3163],[5.2686,4.5166])
out["A_headline_recomputed"]=comb([S32PC,S98U])
out["B_S98only"]=comb([([1.1905,-1.0173],[2.2947,3.1164])])
# capable-star and 11" star sigma distances per product epoch, from (+5,0)+-3 and exact
pp=pd.read_csv(U+"03_pixels/neighbour_census_per_product.csv")
res={}
for tic in [712986249.0,4206077.0,4206064.0]:
    r=pp[pp.TIC==tic][["product","sep","dE","dN"]]
    rows=[]
    for _,q in r.iterrows():
        rows.append(dict(product=q["product"],sep=round(q.sep,2),
          sig_pm3=round(float(np.hypot(q.dE-5,q.dN-0)/3),2),
          sig_exact=round(float(np.hypot((q.dE-5.3185)/2.6597,(q.dN-0.1977)/2.3691)),2)))
    res[str(int(tic))]=rows
out["neighbour_sigma_by_product_epoch"]=res
json.dump(out,open(__file__.replace(".py",".json"),"w"),indent=1)
print(json.dumps(out,indent=1))
