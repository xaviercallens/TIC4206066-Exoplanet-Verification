import pickle, numpy as np, json
p="/Users/paulsmith/space <3/tess-jev/addenda/companion_check_tic4206066/sed_out/BMA.pkl"
class U(pickle.Unpickler):
    def find_class(self,m,n):
        try: return super().find_class(m,n)
        except Exception: return type(n,(),{'__setstate__':lambda self,st: self.__dict__.update(st if isinstance(st,dict) else {'s':st})})
d=U(open(p,'rb')).load()
print(d.keys())
w=d['weights']; print('weights',w)
rng=np.random.default_rng(1)
res={}
for k in ['teff','rad','norm','lum']:
    if k not in d['originals'][list(w)[0]] and k not in d.get('weighted_average',{}): continue
    out={}
    if k in d['weighted_average']:
        wa=np.asarray(d['weighted_average'][k]); ws=np.asarray(d['weighted_samples'][k])
        out['wavg']=list(np.percentile(wa,[16,50,84])); out['wsamp']=list(np.percentile(ws,[16,50,84]))
    per={}
    mix=[]
    for g,wt in w.items():
        o=d['originals'][g]
        if k in o:
            x=np.asarray(o[k]); per[g]=[float(np.median(x)),float(np.std(x))]
            mix.append(rng.choice(x,int(round(wt*400000))))
    if mix:
        m=np.concatenate(mix); out['mixture']=list(np.percentile(m,[16,50,84]))
    out['per_grid']=per
    res[k]=out
    print(k,out)
print({g:list(d['originals'][g].keys()) for g in list(w)[:1]})
json.dump(res,open('v01_ariadne.json','w'),indent=1,default=float)
