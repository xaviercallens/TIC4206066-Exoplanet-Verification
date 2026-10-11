"""Verifier: build an own TRICERATOPS 1.1.0 target (network: TIC, TESScut, Gaia DR3 background). trice_env, cwd=verify/."""
import os, pickle, time, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join("/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad", "v08_cache")
os.makedirs(CACHE, exist_ok=True)
import triceratops.triceratops as tt
t0 = time.time()
T = tt.target(ID=4206066, sectors=np.array([32, 98]), search_radius=10, lightkurve_cache_dir=CACHE)
print("init %.0f s; %d stars; bkg=%s" % (time.time() - t0, len(T.stars), T.trilegal_fname), flush=True)
try:
    aps = T.get_spoc_apertures()
    print("spoc apertures:", [None if a is None else len(a) for a in aps], flush=True)
    T._v_spoc_aps = aps
except Exception as ex:
    print("get_spoc_apertures failed:", repr(ex), flush=True); T._v_spoc_aps = None
pickle.dump(T, open(os.path.join(H, "target_v.pkl"), "wb"))
print(T.stars[["ID", "Tmag", "sep (arcsec)", "mass", "rad", "Teff", "plx"]].head(8).to_string())
print("pix_coords:", [np.round(p[0], 2).tolist() for p in T.pix_coords])
