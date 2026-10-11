"""Single-star SED fit of TIC 4206066 with astroARIADNE (Vines & Jenkins 2022, MNRAS 513, 2719), BMA over the
model grids recommended in the package README. Photometry from ARIADNE's own Librarian (catalogue lookup by Gaia id).
Av fixed to 0 (d = 35.6 pc; README: "Fix Av to zero (nearby stars)"). Run in a separate env (not the project venv).

Usage: python sed_fit.py <out_folder> [nlive] [threads]
"""
import sys, json, pickle
from pathlib import Path
from astroARIADNE.star import Star
from astroARIADNE.fitter import Fitter

out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
nlive = int(sys.argv[2]) if len(sys.argv) > 2 else 500
threads = int(sys.argv[3]) if len(sys.argv) > 3 else 8

s = Star('TIC4206066', 80.859827, -1.326964, g_id=3220388198192519424, Av=0.0, verbose=True)
s.print_mags()
# Removed: SDSS and PS1 are saturated for G = 10.5 (SDSS saturates at r ~ 14, sdss4.org DR17 imaging data quality;
# PS1 at g,r,i ~ 13.5, z ~ 13, y ~ 12, Tonry et al. 2012 ApJ 750, 99). TIC Tmag is computed from Gaia photometry
# (Stassun et al. 2019 AJ 158, 138), so it is not an independent measurement.
REMOVE = ['SDSS_u', 'SDSS_g', 'SDSS_r', 'SDSS_i', 'SDSS_z', 'PS1_g', 'PS1_r', 'PS1_i', 'PS1_z', 'PS1_y', 'TESS']
for k in REMOVE:
    if k in s.filter_names[s.filter_mask]:
        s.remove_mag(k)
s.print_mags()
mags = {str(n): [float(m), float(e)] for n, m, e, u in zip(s.filter_names, s.mags, s.mag_errs, s.used_filters) if u == 1}
json.dump({"mags_used": mags, "plx": s.plx, "plx_e": s.plx_e, "dist": s.dist, "dist_e": s.dist_e,
           "spectroscopic_params": str(s.spectroscopic_params)}, open(out / "star_inputs.json", "w"), indent=1, default=float)

f = Fitter()
f.star = s
f.setup = ['dynesty', nlive, 0.5, 'multi', 'rwalk', threads, False]
f.av_law = 'fitzpatrick'
f.out_folder = str(out)
f.bma = True
f.models = ['phoenix', 'btsettl', 'btnextgen', 'btcond', 'kurucz', 'ck04', 'bosz']
f.n_samples = 100000
f.estimate_age = False
f.prior_setup = {'teff': ('default'), 'logg': ('default'), 'z': ('default'),
                 'dist': ('default'), 'rad': ('default'), 'Av': ('fixed', 0)}
f.initialize()
f.fit_bma()
print("prior_sources", getattr(f, "prior_sources", None))
