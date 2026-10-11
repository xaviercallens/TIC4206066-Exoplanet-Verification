"""Per-band residuals of the ARIADNE BMA fit (sed_out/BMA.pkl): observed flux vs model flux at the median parameters,
for each fitted grid and for the BMA-averaged parameters (evaluated on the highest-weight grid). Residual in mag =
-2.5 log10(obs / model); sigma uses the catalogue error only (the fitted per-band noise terms are listed separately).
An unresolved companion shows up as a positive (too bright) residual growing towards the red / infrared.
"""
import sys, json, pickle
from pathlib import Path
import numpy as np
from astroARIADNE.fitter import Fitter
from astroARIADNE.sed_library import model_grid

HERE = Path(__file__).resolve().parent
d = pickle.load(open(HERE / "sed_out/BMA.pkl", "rb"))
s = d["star"]; m = s.filter_mask
filts, wave, obs, err = s.filter_names[m], s.wave[m], s.flux[m], s.flux_er[m]
f = Fitter(); f.star = s
res = {"filters": list(map(str, filts)), "wave_um": wave.tolist(), "obs_flux": obs.tolist(), "obs_err": err.tolist(),
       "weights": {k: float(v) for k, v in d["weights"].items()}, "per_model": {}, "noise_terms_bma": {}}
for k in d["weighted_average"]:
    if "noise" in str(k): res["noise_terms_bma"][str(k)] = float(np.median(d["weighted_average"][k]))
best = max(d["weights"], key=d["weights"].get)
def run(label, grid, pars):
    intp = f.load_interpolator(grid)
    theta = [pars["teff"], pars["logg"], pars["z"], pars["dist"], pars["rad"], 0.0]
    model = model_grid(theta, filts, wave, intp, False, d["av_law"])
    rmag = -2.5 * np.log10(obs / model)
    nsig = (obs - model) / err
    res["per_model"][label] = {"grid": grid, "params": {a: float(b) for a, b in pars.items()},
                               "resid_mag": rmag.tolist(), "resid_sigma": nsig.tolist(),
                               "chi2_catalog_errors": float(np.sum(nsig ** 2)), "n": int(len(obs))}
    print(f"{label:>12} teff {pars['teff']:.0f} rad {pars['rad']:.3f}  resid(mag): " +
          " ".join(f"{a}:{b:+.3f}" for a, b in zip(filts, rmag)))
for g in d["originals"]:
    o = d["originals"][g]
    run(g, g, {p: float(np.median(o[p])) for p in ["teff", "logg", "z", "dist", "rad"]})
run("BMA", best, {p: float(np.median(d["weighted_average"][p])) for p in ["teff", "logg", "z", "dist", "rad"]})
json.dump(res, open(HERE / "out/sed_residuals.json", "w"), indent=1)

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True)
b = res["per_model"]["BMA"]
lam = np.array(wave)
ax[0].errorbar(lam, lam * obs, yerr=lam * err, fmt="o", color="#c0392b", label="observed (λF_λ)")
intp = f.load_interpolator(best); pb = b["params"]
ax[0].plot(lam, lam * model_grid([pb["teff"], pb["logg"], pb["z"], pb["dist"], pb["rad"], 0.0], filts, wave, intp, False, d["av_law"]),
           "s", mfc="none", mec="#3b6ea5", ms=9, label=f"model ({best}, BMA medians)")
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_ylabel("λF_λ (erg s⁻¹ cm⁻²)"); ax[0].legend()
ax[0].set_title(f"TIC 4206066 SED: Teff {pb['teff']:.0f} K, R {pb['rad']:.3f} R☉, Av = 0")
for lab, c in [("phoenix", "#7f8c8d"), ("btsettl", "#3b6ea5"), ("btnextgen", "#27ae60"), ("btcond", "#8e44ad"), ("BMA", "#c0392b")]:
    if lab in res["per_model"]:
        ax[1].plot(lam, res["per_model"][lab]["resid_mag"], "o-", color=c, label=lab, ms=4 if lab != "BMA" else 7)
ax[1].axhline(0, color="k", lw=0.8); ax[1].set_xscale("log"); ax[1].set_xlabel("wavelength (µm)")
ax[1].set_ylabel("obs − model (mag)\n+ = observed brighter"); ax[1].legend(fontsize=8, ncol=5)
for x, n in zip(lam, filts): ax[1].annotate(str(n).replace("GaiaDR2v2_", "").replace("_MvB", "").replace("2MASS_", "").replace("WISE_RSR_", ""), (x, ax[1].get_ylim()[0]), fontsize=7, rotation=90, va="bottom")
fig.tight_layout(); fig.savefig(HERE / "plots/sed_residuals.png", dpi=130)
