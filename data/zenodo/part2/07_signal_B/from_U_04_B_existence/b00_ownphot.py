"""Own TESScut aperture photometry (two aperture sizes) for S6, S32, S98; cached to out/ownphot_s{S}.npz."""
import numpy as np, json, sys
import common_b as c
res = {}
for s in [6, 32, 98]:
    for name, thr in [("ap12", 0.12), ("ap05", 0.05)]:
        p = c.tesscut_phot(s, ap_thresh=thr)
        np.savez(f"out/ownphot_s{s}_{name}.npz", t=p["t"], flux=p["flux"], bkg=p["bkg"], cx=p["cx"], cy=p["cy"], q=p["q"], ap=p["ap"])
        res[f"s{s}_{name}"] = dict(peak=p["peak"], n_ap=p["n_ap"], ap_pixels=np.argwhere(p["ap"]).tolist())
        print(s, name, p["peak"], p["n_ap"], np.nanmedian(p["flux"]), flush=True)
json.dump(res, open("out/ownphot_apertures.json", "w"), indent=1)
