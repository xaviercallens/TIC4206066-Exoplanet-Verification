"""C021: classify the 16 non-target Gaia DR3 sources within 60 arcsec (read-only raw cone file)."""
import json, os
A = "/Users/paulsmith/space <3/tess-jev/addenda/archival_extra_tic4206066/highres_imaging/raw/gaia_dr3_cone_60arcsec.json"
f = lambda x: None if x in (None, "", "--", "nan") else float(x)
rows = json.load(open(A)); tgt = [r for r in rows if f(r["sep_arcsec"]) < 0.1][0]
plx_t, G_t = f(tgt["parallax"]), f(tgt["phot_g_mean_mag"])
out = []
for r in sorted(rows, key=lambda r: f(r["sep_arcsec"])):
    if r is tgt: continue
    plx, pe, G = f(r["parallax"]), f(r["parallax_error"]), f(r["phot_g_mean_mag"])
    out.append(dict(source_id=r["source_id"], sep=round(f(r["sep_arcsec"]), 2), G=round(G, 2),
                    solution="5/6-param" if plx is not None else "2-param",
                    n_sigma_from_target_plx=None if plx is None else round((plx_t - plx) / pe, 1),
                    max_depth_ppm_full_eclipse_Gband=round(1e6 * 10 ** (-0.4 * (G - G_t)), 1)))
summ = dict(n=len(out), n_2param=sum(o["solution"] == "2-param" for o in out),
            min_nsigma_plx=min(o["n_sigma_from_target_plx"] for o in out if o["n_sigma_from_target_plx"] is not None),
            max_depth_2param_ppm=max(o["max_depth_ppm_full_eclipse_Gband"] for o in out if o["solution"] == "2-param"))
json.dump(dict(summary=summ, rows=out), open(os.path.dirname(os.path.abspath(__file__)) + "/t4_04_gaia60_census.json", "w"), indent=1)
print(json.dumps(summ, indent=1))
