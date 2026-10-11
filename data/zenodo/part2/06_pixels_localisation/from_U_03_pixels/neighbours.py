"""G4: which Gaia DR3 (and TIC-only) sources within 180" can produce the observed depths?

For every light-curve product used in the manuscript the actual photometric aperture is used:
  S6  TESS-SPOC 30-min  (APERTURE extension of the TESS-SPOC TP file, bit 2 = optimal aperture)
  S32 SPOC 2-min        (APERTURE extension of the SPOC TPF)
  S32 TESS-SPOC 10-min  (APERTURE extension)       [not used in the manuscript; for completeness]
  S98 QLP 200-s         (BESTAP = '3x3 square' in the QLP header; centred on the target pixel of the FFI grid;
                         QLP >= S94 uses TGLC decontaminated photometry, i.e. neighbours' static flux removed)
Each star's flux fraction inside the aperture comes from the TESS PRF (MAST start_s0004) at its Gaia DR3 position
propagated to the sector epoch, after a global WCS shift fitted with the out-of-transit scene model (G2).
TESS magnitudes: TIC 8.2 (Gaia id match); otherwise Stassun et al. 2019 eqs. 1-2.
Required eclipse depth of star j:  d_j = delta_obs * F_ref / F_j,ap  with F_ref = CROWDSAP * F_tot,ap for PDCSAP
(crowding-corrected) and F_ref = F_target,ap for the decontaminated QLP/TGLC photometry.
"""
from pixloc import *
import pandas as pd

DEPTHS = {"A": [500e-6, 400e-6], "B": [1000e-6, 800e-6]}   # nominal, conservative-low


def aperture_products():
    prods = []
    for key, kind in [("S6_tspoc_tp", "pdcsap"), ("S32_tpf", "pdcsap"), ("S32_tspoc_tp", "pdcsap")]:
        ds = Dataset(key)
        ap = (ds.cube.aperture & 2) > 0
        with fits.open(FILES[key]) as h:
            pass
        lc_key = {"S6_tspoc_tp": "S6_tspoc_lc", "S32_tpf": "S32_spoc_lc", "S32_tspoc_tp": "S32_tspoc_lc"}[key]
        hl = fits.getheader(FILES[lc_key], 1)
        prods.append(dict(name=key, ds=ds, ap=ap, kind=kind, crowdsap=float(hl["CROWDSAP"]),
                          flfrcsap=float(hl["FLFRCSAP"])))
    ds = Dataset("S98_ffi", "clean")
    sf = ds.scene_fit()
    for size, label in [(3, "QLP98_3x3_best"), (1, "QLP98_1x1_small"), (5, "QLP98_5x5_large")]:
        ap = np.zeros((ds.cube.ny, ds.cube.nx), bool)
        xc, yc = int(round(sf["x_t"])), int(round(sf["y_t"]))
        h = size // 2
        ap[yc - h:yc + h + 1, xc - h:xc + h + 1] = True
        prods.append(dict(name=label, ds=ds, ap=ap, kind="decontaminated", crowdsap=None, flfrcsap=None))
    return prods


def census():
    rows = []
    summary = {}
    for p in aperture_products():
        ds = p["ds"]
        sf = ds.scene_fit()
        shape = (ds.cube.ny, ds.cube.nx)
        T = ds.gaia["Tmag_used"].values
        flux = 10 ** (-0.4 * (T - T[ds.it]))
        frac = np.zeros(len(T))
        for j in range(len(T)):
            xj, yj = ds.g_x[j] + sf["dx"], ds.g_y[j] + sf["dy"]
            if xj < -8 or yj < -8 or xj > shape[1] + 7 or yj > shape[0] + 7:
                continue
            frac[j] = (ds.prf.image(xj, yj, shape) * p["ap"]).sum()
        Fap = flux * frac
        Ftot = Fap.sum()
        crowd_model = Fap[ds.it] / Ftot
        if p["kind"] == "pdcsap":
            Fref = p["crowdsap"] * Ftot
        else:
            Fref = Fap[ds.it]
        summary[p["name"]] = dict(n_ap_pix=int(p["ap"].sum()), target_frac_in_ap=float(frac[ds.it]),
                                  crowding_model=float(crowd_model), CROWDSAP_header=p["crowdsap"],
                                  FLFRCSAP_header=p["flfrcsap"], scene_shift_px=[sf["dx"], sf["dy"]],
                                  F_tot_over_F_target=float(Ftot / Fap[ds.it]))
        for j in range(len(T)):
            if j == ds.it:
                continue
            e_, n_ = offsets_arcsec(ds.g_ra[ds.it], ds.g_dec[ds.it], ds.g_ra[j], ds.g_dec[j])
            r = float(np.hypot(e_[0], n_[0]))
            if r > 180:
                continue
            row = dict(product=p["name"], source_id=int(ds.gaia["source_id"].values[j]),
                       TIC=(int(ds.gaia["ID"].values[j]) if np.isfinite(ds.gaia["ID"].values[j]) else None),
                       Tmag=float(T[j]), Tmag_src=ds.gaia["Tmag_src"].values[j], G=float(ds.gaia["phot_g_mean_mag"].values[j]),
                       dE=float(e_[0]), dN=float(n_[0]), sep=r, frac_in_ap=float(frac[j]), flux_rel=float(flux[j]),
                       Fj_over_Fref=float(Fap[j] / Fref) if Fref > 0 else np.nan)
            for sig in "AB":
                for k, d in enumerate(DEPTHS[sig]):
                    row[f"dreq_{sig}{'_lo' if k else ''}"] = d * Fref / Fap[j] if Fap[j] > 0 else np.inf
            rows.append(row)
    return pd.DataFrame(rows), summary


if __name__ == "__main__":
    df, summary = census()
    df.to_csv(os.path.join(OUT, "neighbour_census_per_product.csv"), index=False)
    # per star: most favourable product (smallest required depth)
    agg = df.groupby("source_id").agg(TIC=("TIC", "first"), Tmag=("Tmag", "first"), Tmag_src=("Tmag_src", "first"),
                                      G=("G", "first"), dE=("dE", "first"), dN=("dN", "first"), sep=("sep", "first"),
                                      dreq_A_min=("dreq_A", "min"), dreq_A_lo_min=("dreq_A_lo", "min"),
                                      dreq_B_min=("dreq_B", "min"), dreq_B_lo_min=("dreq_B_lo", "min"),
                                      best_product_A=("dreq_A", lambda s: df.loc[s.idxmin(), "product"]))
    # per-product required depth for the three manuscript products
    for prod in ["S6_tspoc_tp", "S32_tpf", "QLP98_3x3_best"]:
        sub = df[df["product"] == prod].set_index("source_id")
        agg[f"dreq_A_{prod}"] = sub["dreq_A"]
        agg[f"dreq_B_{prod}"] = sub["dreq_B"]
    agg = agg.sort_values("sep")
    agg.to_csv(os.path.join(OUT, "neighbour_census.csv"))
    cap = {}
    for sig in "AB":
        for thr in [1.0, 0.5]:
            for col in [f"dreq_{sig}_min", f"dreq_{sig}_lo_min"]:
                c = agg[agg[col] <= thr]
                cap[f"{col}<={thr}"] = dict(n=int(len(c)), nearest_sep=float(c["sep"].min()) if len(c) else None,
                                            stars=c[["TIC", "Tmag", "sep", "dE", "dN", col]].reset_index().to_dict("records"))
    # the manuscript-relevant rule: capable in the products the manuscript used (all three must allow it?)
    dump(dict(products=summary, capable=cap), "neighbour_census_summary.json")
    print(json.dumps(jsonable(summary), indent=1))
    for k, v in cap.items():
        print(k, v["n"], v["nearest_sep"])
    pd.set_option("display.width", 250)
    print(agg[agg["dreq_A_lo_min"] <= 2][["TIC", "Tmag", "sep", "dE", "dN", "dreq_A_min", "dreq_A_lo_min", "best_product_A",
                                          "dreq_A_S6_tspoc_tp", "dreq_A_S32_tpf", "dreq_A_QLP98_3x3_best", "dreq_B_min"]].to_string())
