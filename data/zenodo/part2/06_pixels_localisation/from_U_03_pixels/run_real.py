"""Real-data difference images and PRF localisation for all sectors/products/variants (audit 03_pixels)."""
import sys, time
from pixloc import *

CONFIGS = [
    # name, key, variant, kwargs
    ("S32_tpf", "S32_tpf", "default", {}),
    ("S32_tpf_clean", "S32_tpf", "clean", {}),
    ("S32_ffi", "S32_ffi", "default", {}),
    ("S98_clean", "S98_ffi", "clean", {}),
    ("S98_clean29", "S98_ffi", "clean29", {}),
    ("S98_clean29_W3", "S98_ffi", "clean29", dict(W_hours=3.0)),
    ("S98_clean29_W6", "S98_ffi", "clean29", dict(W_hours=6.0)),
    ("S98_clean29_deg2", "S98_ffi", "clean29", dict(deg=2)),
    ("S98_clean29_fit7", "S98_ffi", "clean29", dict(fit_half=3)),
    ("S98_all", "S98_ffi", "default", {}),
    ("S98_clean_q30", "S98_ffi", "clean_q30", {}),
    ("S6_ffi", "S6_ffi", "default", {}),
    ("S6_step", "S6_ffi", "default", dict(dump_step=True)),
    ("S6_dumpx0.5h", "S6_ffi", "default", dict(dump_exclude_h=0.5)),
    ("S6_dumpx1h", "S6_ffi", "default", dict(dump_exclude_h=1.0)),
    # baseline-model sensitivity
    ("S32_tpf_W3", "S32_tpf", "default", dict(W_hours=3.0)),
    ("S32_tpf_W6", "S32_tpf", "default", dict(W_hours=6.0)),
    ("S32_tpf_deg2", "S32_tpf", "default", dict(deg=2)),
    ("S98_clean_W3", "S98_ffi", "clean", dict(W_hours=3.0)),
    ("S98_clean_W6", "S98_ffi", "clean", dict(W_hours=6.0)),
    ("S98_clean_deg2", "S98_ffi", "clean", dict(deg=2)),
    ("S98_clean_fit7", "S98_ffi", "clean", dict(fit_half=3)),
    ("S32_ffi_fit7", "S32_ffi", "default", dict(fit_half=3)),
    # TESS-SPOC FFI-based TPs (same photons as TESScut; carry SPOC POS_CORR)
    ("S6_tspoc", "S6_tspoc_tp", "default", {}),
    ("S6_tspoc_dumpx0.5h", "S6_tspoc_tp", "default", dict(dump_exclude_h=0.5)),
    ("S32_tspoc", "S32_tspoc_tp", "default", {}),
]
MIN_SIDE = 0.25
NOMINAL_BOXDEPTH = dict(A=400e-6, B=900e-6)


def analyse(ds, sig, sf, events=None, F=None, want_nb=True):
    if events is None:
        events = [e for e in ds.event_windows(sig, min_side_frac=MIN_SIDE) if e["ok"]]
    Dm, Vm, Ds, Vs, chi2r = ds.combine(events, sig, F)
    infl = max(1.0, float(np.median(ds.fit_region(chi2r)))) if len(events) > 1 else 1.0
    Ve = Vm * infl
    r = ds.fit_position(Dm, Ve)
    sig_reg = np.sqrt(ds.fit_region(Ve))
    off, J = ds.sky_offset(r["x"], r["y"], sf["x_t"], sf["y_t"])
    c2t, at, bt, sat = ds.chi2_at(Dm, sig_reg, sf["x_t"], sf["y_t"])
    out = dict(n_events=len(events), event_tc=[e["tc"] for e in events], err_inflation=infl, fit=r,
               off_E=float(off[0]), off_N=float(off[1]), J=J,
               chi2_target=c2t, dchi2_target=c2t - r["chi2"], amp_target=at, sig_amp_target=sat,
               depth_target_ppm=1e6 * at / sf["target_flux"], depth_free_ppm=1e6 * r["a"] / sf["target_flux"],
               snr_target=at / sat)
    if want_nb:
        nb = []
        dxs, dys = sf["dx"], sf["dy"]
        for j in range(len(ds.gaia)):
            if j == ds.it:
                continue
            xj, yj = ds.g_x[j] + dxs, ds.g_y[j] + dys
            e_, n_ = offsets_arcsec(ds.g_ra[ds.it], ds.g_dec[ds.it], ds.g_ra[j], ds.g_dec[j])
            if np.hypot(e_[0], n_[0]) > 180:
                continue
            c2j = ds.chi2_at(Dm, sig_reg, xj, yj)
            nb.append(dict(source_id=int(ds.gaia["source_id"].values[j]), dE=float(e_[0]), dN=float(n_[0]),
                           Tmag=float(ds.gaia["Tmag_used"].values[j]), dchi2=c2j[0] - r["chi2"],
                           amp=c2j[1], depth_if_source=c2j[1] / (sf["target_flux"] * 10 ** (-0.4 * (
                               ds.gaia["Tmag_used"].values[j] - ds.gaia["Tmag_used"].values[ds.it])))))
        out["neighbours"] = nb
    return out, Dm, Ve, Ds


if __name__ == "__main__":
    only = sys.argv[1:] if len(sys.argv) > 1 else None
    res, imgs = {}, {}
    for name, key, var, kw in CONFIGS:
        if only and name not in only:
            continue
        t0 = time.time()
        ds = Dataset(key, var, **kw)
        sf = ds.scene_fit()
        res[name] = dict(key=key, variant=var, kwargs=kw, n_cadences=int(len(ds.t)), epoch=ds.epoch,
                         target_wcs_pix=[ds.xt_wcs, ds.yt_wcs], scene=sf, dump_times=ds.dump_times,
                         fit_region=[ds.ys.start, ds.ys.stop, ds.xs.start, ds.xs.stop])
        for sig in "AB":
            allev = ds.event_windows(sig, min_side_frac=MIN_SIDE)
            res[name][sig + "_all_events"] = [dict(tc=e["tc"], ok=e["ok"], in_sum=e["in_sum"],
                                                  n_in_exp=e["n_in_exp"], pre=e["pre"], post=e["post"])
                                             for e in allev]
            if sum(e["ok"] for e in allev) == 0:
                continue
            out, Dm, Ve, Ds = analyse(ds, sig, sf)
            res[name][sig] = out
            imgs[f"{name}_{sig}_D"] = Dm
            imgs[f"{name}_{sig}_V"] = Ve
            f = out["fit"]
            print(f"{name:16s} {sig} nev={out['n_events']:2d} off=({out['off_E']:+6.2f},{out['off_N']:+6.2f})\" "
                  f"depth_t={out['depth_target_ppm']:6.0f}ppm S/N_t={out['snr_target']:5.1f} "
                  f"dchi2_t={out['dchi2_target']:6.2f} chi2={f['chi2']:.1f}/{f['dof']} infl={out['err_inflation']:.2f}",
                  flush=True)
        imgs[f"{name}_med"] = ds.medimg
        print(f"   {name} done in {time.time()-t0:.1f}s; scene shift ({sf['dx']:+.3f},{sf['dy']:+.3f}) px", flush=True)
    tag = "" if not only else "_" + "_".join(only)
    dump(res, f"real_results{tag}.json")
    np.savez_compressed(os.path.join(OUT, f"real_images{tag}.npz"), **imgs)
