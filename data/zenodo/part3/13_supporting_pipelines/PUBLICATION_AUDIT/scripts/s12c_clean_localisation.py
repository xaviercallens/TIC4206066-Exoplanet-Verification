"""Step 12c: localisation robustness.
(1) S98 PRF localisation restricted to CLEAN intervals (QLP bits 2^29/2^30 and SPOC straylight all clear) for the
    real transits AND for the injections (so both sample the same kind of data).
(2) S32 SPOC 2-min: MOM_CENTR1/2 and POS_CORR1/2 in-transit shifts with an EMPIRICAL null (random epoch sets) instead of
    the white-noise error used before; residual centroid = MOM_CENTR - POS_CORR.
(3) Sky-frame combination of the per-sector PRF offsets (S6 FFI, S32 2-min TPF, S98 FFI clean) using the injection
    means and covariances -> combined offset and chi2 (2 dof) against 'on target'.
"""
import glob
import json
import os
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.dirname(__file__))
import lcio  # noqa
import prf_tools as pt  # noqa
import s12_pixels as s12  # noqa
from audit_common import AUD, save_json  # noqa

rng = np.random.default_rng(2026)


def qlp_clean_mask(t):
    q = fits.open(os.path.join(lcio.MAST, lcio.PRODUCTS["QLP_S98_SYS"][0]))[1].data
    tq = np.asarray(q["TIME"], float)
    qq = np.asarray(q["QUALITY"]).astype(np.int64)
    bad = ((qq & lcio.FROZEN_QLP_MASK) != 0) | ~np.isfinite(tq)
    idx = np.clip(np.searchsorted(tq, t), 1, len(tq) - 1)
    near = np.where(np.abs(tq[idx] - t) < np.abs(tq[idx - 1] - t), idx, idx - 1)
    return ~bad[near] & (np.abs(tq[near] - t) < 2.0 / 1440)


def localise(name, sig, clean=False, nsets=300, seed=5):
    C = s12.load_cube(name)
    if clean:
        keep = qlp_clean_mask(C["t"])
        for k in ("t", "cube", "bkg", "pc1", "pc2"):
            C[k] = C[k][keep]
    shape = C["cube"].shape[1:]
    sc0 = pt.scene(C["sources"], C["wcs"], C["tmid"], s12.TARGET_ID)
    tg0 = [s for s in sc0 if s["is_target"]][0]
    prf = pt.get_prf(C["cam"], C["ccd"], C["sector"], C["col0"] + tg0["x"], C["row0"] + tg0["y"])
    astro, srcs, _ = s12.astrometric_solution(C, prf)
    tg = [s for s in srcs if s["is_target"]][0]
    e = s12.EPH[sig]
    other = "B" if sig == "A" else "A"
    okc = ~lcio.in_transit(C["t"], other, 1.0, s12.EPH)
    eps = s12.epochs_for(C, sig)
    tcs = [x[1] for x in eps]
    D = s12.diff_image(C, C["cube"], tcs, e["T14"], okc).mean(0)
    allok = np.ones(len(C["t"]), bool)
    fsets = s12.random_epoch_sets(C, sig, len(tcs), 300, np.random.default_rng(seed))
    fakes = np.array([s12.diff_image(C, C["cube"], s, e["T14"], allok).mean(0) for s in fsets])
    sig_map = np.std(fakes, 0)
    G = s12.GridFitter(prf, shape, tg["x"], tg["y"])
    fr = G.fit(D, sig_map)
    ft = G.fit_fixed(D, sig_map, tg["x"], tg["y"])
    depth = ft["A"] / astro["target_flux"]
    isets = s12.random_epoch_sets(C, sig, len(tcs), nsets, np.random.default_rng(seed + 1))
    offs = []
    for s in isets:
        ci = s12.inject_on_target(C["cube"] + C["bkg"][:, None, None], C["bkg"], C["t"], s, e["T14"], depth, tg["x"], tg["y"], 30.0) \
            - C["bkg"][:, None, None]
        Di = s12.diff_image(C, ci, s, e["T14"], allok).mean(0)
        f = G.fit(Di, sig_map)
        offs.append((f["x"] - tg["x"], f["y"] - tg["y"]))
    o = np.array(offs)
    real = np.array([fr["x"] - tg["x"], fr["y"] - tg["y"]])
    mu, cov = o.mean(0), np.cov(o.T)
    # pixel -> sky Jacobian at the target (arcsec per px)
    w = C["wcs"]
    p0 = w.all_pix2world([[tg["x"], tg["y"]]], 0)[0]
    px = w.all_pix2world([[tg["x"] + 1, tg["y"]]], 0)[0]
    py = w.all_pix2world([[tg["x"], tg["y"] + 1]], 0)[0]
    cosd = np.cos(np.radians(p0[1]))
    J = np.array([[(px[0] - p0[0]) * cosd * 3600, (py[0] - p0[0]) * cosd * 3600],
                  [(px[1] - p0[1]) * 3600, (py[1] - p0[1]) * 3600]])      # [dRA, dDec] = J @ [dx, dy]
    return dict(name=name, signal=sig, clean=clean, n_epochs=len(tcs), epochs=[int(x[0]) for x in eps], depth_ppm=float(depth * 1e6),
                real_px=real.tolist(), inj_mean_px=mu.tolist(), inj_cov_px=cov.tolist(), n_inj=len(o),
                z=[float((real[i] - mu[i]) / np.sqrt(cov[i, i])) for i in (0, 1)],
                chi2_2dof=float((real - mu) @ np.linalg.inv(cov) @ (real - mu)),
                frac_inj_further=float(np.mean(np.hypot(o[:, 0] - mu[0], o[:, 1] - mu[1]) >= np.hypot(*(real - mu)))),
                J_arcsec_per_px=J.tolist(), real_sky_arcsec=(J @ (real - mu)).tolist(), sky_cov_arcsec2=(J @ cov @ J.T).tolist())


def spoc2min_centroids(sig, nnull=2000):
    d = fits.open(os.path.join(lcio.MAST, lcio.PRODUCTS["SPOC2min_S32_PDC"][0]))[1].data
    t = np.asarray(d["TIME"], float)
    q = np.asarray(d["QUALITY"]).astype(np.int64)
    cols = {c: np.asarray(d[c], float) for c in ("MOM_CENTR1", "MOM_CENTR2", "POS_CORR1", "POS_CORR2")}
    ok = np.isfinite(t) & ((q & lcio.DEFAULT_BITMASK) == 0) & np.all([np.isfinite(v) for v in cols.values()], axis=0)
    t = t[ok]
    cols = {k: v[ok] for k, v in cols.items()}
    cols["RES1"] = cols["MOM_CENTR1"] - cols["POS_CORR1"]
    cols["RES2"] = cols["MOM_CENTR2"] - cols["POS_CORR2"]
    e = s12.EPH[sig]
    other = "B" if sig == "A" else "A"
    okc = ~lcio.in_transit(t, other, 1.0, s12.EPH)
    T14 = e["T14"]

    def stat(tcs):
        out = {k: [] for k in cols}
        for tc in tcs:
            i = okc & (np.abs(t - tc) < 0.4 * T14)
            o = okc & (np.abs(t - tc) > 0.75 * T14) & (np.abs(t - tc) < 2.5 * T14)
            if i.sum() < 5 or o.sum() < 10:
                continue
            X = np.vstack([np.ones(o.sum()), t[o] - tc]).T
            for k, v in cols.items():
                out[k].append(v[i].mean() - np.linalg.lstsq(X, v[o], rcond=None)[0][0])
        return {k: float(np.mean(v)) for k, v in out.items()}, len(out["RES1"])

    n0 = int(np.floor((t.min() - e["T0"]) / e["P"]))
    n1 = int(np.ceil((t.max() - e["T0"]) / e["P"]))
    tcs = [e["T0"] + n * e["P"] for n in range(n0, n1 + 1) if t.min() < e["T0"] + n * e["P"] < t.max()]
    real, nreal = stat(tcs)
    bad = lcio.in_transit(t, "A", 3.0, s12.EPH) | lcio.in_transit(t, "B", 3.0, s12.EPH)
    cand = t[~bad]
    nulls = {k: [] for k in cols}
    for _ in range(nnull):
        cs = rng.choice(cand, nreal, replace=False)
        v, n = stat(cs)
        if n == nreal:
            for k in cols:
                nulls[k].append(v[k])
    res = {}
    for k in cols:
        a = np.array(nulls[k])
        res[k] = dict(real_mpx=real[k] * 1e3, null_mean_mpx=float(a.mean() * 1e3), null_sd_mpx=float(a.std() * 1e3),
                      z=float((real[k] - a.mean()) / a.std()), p_two_sided=float(np.mean(np.abs(a - a.mean()) >= abs(real[k] - a.mean()))))
    return dict(signal=sig, n_transits=nreal, n_null=len(nulls["RES1"]), columns=res)


if __name__ == "__main__":
    out = {"clean_S98": {}, "spoc2min_centroids": {}, "combined": {}, "combined_without_S6": {}, "per_set": {}}
    per = {}
    for sig in ("A", "B"):
        r = localise("S98_ffi", sig, clean=True)
        out["clean_S98"][sig] = r
        print(f"S98 clean {sig}: n={r['n_epochs']} {r['epochs']} real {np.round(r['real_px'], 3)} inj {np.round(r['inj_mean_px'], 3)} "
              f"z {np.round(r['z'], 2)} chi2(2) {r['chi2_2dof']:.1f} frac_further {r['frac_inj_further']:.3f} sky {np.round(r['real_sky_arcsec'], 1)}", flush=True)
        c = spoc2min_centroids(sig)
        out["spoc2min_centroids"][sig] = c
        print(f"S32 2-min centroids {sig}: " + "; ".join(f"{k} {v['real_mpx']:+.3f} (null {v['null_mean_mpx']:+.3f}±{v['null_sd_mpx']:.3f}, z {v['z']:+.2f}, p {v['p_two_sided']:.3f})" for k, v in c["columns"].items()), flush=True)
        per[sig] = [r]
        for nm in ("S6_ffi", "S32_tpf2min"):
            rr = localise(nm, sig, clean=False)
            per[sig].append(rr)
            print(f"{nm} {sig}: real {np.round(rr['real_px'], 3)} z {np.round(rr['z'], 2)} chi2(2) {rr['chi2_2dof']:.1f} sky {np.round(rr['real_sky_arcsec'], 1)}", flush=True)
        out["per_set"][sig] = per[sig]
        # same combination without S6 (every S6 A transit lies 0.4-3.5 h from a momentum dump; see s12d)
        Si2, Sx2 = np.zeros((2, 2)), np.zeros(2)
        for rr in per[sig]:
            if rr["name"].startswith("S6") or rr["n_epochs"] == 0:
                continue
            Ci = np.linalg.inv(np.array(rr["sky_cov_arcsec2"]))
            Si2 += Ci
            Sx2 += Ci @ np.array(rr["real_sky_arcsec"])
        C2 = np.linalg.inv(Si2)
        m2 = C2 @ Sx2
        out["combined_without_S6"][sig] = dict(mean_offset_sky_arcsec=m2.tolist(), sigma_arcsec=[float(np.sqrt(C2[0, 0])), float(np.sqrt(C2[1, 1]))],
                                               chi2_mean_2dof=float(m2 @ Si2 @ m2))
        print(f"COMBINED without S6 {sig}: ({m2[0]:+.1f}, {m2[1]:+.1f})\" ± ({np.sqrt(C2[0,0]):.1f}, {np.sqrt(C2[1,1]):.1f}); chi2 {m2 @ Si2 @ m2:.1f}", flush=True)
        # combination in the sky frame (offsets relative to injection means)
        Sinv = np.zeros((2, 2))
        Sx = np.zeros(2)
        chi = 0.0
        for rr in per[sig]:
            if rr["n_epochs"] == 0:
                continue
            Cs = np.array(rr["sky_cov_arcsec2"])
            v = np.array(rr["real_sky_arcsec"])
            Ci = np.linalg.inv(Cs)
            Sinv += Ci
            Sx += Ci @ v
            chi += float(v @ Ci @ v)
        Cc = np.linalg.inv(Sinv)
        m = Cc @ Sx
        out["combined"][sig] = dict(sets=[rr["name"] + ("_clean" if rr["clean"] else "") for rr in per[sig]],
                                    mean_offset_sky_arcsec=m.tolist(), cov_arcsec2=Cc.tolist(),
                                    sigma_arcsec=[float(np.sqrt(Cc[0, 0])), float(np.sqrt(Cc[1, 1]))],
                                    chi2_mean_2dof=float(m @ Sinv @ m), chi2_total=chi, dof_total=2 * len(per[sig]))
        cb = out["combined"][sig]
        print(f"COMBINED {sig}: offset (dRA, dDec) = ({m[0]:+.1f}, {m[1]:+.1f})\" ± ({cb['sigma_arcsec'][0]:.1f}, {cb['sigma_arcsec'][1]:.1f}); "
              f"chi2 of mean (2 dof) {cb['chi2_mean_2dof']:.1f}; total chi2 {chi:.1f}/{cb['dof_total']}", flush=True)
        save_json(out, "out/s12c_clean_localisation.json")
