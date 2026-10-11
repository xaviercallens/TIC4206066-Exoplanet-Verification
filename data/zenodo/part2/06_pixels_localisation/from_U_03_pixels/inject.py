"""G3 injection-recovery calibration of the difference-image localisation (audit 03_pixels).

On-target injections are model-free: the real pixels near the target (median out-of-transit image, pixels within
R_INJ px of the target) are multiplied by (1 - depth * box(t)) at random out-of-transit epochs, i.e.
F'(t) = F(t) - depth * box(t) * M.  The real events' cadence patterns are transplanted onto the fake epochs
(same relative cadences, same other-signal masks, same acceptance rules, same number of events).
Off-target injections use the same real target image shifted by an integer number of pixels (no PRF model),
so neither kind of injection uses the PRF model that is used for fitting.

usage: python inject.py JOBNAME [JOBNAME ...]   (see JOBS)
"""
import sys, time
from run_real import *

R_INJ = 3.5
JOBS = {
    # name: (config name in run_real.CONFIGS, signal, ntrial, shift(px), seed)
    "S32_tpf_A": ("S32_tpf", "A", 300, (0, 0), 1),
    "S32_tpf_B": ("S32_tpf", "B", 300, (0, 0), 2),
    "S32_ffi_A": ("S32_ffi", "A", 300, (0, 0), 3),
    "S32_ffi_B": ("S32_ffi", "B", 300, (0, 0), 4),
    "S98_clean_A": ("S98_clean", "A", 300, (0, 0), 5),
    "S98_clean_B": ("S98_clean", "B", 300, (0, 0), 6),
    "S98_all_A": ("S98_all", "A", 300, (0, 0), 7),
    "S98_all_B": ("S98_all", "B", 300, (0, 0), 8),
    "S6_ffi_A": ("S6_ffi", "A", 300, (0, 0), 9),
    "S6_ffi_B": ("S6_ffi", "B", 300, (0, 0), 10),
    "S6_step_A": ("S6_step", "A", 300, (0, 0), 11),
    "S6_dumpx0.5h_A": ("S6_dumpx0.5h", "A", 300, (0, 0), 12),
    "S32_tpf_deg2_A": ("S32_tpf_deg2", "A", 300, (0, 0), 13),
    "S98_clean_deg2_A": ("S98_clean_deg2", "A", 300, (0, 0), 14),
    "S98_clean29_A": ("S98_clean29", "A", 300, (0, 0), 15),
    "S98_clean29_B": ("S98_clean29", "B", 300, (0, 0), 16),
    "S98_clean29_deg2_A": ("S98_clean29_deg2", "A", 300, (0, 0), 17),
    "S98_clean29_W3_A": ("S98_clean29_W3", "A", 300, (0, 0), 18),
    "S98_clean29_W6_A": ("S98_clean29_W6", "A", 300, (0, 0), 19),
    "S32_tpf_W3_A": ("S32_tpf_W3", "A", 300, (0, 0), 20),
    "S32_tpf_W6_A": ("S32_tpf_W6", "A", 300, (0, 0), 26),
    "S6_dumpx1h_A": ("S6_dumpx1h", "A", 300, (0, 0), 27),
    "S98_clean29_A_shift+2+0": ("S98_clean29", "A", 80, (2, 0), 28),
    # off-target (shifted real image) checks
    "S32_tpf_A_shift+2+0": ("S32_tpf", "A", 80, (2, 0), 21),
    "S32_tpf_A_shift+0+2": ("S32_tpf", "A", 80, (0, 2), 22),
    "S98_clean_A_shift+2+0": ("S98_clean", "A", 80, (2, 0), 23),
    "S98_clean_A_shift+0-2": ("S98_clean", "A", 80, (0, -2), 24),
    "S98_clean_A_shift+1+0": ("S98_clean", "A", 80, (1, 0), 25),
}


def make_dataset(cfgname):
    for name, key, var, kw in CONFIGS:
        if name == cfgname:
            return Dataset(key, var, **kw)
    raise KeyError(cfgname)


def run_job(job, out_path):
    cfg, sig, ntrial, shift, seed = JOBS[job]
    rng = np.random.default_rng(seed)
    ds = make_dataset(cfg)
    sf = ds.scene_fit()
    e = EPH[sig]
    T14 = e["T14"]
    real_events = [ev for ev in ds.event_windows(sig, min_side_frac=MIN_SIDE) if ev["ok"]]
    real, _, _, _ = analyse(ds, sig, sf, events=real_events, want_nb=False)
    # injected depth: measured on-target box depth of the real signal if positive, else nominal
    depth = real["depth_target_ppm"] * 1e-6
    depth_src = "measured"
    if not (depth > 50e-6):
        depth, depth_src = NOMINAL_BOXDEPTH[sig], "nominal"
    # model-free template image
    yy, xx = np.indices(ds.medimg.shape)
    M = ds.medimg * (np.hypot(xx - sf["x_t"], yy - sf["y_t"]) <= R_INJ)
    sx, sy = shift
    if sx or sy:
        M = np.roll(M, (sy, sx), axis=(0, 1))
        if sy > 0: M[:sy] = 0
        if sy < 0: M[sy:] = 0
        if sx > 0: M[:, :sx] = 0
        if sx < 0: M[:, sx:] = 0
    true_off, _ = ds.sky_offset(sf["x_t"] + sx, sf["y_t"] + sy, sf["x_t"], sf["y_t"])
    # forbidden cadences: real transits of both signals (+/- 30 min)
    same_real = in_any_transit(ds.t, ds.dt, e, PAD_OTHER)
    om = ds.other_mask(sig)
    t = ds.t
    tmin, tmax = t.min(), t.max()
    patterns = [np.sort(t[ev["idx"]] - ev["tc"]) for ev in real_events]
    trials = []
    t0 = time.time()
    n_fail = 0
    for k in range(ntrial):
        fake = []
        used_tc = []
        for pat, ev in zip(patterns, real_events):
            for attempt in range(3000):
                tc = rng.uniform(tmin, tmax)
                if any(abs(tc - u) < 2 * (T14 / 2 + ds.W) for u in used_tc):
                    continue
                sel = (np.abs(t - tc) <= T14 / 2 + ds.W) & ~om & ~same_real
                idx = np.where(sel)[0]
                if len(idx) < 0.9 * len(ev["idx"]):
                    continue
                rel = t[idx] - tc
                j = np.clip(np.searchsorted(pat, rel), 1, len(pat) - 1)
                near = np.minimum(np.abs(pat[j] - rel), np.abs(pat[j - 1] - rel)) < ds.dt / 2
                idx = idx[near]
                if len(idx) < 0.9 * len(ev["idx"]):
                    continue
                cand = ds.event_windows(sig, centres=[tc], min_side_frac=MIN_SIDE)
                if not cand or not cand[0]["ok"]:
                    continue
                # the event must not contain the real in-transit cadences of this signal
                if same_real[idx].any():
                    continue
                frac = box_fraction(t[idx], ds.dt, tc, T14)
                if frac.sum() < 0.6 * T14 / ds.dt:
                    continue
                Fsel = ds.F[idx] - depth * frac[:, None, None] * M[None]
                fake.append(dict(tc=tc, idx=idx, Fsel=Fsel))
                used_tc.append(tc)
                break
            else:
                n_fail += 1
        if len(fake) < len(real_events):
            continue
        out, Dm, Ve, _ = analyse(ds, sig, sf, events=fake, want_nb=False)
        rec = dict(off_E=out["off_E"], off_N=out["off_N"], dchi2_target=out["dchi2_target"],
                   snr_target=out["snr_target"], depth_target_ppm=out["depth_target_ppm"],
                   chi2=out["fit"]["chi2"], infl=out["err_inflation"])
        if sx or sy:
            sig_reg = np.sqrt(ds.fit_region(Ve))
            c2s = ds.chi2_at(Dm, sig_reg, sf["x_t"] + sx, sf["y_t"] + sy)[0]
            rec["dchi2_truepos"] = c2s - out["fit"]["chi2"]
        trials.append(rec)
        if (k + 1) % 25 == 0:
            print(f"{job}: {k+1}/{ntrial} trials, {time.time()-t0:.0f}s", flush=True)
            dump(dict(job=job, partial=True, trials=trials), out_path)
    res = dict(job=job, config=cfg, signal=sig, shift_px=shift, true_offset_arcsec=true_off, depth_injected=depth,
               depth_source=depth_src, R_inj_px=R_INJ, n_events=len(real_events), seed=seed,
               real=dict(off_E=real["off_E"], off_N=real["off_N"], dchi2_target=real["dchi2_target"],
                         snr_target=real["snr_target"], depth_target_ppm=real["depth_target_ppm"]),
               n_trials=len(trials), n_event_placement_failures=n_fail, trials=trials, partial=False,
               runtime_s=time.time() - t0)
    dump(res, out_path)
    return res


if __name__ == "__main__":
    os.makedirs(os.path.join(OUT, "inj"), exist_ok=True)
    for job in sys.argv[1:]:
        p = os.path.join("inj", f"inj_{job}.json")
        if os.path.exists(os.path.join(OUT, p)):
            try:
                if not json.load(open(os.path.join(OUT, p))).get("partial", True):
                    print("skip", job); continue
            except Exception:
                pass
        run_job(job, p)
        print("done", job, flush=True)
