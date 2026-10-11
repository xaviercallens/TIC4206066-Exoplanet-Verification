"""Injection-recovery calibration of per-sector and per-event depth errors for signal B.

Usage: python s03_injections.py <cfgname> <mode: sector|event> <ntrials> [seed]
  sector : the whole B transit train is injected with a random phase offset (0.08-0.92 P from the real events), i.e. the
           real event spacing is kept and the number of covered events follows the real gaps; the sector depth is
           re-measured with exactly the same pipeline (re-detrending included).
  event  : a single B transit is injected at a random cadence >= 0.6 d from any real B event.
The injected depth is 1000 ppm (template of cfg). For each trial the normalised residual z = (D_rec - D_inj)/eD_white is
stored; kappa = std(z) is the factor by which white errors must be inflated (red noise + detrending), per sector.
Output: out/s03_inj_<cfgname>_<mode>.json
"""
import sys, time
from b_run import *
from s_configs import CONFIGS, template_for

cfgname = sys.argv[1]; mode = sys.argv[2]; ntr = int(sys.argv[3]); seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
cfg = CONFIGS[cfgname]; tpl = template_for(cfg)
rng = np.random.default_rng(seed)
D_INJ = 1000.0
res = dict(cfg=cfgname, mode=mode, ntrials=ntr, D_inj=D_INJ, template=tpl.shape, sectors={})
t0 = time.time()
for sec in (6, 32, 98):
    lc = get_lc(cfg["products"][sec], cfg.get("s98_mask", MASK_S98_BASE), cfg.get("spoc_mask", 16575))
    real_tc = EPH_B["T0"] + EPH_B["P"] * np.arange(-400, 20)
    rows = []
    tries = 0
    while len(rows) < ntr and tries < ntr * 5:
        tries += 1
        if mode == "sector":
            dlt = rng.uniform(0.08, 0.92) * EPH_B["P"]
            eph = dict(P=EPH_B["P"], T0=EPH_B["T0"] + dlt, T14=EPH_B["T14"])
        else:
            tc = rng.choice(lc["t"])
            if np.min(np.abs(real_tc - tc)) < 0.6:
                continue
            eph = dict(P=1.0e4, T0=tc, T14=EPH_B["T14"])
        r = inject_measure(lc, tpl, cfg, eph, D_INJ)
        if r is None:
            continue
        rows.append(dict(T0=eph["T0"], D=r["D"], eD=r["eD_white"], chi2_red=r["chi2_red"], n_ev=len(r["epochs"]), cov=r["cov"],
                         De=r["D_event"], eDe=r["eD_event_white"], A_in=r["A_in_window"]))
    D = np.array([x["D"] for x in rows]); eD = np.array([x["eD"] for x in rows]); z = (D - D_INJ) / eD
    ze = np.concatenate([(np.array(x["De"]) - D_INJ) / np.array(x["eDe"]) for x in rows])
    cove = np.concatenate([np.array(x["cov"]) for x in rows])
    full = cove >= 0.9
    def rstd(x):
        return float(1.4826 * np.median(np.abs(x - np.median(x))))
    # bootstrap over trials for kappa uncertainty (trials are not fully independent; see REPORT)
    bs = [np.std(rng.choice(z, len(z))) for _ in range(500)]
    res["sectors"][str(sec)] = dict(
        product=cfg["products"][sec], n=len(rows),
        mean_bias_ppm=float(np.mean(D - D_INJ)), median_bias_ppm=float(np.median(D - D_INJ)),
        sd_D_ppm=float(np.std(D)), median_eD_white=float(np.median(eD)),
        kappa=float(np.std(z)), kappa_robust=rstd(z), kappa_bs_sd=float(np.std(bs)),
        kappa_event=float(np.std(ze)), kappa_event_robust=rstd(ze), kappa_event_full=float(np.std(ze[full])) if full.sum() > 5 else None,
        n_events_hist={str(k): int(v) for k, v in zip(*np.unique([x["n_ev"] for x in rows], return_counts=True))},
        sd_event_D_full_ppm=float(np.std(np.concatenate([np.array(x["De"])[np.array(x["cov"]) >= 0.9] for x in rows]))) if full.sum() > 5 else None,
        rows=rows)
    s = res["sectors"][str(sec)]
    print(f"S{sec} {cfg['products'][sec]}: n={s['n']} bias={s['mean_bias_ppm']:.0f} sdD={s['sd_D_ppm']:.0f} med eD={s['median_eD_white']:.0f} "
          f"kappa={s['kappa']:.2f} (robust {s['kappa_robust']:.2f}, bs {s['kappa_bs_sd']:.2f}) kappa_ev={s['kappa_event']:.2f} "
          f"sd_ev_full={s['sd_event_D_full_ppm']} nev={s['n_events_hist']} t={time.time()-t0:.0f}s", flush=True)
jdump(res, f"s03_inj_{cfgname}_{mode}.json")
print("done", time.time() - t0)
