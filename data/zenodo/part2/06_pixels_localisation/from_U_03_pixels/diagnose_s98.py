"""Diagnostics of the S98/S32 localisation of signal A: jackknife over events, and a variant that mimics the
baseline choice described in the package's method notes (baseline 0.75-2.5 T14 from mid-transit, one-sided events
allowed with a constant baseline), implemented here from that verbal description only (no code reused)."""
from run_real import *


class PriorLike(Dataset):
    def event_windows(self, sig, centres=None, min_in_frac=0.6, min_side_frac=0.3):
        e = EPH[sig]; T14 = e["T14"]
        if centres is None:
            _, centres = transit_centres(e, self.t.min() - 1, self.t.max() + 1)
        om = self.other_mask(sig)
        out = []
        for tc in centres:
            dtc = np.abs(self.t - tc)
            intr = dtc <= T14 / 2
            base = (dtc >= 0.75 * T14) & (dtc <= 2.5 * T14)
            idx = np.where((intr | base) & ~om)[0]
            if len(idx) < 5:
                continue
            frac = box_fraction(self.t[idx], self.dt, tc, T14)
            nb = np.sum(base[idx])
            ok = frac.sum() >= 0.6 * T14 / self.dt and nb >= 4
            out.append(dict(tc=tc, idx=idx, ok=bool(ok), in_sum=float(frac.sum()), n_in_exp=T14 / self.dt,
                            pre=int(np.sum(self.t[idx] < tc - T14 / 2)), post=int(np.sum(self.t[idx] > tc + T14 / 2))))
        return out

    def event_diff(self, idx, tc, sig, F=None, Fsel=None):
        if F is None:
            F = self.F
        t = self.t[idx]
        T14 = EPH[sig]["T14"]
        f = box_fraction(t, self.dt, tc, T14)
        pre = np.sum(t < tc - T14 / 2); post = np.sum(t > tc + T14 / 2)
        cols = [np.ones_like(t)]
        if pre >= 2 and post >= 2:
            cols.append(t - tc)
        cols.append(f)
        X = np.vstack(cols).T
        Y = (F[idx] if Fsel is None else Fsel).reshape(len(idx), -1)
        XtX_inv = np.linalg.pinv(X.T @ X)
        beta = XtX_inv @ (X.T @ Y)
        resid = Y - X @ beta
        s2 = (resid ** 2).sum(0) / max(len(idx) - X.shape[1], 1)
        return (-beta[-1]).reshape(F.shape[1:]), (s2 * XtX_inv[-1, -1]).reshape(F.shape[1:])


res = {}
for name, cls, key, var in [("S98_clean29", Dataset, "S98_ffi", "clean29"), ("S98_clean29_priorlike", PriorLike, "S98_ffi", "clean29"),
                            ("S32_tpf", Dataset, "S32_tpf", "default"), ("S32_tpf_priorlike", PriorLike, "S32_tpf", "default"),
                            ("S6_ffi_priorlike", PriorLike, "S6_ffi", "default")]:
    ds = cls(key, var)
    sf = ds.scene_fit()
    for sig in "AB":
        evs = [e for e in ds.event_windows(sig, min_side_frac=MIN_SIDE) if e["ok"]]
        out, Dm, Ve, Ds = analyse(ds, sig, sf, events=evs, want_nb=False)
        jk = []
        if len(evs) > 2:
            for k in range(len(evs)):
                o2, *_ = analyse(ds, sig, sf, events=[e for i, e in enumerate(evs) if i != k], want_nb=False)
                jk.append(dict(dropped_tc=evs[k]["tc"], off_E=o2["off_E"], off_N=o2["off_N"]))
            jE = np.array([j["off_E"] for j in jk]); jN = np.array([j["off_N"] for j in jk])
            n = len(jk)
            jk_sig = [float(np.sqrt((n - 1) / n * ((jE - jE.mean()) ** 2).sum())),
                      float(np.sqrt((n - 1) / n * ((jN - jN.mean()) ** 2).sum()))]
        else:
            jk_sig = None
        res[f"{name}_{sig}"] = dict(n_events=len(evs), event_tc=[e["tc"] for e in evs], off=[out["off_E"], out["off_N"]],
                                   snr_target=out["snr_target"], dchi2_target=out["dchi2_target"],
                                   jackknife=jk, jackknife_sigma=jk_sig)
        print(f"{name:24s} {sig} nev={len(evs):2d} off=({out['off_E']:+6.1f},{out['off_N']:+6.1f}) S/N={out['snr_target']:.1f} "
              f"jackknife sigma={None if jk_sig is None else np.round(jk_sig, 1)} "
              f"range E[{min([j['off_E'] for j in jk], default=np.nan):+.1f},{max([j['off_E'] for j in jk], default=np.nan):+.1f}] "
              f"N[{min([j['off_N'] for j in jk], default=np.nan):+.1f},{max([j['off_N'] for j in jk], default=np.nan):+.1f}]",
              flush=True)
dump(res, "diagnose_s98.json")
