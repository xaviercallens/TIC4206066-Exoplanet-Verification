"""Step 2: S98-only BLS and TLS search (A removed), rank/SDE of B, aliases, odd/even, secondary."""
import numpy as np, json, sys, time
import common_b as c
import s98data as sd


def main_bls():

    out = {}
    for var in ["QLP_DET_bw075", "QLP_DET_raw", "QLP_SYSRM_bw075", "OWN_bhalves_bw075", "OWN_bkgcut_bw075"]:
        t, f = sd.get(var)
        T = t.max() - t.min()
        P = sd.period_grid(T)
        r, pw, snr = sd.bls_spectrum(t, f, P)
        S = sd.sde(pw)
        # distinct peaks: local maxima sorted, separated by > 0.5% in period and not harmonics-merged
        order = np.argsort(pw)[::-1]
        peaks = []
        for i in order:
            if all(abs(P[i] - P[j]) / P[j] > 0.005 for j in peaks):
                peaks.append(i)
            if len(peaks) >= 15: break
        near = np.abs(P - 11.133) < 0.03
        iB = np.where(near)[0][np.argmax(pw[near])]
        rankB = int(np.sum(pw[peaks] > pw[iB]) + 1) if iB not in peaks else peaks.index(iB) + 1
        def alias(Pa, tol=0.003):
            m = np.abs(P - Pa) / Pa < tol
            if not m.any(): return None
            i = np.where(m)[0][np.argmax(pw[m])]
            return dict(P=float(P[i]), SDE=float(S[i]), depth_ppm=float(r.depth[i] * 1e6), snr=float(snr[i]))
        res = dict(n_points=int(len(t)), baseline_d=float(T), n_periods=int(len(P)),
                   B=dict(P=float(P[iB]), SDE=float(S[iB]), power=float(pw[iB]), depth_ppm=float(r.depth[iB] * 1e6), snr=float(snr[iB]),
                          dur_h=float(r.duration[iB] * 24), T0=float(r.transit_time[iB]), rank=rankB),
                   global_max=dict(P=float(P[np.argmax(pw)]), SDE=float(S.max()), snr=float(snr[np.argmax(pw)])),
                   max_snr_any_period=float(np.nanmax(snr)), P_at_max_snr=float(P[np.nanargmax(snr)]),
                   top_peaks=[dict(P=float(P[i]), SDE=float(S[i]), depth_ppm=float(r.depth[i] * 1e6), snr=float(snr[i]), dur_h=float(r.duration[i] * 24)) for i in peaks[:10]],
                   aliases={"P/2": alias(P[iB] / 2), "P/3": alias(P[iB] / 3), "3P/2": alias(1.5 * P[iB]), "2P/3": alias(2 * P[iB] / 3)})
        # 2P (22.27 d) lies outside the 1-20 d search grid: evaluate a local grid and express it in units of the 1-20 d spectrum
        P2 = np.linspace(2 * P[iB] * 0.997, 2 * P[iB] * 1.003, 200)
        r2, pw2, snr2 = sd.bls_spectrum(t, f, P2)
        i2 = np.argmax(pw2)
        res["aliases"]["2P"] = dict(P=float(P2[i2]), SDE=float((pw2[i2] - np.mean(pw)) / np.std(pw)), depth_ppm=float(r2.depth[i2] * 1e6), snr=float(snr2[i2]),
                                    note="SDE computed with the mean/std of the 1-20 d spectrum")
        # odd/even + secondary at the BLS ephemeris of B
        PB, T0B = P[iB], r.transit_time[iB]
        ep = np.round((t - T0B) / PB).astype(int)
        ph = t - (T0B + ep * PB)
        def box(sel):
            inn = sel & (np.abs(ph) < c.T14_B / 2)
            oot = sel & (np.abs(ph) > c.T14_B / 2) & (np.abs(ph) < 3 * c.T14_B)
            if inn.sum() < 5: return None
            d = np.mean(f[oot]) - np.mean(f[inn]); e = np.std(f[oot]) * np.sqrt(1 / inn.sum() + 1 / oot.sum())
            return dict(depth_ppm=float(d * 1e6), err_ppm=float(e * 1e6), n_in=int(inn.sum()), epochs=sorted(set(ep[inn].tolist())))
        res["odd"] = box(ep % 2 == 1); res["even"] = box(ep % 2 == 0)
        if res["odd"] and res["even"]:
            res["odd_even_diff_sigma"] = float((res["odd"]["depth_ppm"] - res["even"]["depth_ppm"]) / np.hypot(res["odd"]["err_ppm"], res["even"]["err_ppm"]))
        # secondary: box of T14 at phases 0.1..0.9 (step T14/2); empirical sigma from the same scan
        phase = ((t - T0B) / PB) % 1.0
        sec = []
        for p0 in np.arange(0.1, 0.9001, (c.T14_B / 2) / PB):
            dphi = (phase - p0 + 0.5) % 1.0 - 0.5
            inn = np.abs(dphi * PB) < c.T14_B / 2
            oot = (np.abs(dphi * PB) > c.T14_B / 2) & (np.abs(dphi * PB) < 3 * c.T14_B)
            if inn.sum() < 5: continue
            sec.append((p0, (np.mean(f[oot]) - np.mean(f[inn])) * 1e6, np.std(f[oot]) / np.sqrt(inn.sum()) * 1e6))
        sec = np.array(sec)
        emp = c.robust_std(sec[:, 1])
        i05 = np.argmin(np.abs(sec[:, 0] - 0.5))
        res["secondary"] = dict(phase05_depth_ppm=float(sec[i05, 1]), phase05_err_white_ppm=float(sec[i05, 2]), empirical_sd_ppm=float(emp),
                                max_depth_ppm=float(sec[:, 1].max()), phase_at_max=float(sec[np.argmax(sec[:, 1]), 0]),
                                max_over_emp=float(sec[:, 1].max() / emp), n_phases=int(len(sec)))
        out[var] = res
        print(var, json.dumps({k: res[k] for k in ("B", "global_max", "odd", "even", "secondary")}, default=float), flush=True)
        print("   top peaks:", [(round(p["P"], 4), round(p["SDE"], 2)) for p in res["top_peaks"]], flush=True)
        print("   aliases:", {k: (v and (round(v["P"], 4), round(v["SDE"], 2))) for k, v in res["aliases"].items()}, flush=True)
    json.dump(out, open("out/b02_s98_bls.json", "w"), indent=1)



def main_tls():
    # TLS on QLP_DET and OWN_bhalves
    from transitleastsquares import transitleastsquares
    tls_out = {}
    for var in ["QLP_DET_bw075", "OWN_bhalves_bw075"]:
        t, f = sd.get(var)
        t0 = time.time()
        import contextlib, io
        model = transitleastsquares(t, f, dy=np.full_like(f, c.robust_std(f)))
        with contextlib.redirect_stdout(io.StringIO()):
          rr = model.power(period_min=1.0, period_max=20.0, oversampling_factor=3, duration_grid_step=1.1, use_threads=1,
                         show_progress_bar=False, R_star=0.62, M_star=0.61, R_star_min=0.5, R_star_max=0.75, M_star_min=0.5, M_star_max=0.75)
        Ps = np.asarray(rr.periods); Sd = np.asarray(rr.power)
        near = np.abs(Ps - 11.133) < 0.03
        iB = np.where(near)[0][np.argmax(Sd[near])]
        order = np.argsort(Sd)[::-1]; peaks = []
        for i in order:
            if all(abs(Ps[i] - Ps[j]) / Ps[j] > 0.005 for j in peaks): peaks.append(i)
            if len(peaks) >= 10: break
        rankB = int(np.sum(Sd[peaks] > Sd[iB]) + 1)
        tls_out[var] = dict(best_P=float(rr.period), best_SDE=float(rr.SDE), B_P=float(Ps[iB]), B_SDE=float(Sd[iB]), B_rank=rankB,
                            top=[(float(Ps[i]), float(Sd[i])) for i in peaks], runtime_s=time.time() - t0,
                            odd_even_mismatch=float(rr.odd_even_mismatch) if rr.period and abs(rr.period - 11.133) < 0.05 else None)
        print("TLS", var, tls_out[var], flush=True)
    json.dump(tls_out, open("out/b02_s98_tls.json", "w"), indent=1)


if __name__ == "__main__":
    import sys
    if "tls" not in sys.argv:
        main_bls()
    if "bls" not in sys.argv:
        main_tls()
