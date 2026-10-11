"""Audit 2, step 2c/3: injection-recovery of the fixed-shape template at random out-of-transit times.

Mirrors the real measurement exactly: injection into the normalised raw flux (before detrending), local wotan
detrending with the injected event masked like a real event, same window, same B removal, same clipping, same
coverage / point-count rules, same search grid (+-144 min, and the +-60 min sub-range), same selection flags.
Injection times are uniform over the sector; the only exclusion is that the fit window must not contain a real
A transit (|t_inj - t_A| > W + T14/2 + 1 h). Bad data are NOT avoided (quality-flagged cadences are simply absent,
as for the real events). Depth = the sector's weighted-mean real depth (t3 base).

Usage: python t4_injections.py <sector> <n_accept> [seed]
"""
import sys, json, time
import numpy as np, pandas as pd
import t2_common as C

sector = int(sys.argv[1]); n_acc = int(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else sector
rng = np.random.default_rng(1000 + seed)
cfg = C.CFG
tpl = C.Template(); T14 = tpl.t14
ev = pd.read_csv(f"{C.OUT}/t3_events_base.csv")
meta = json.load(open(f"{C.OUT}/t3_meta_base.json"))
P, T0 = meta["ephem_window"]["P"], meta["ephem_window"]["T0"]
e = ev[(ev.sector == sector) & ev.sel]   # sector depth from the selected (timing-set) events
w = 1 / e.eD_white ** 2
Dinj = float(np.sum(w * e.D) / w.sum())
lc = C.load_lc(sector)
t, f0 = lc["t"], lc["f"]
_, tcsA = C.a_epochs(t, P, T0)
excl = cfg["W"] + T14 / 2 + 1 / 24
rows = []; tries = 0; t_start = time.time()
out = f"{C.OUT}/t4_inj_S{sector}.csv"
while len(rows) < n_acc and tries < 50 * n_acc:
    tries += 1
    ti = rng.uniform(t.min(), t.max())
    if np.min(np.abs(ti - tcsA)) < excl:
        continue
    ch = np.abs(t - ti) < cfg["chunk"] + 0.05
    f = f0.copy()
    f[ch] = f0[ch] * (1.0 - Dinj * tpl(t[ch] - ti, lc["cad"]))
    lcm = dict(lc, f=f)
    prep = C.prepare_window(lcm, ti, tcsA, T14, extra_mask=ti)
    if prep is None:
        continue
    tt, yy, cov = prep
    if cov < cfg["cov_set"] or len(tt) < C.n_min(lc):
        continue
    r = C.measure(tt, yy, ti, lc["cad"], tpl, rng, n_boot=100, pb_max=150)
    r.pop("chi_curve")
    r.update(sector=sector, t_inj=float(ti), coverage=cov, D_inj=Dinj)
    rows.append(r)
    if len(rows) % 100 == 0:
        pd.DataFrame(rows).to_csv(out, index=False)
        print(f"S{sector}: {len(rows)} accepted / {tries} tries, {time.time()-t_start:.0f} s", flush=True)
df = pd.DataFrame(rows)
df["sel"] = (df.coverage >= cfg["cov_tim"]) & (df.snr >= cfg["snr_tim"]) & (~df.edge)
df.to_csv(out, index=False)
print("done", sector, len(df), tries, "Dinj", Dinj)
