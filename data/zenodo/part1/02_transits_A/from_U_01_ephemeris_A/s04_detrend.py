"""Step 2: detrend each sector (wotan biweight, A and B transits masked) for several variants; save npz.
Also a synthetic test that wotan's `mask` really excludes masked points from the trend."""
import numpy as np, json, os
from a_common import *
os.makedirs('lc', exist_ok=True)

# --- wotan mask sanity test ---
rng = np.random.default_rng(1)
t = np.arange(0, 10, 200 / 86400.); f = 1 + 3e-4 * np.sin(2 * np.pi * t / 4) + rng.normal(0, 8e-4, t.size)
inj = np.where(np.abs(t - 5) < T14A / 2, 1 - 5e-4, 1.0); fi = f * inj
m = np.abs(t - 5) < T14A
flat_m, _ = detrend(t, fi, m); flat_n, _ = detrend(t, fi, np.zeros_like(m))
w = np.abs(t - 5) < T14A / 2
test = dict(true_depth_ppm=500, depth_masked_ppm=float((1 - flat_m[w].mean()) * 1e6 - (1 - flat_m[~m].mean()) * 1e6),
            depth_unmasked_ppm=float((1 - flat_n[w].mean()) * 1e6 - (1 - flat_n[~m].mean()) * 1e6))
print('wotan mask test', test)

variants = {}
for key in [6, 32, '32ffi', 98]:
    for win in [0.5, 0.75, 1.0]:
        opts = [('clean', 'SYS_RM_FLUX')]
        if key == 98:
            opts += [('keepSL', 'SYS_RM_FLUX'), ('no30', 'SYS_RM_FLUX'), ('clean', 'SAP_FLUX'), ('clean', 'DET_FLUX')]
        for msk, col in opts:
            d = prepare(key, window=win, s98_mask=msk, s98_flux=col)
            tag = f'S{key}_w{win}_{msk}_{col}' if key == 98 else f'S{key}_w{win}'
            ok = d['ok']
            np.savez(f'lc/{tag}.npz', t=d['t'][ok], flat=d['flat'][ok], f=d['f'][ok], trend=d['trend'][ok],
                     inB=d['inB'][ok], exptime=d['exptime'], dump_t=d['dump_t'], sig_pp=d['sig_pp'], key=str(key))
            variants[tag] = dict(n=int(ok.sum()), sig_pp_ppm=round(d['sig_pp'] * 1e6, 1), rms_oot_ppm=round(d['rms_oot'] * 1e6, 1))
            print(tag, variants[tag])
json.dump(dict(wotan_mask_test=test, variants=variants), open('out_s04_detrend.json', 'w'), indent=1)
