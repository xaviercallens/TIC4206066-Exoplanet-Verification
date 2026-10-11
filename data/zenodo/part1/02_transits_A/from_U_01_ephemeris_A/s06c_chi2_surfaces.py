"""Chi2(tc) surfaces for the problematic events (S32 E-573, E-566, S98 E-1, E10) in the baseline method and in
the S32 10-min FFI variant, to document bimodality."""
import numpy as np, json
from a_common import *
from a_timing import template, tpl_eval, segment_window, WIN_HALF
from a_global import b_in
x = np.load('template_params.npy'); T0g, Pg = x[0], x[1]
out = {}
for tag, sec, ep in [('S32_w0.75', 32, -573), ('S32ffi_w0.75', 32, -573), ('S32_w0.75', 32, -566), ('S32ffi_w0.75', 32, -566),
                     ('S98_w0.75_clean_SYS_RM_FLUX', 98, -1), ('S98_w0.75_clean_SYS_RM_FLUX', 98, 10)]:
    z = np.load(f'lc/{tag}.npz'); key = z['key'].item(); key = int(key) if key.isdigit() else key
    t = z['t']; y = z['f']; m = ~b_in(t); t = t[m]; y = y[m]; sig = float(z['sig_pp'])
    tp = T0g + ep * Pg; tpl = template(key, x)
    sel = segment_window(t, tp); tt = t[sel]; yy = y[sel]
    xs = (tt - tp) / WIN_HALF; Pm = np.vstack([np.ones_like(xs), xs])
    grid = np.arange(-144, 144.1, 1.0)
    chi = []
    for g in grid:
        mm = tpl_eval(tpl, tt - tp - g / 1440.)
        X = (Pm * mm[None, :]).T; b, *_ = np.linalg.lstsq(X, yy, rcond=None); r = yy - X @ b; chi.append(r @ r / sig ** 2)
    chi = np.array(chi); chi -= chi.min()
    from scipy.signal import argrelmin
    mins = argrelmin(chi, order=5)[0]
    out[f'{tag}|E{ep}'] = dict(local_minima=[(float(grid[i]), float(chi[i])) for i in mins if chi[i] < 15],
                               chi2_at_0=float(chi[np.argmin(np.abs(grid))]))
    print(tag, ep, out[f'{tag}|E{ep}'])
json.dump(out, open('out_s06c_chi2_surfaces.json', 'w'), indent=1)
