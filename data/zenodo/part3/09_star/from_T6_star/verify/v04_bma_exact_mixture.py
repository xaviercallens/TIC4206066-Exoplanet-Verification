"""Verifier check: exact weighted-mixture percentiles (no resampling) of the ARIADNE per-grid posteriors in
companion_check_tic4206066/sed_out/BMA.pkl (read-only), plus a check that 'lum' in weighted_average/weighted_samples
equals 4 pi R^2 sigma T^4 of the separately averaged teff/rad arrays (loss of T-R correlation)."""
import json, pickle, numpy as np
from pathlib import Path
H = Path(__file__).resolve().parent
P = Path("/Users/paulsmith/space <3/tess-jev/addenda/companion_check_tic4206066/sed_out/BMA.pkl")
class _S:
    def __init__(self, *a, **k): pass
    def __setstate__(self, s): self.__dict__["_s"] = s
class U(pickle.Unpickler):
    def find_class(self, m, n):
        try: return super().find_class(m, n)
        except Exception: return type(n, (_S,), {})
b = U(open(P, "rb")).load()
w = {k: float(v) for k, v in b["weights"].items()}
def wpct(par):
    xs, ws = [], []
    for g, o in b["originals"].items():
        a = np.asarray(o[par], float); xs.append(a); ws.append(np.full(a.size, w[g] / a.size))
    x = np.concatenate(xs); ww = np.concatenate(ws); i = np.argsort(x); x, ww = x[i], np.cumsum(ww[i]); ww /= ww[-1]
    q = lambda p: float(np.interp(p, ww, x))
    lo, med, hi = q(0.158655), q(0.5), q(0.841345)
    return dict(median=med, plus=hi - med, minus=med - lo, n_per_grid={g: int(len(o[par])) for g, o in b["originals"].items()})
out = {"weights": w, "exact_mixture": {p: wpct(p) for p in ("teff", "rad", "lum", "z")}}
sb = 5.670374419e-5; Rs = 6.957e10; Ls = 3.828e33
for k in ("weighted_average", "weighted_samples"):
    T = np.asarray(b[k]["teff"], float); R = np.asarray(b[k]["rad"], float); L = np.asarray(b[k]["lum"], float)
    Lc = 4 * np.pi * (R * Rs) ** 2 * sb * T ** 4 / Ls
    out[k + "_lum_equals_RT4_of_arrays"] = dict(max_rel_diff=float(np.max(np.abs(Lc / L - 1))), corr_T_R=float(np.corrcoef(T, R)[0, 1]))
for g, o in b["originals"].items():
    out.setdefault("per_grid_corr_T_R", {})[g] = float(np.corrcoef(np.asarray(o["teff"], float), np.asarray(o["rad"], float))[0, 1])
json.dump(out, open(H / "v04_bma_exact_mixture.json", "w"), indent=1); print(json.dumps(out, indent=1))
