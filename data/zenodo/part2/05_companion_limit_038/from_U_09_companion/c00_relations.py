"""Audit 9: empirical relations used throughout (own implementation; coefficients verified in the primary papers):
 * Mann et al. 2019 (ApJ 871, 63), Table 6, n=5: log10 M = sum a_i (M_Ks - 7.5)^i, valid 4.0 < M_Ks < 11.0
   (0.075-0.70 Msun). a = [-0.642, -0.208, -8.43e-4, 7.87e-3, 1.42e-4, -2.13e-4].
 * Mann et al. 2015 (ApJ 804, 64), Table 1: R = 1.9515 - 0.3520 M_Ks + 0.01680 M_Ks^2 (valid 4.6 < M_Ks < 9.8).
 * Pecaut & Mamajek 2013 / Mamajek table v2022.04.16 (EEM_dwarf_UBVIJHK_colors_Teff.txt, downloaded 2026-09-25):
   M_G, Bp-Rp, G-Rp, Teff as a function of M_Ks (interpolated).
 * Stassun et al. 2019 (AJ 158, 138) Eq. 1: T = G - 0.00522555 c^3 + 0.0891337 c^2 - 0.633923 c + 0.0324473,
   c = Bp-Rp (valid -0.2 < c < 3.5; extrapolated beyond, degraded for M dwarfs).
"""
from pathlib import Path
import numpy as np

EEM = Path("/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/c09docs/EEM_dwarf.txt")
A19 = np.array([-0.642, -0.208, -8.43e-4, 7.87e-3, 1.42e-4, -2.13e-4])

def mass_mann19(mks):
    x = np.asarray(mks) - 7.5
    return 10 ** sum(a * x ** i for i, a in enumerate(A19))

def radius_mann15(mks):
    m = np.asarray(mks)
    return 1.9515 - 0.3520 * m + 0.01680 * m ** 2

_mk_grid = np.linspace(3.5, 11.5, 8001)
_m_grid = mass_mann19(_mk_grid)            # monotonically decreasing with M_Ks
def mks_of_mass(m):
    return np.interp(np.asarray(m), _m_grid[::-1], _mk_grid[::-1])

def load_eem():
    rows, cols = [], None
    for line in EEM.read_text().splitlines():
        if line.startswith("#SpT") and cols is None:
            cols = line[1:].split()
        elif cols is not None and line[:1] in "OBAFGKML" and line.split()[0].endswith("V"):
            parts = line.split()
            if len(parts) >= len(cols):
                rows.append(parts[:len(cols)])
        elif cols is not None and rows and line.startswith("#SpT"):
            break
    def col(name):
        i = cols.index(name)
        return np.array([np.nan if r[i] in ("...", "....") else float(r[i].rstrip(":")) for r in rows])
    spt = [r[0] for r in rows]
    return spt, {c: col(c) for c in ["Teff", "R_Rsun", "M_G", "Bp-Rp", "G-Rp", "M_Ks", "Msun", "M_J"]}

SPT, EEMC = load_eem()

def _interp_vs_mks(name, mks):
    x, y = EEMC["M_Ks"], EEMC[name]
    ok = np.isfinite(x) & np.isfinite(y)
    o = np.argsort(x[ok])
    return np.interp(np.asarray(mks), x[ok][o], y[ok][o])

def tmag_stassun(G, bprp):
    c = np.asarray(bprp)
    return G - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473

def companion_props(m):
    """For dwarf mass m [Msun]: M_Ks, R (Mann15), Teff, M_G, Bp-Rp, M_RP, M_T (Stassun Eq.1) from the EEM sequence."""
    mks = mks_of_mass(m)
    R = radius_mann15(mks)
    Teff = _interp_vs_mks("Teff", mks)
    MG = _interp_vs_mks("M_G", mks)
    bprp = _interp_vs_mks("Bp-Rp", mks)
    grp = _interp_vs_mks("G-Rp", mks)
    return dict(M=np.asarray(m), MKs=mks, R=R, Teff=Teff, MG=MG, bprp=bprp, MRP=MG - grp, MT=tmag_stassun(MG, bprp))

if __name__ == "__main__":
    print("target check: M(4.94)=%.4f  R(4.94)=%.4f" % (mass_mann19(4.94), radius_mann15(4.94)))
    for m in [0.08, 0.1, 0.15, 0.2, 0.3, 0.38, 0.45, 0.5, 0.54, 0.606]:
        p = companion_props(m)
        print({k: round(float(v), 3) for k, v in p.items()})
