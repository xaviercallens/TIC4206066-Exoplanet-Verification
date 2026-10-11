"""Audit 2 (timing of signal A, TIC 4206066): shared code.

Written from scratch for this audit. No code or derived table from earlier analyses is imported or read.
Inputs: raw MAST light curves downloaded by the data agent of this audit (00_data/mast), sha256-checked against
00_data/DATA_MANIFEST.json at load time.

Conventions: all times BTJD = BJD_TDB - 2457000 (all three products declare TIMESYS=TDB, TIMEREF=SOLARSYSTEM,
BJDREFI=2457000; checked in load_lc).
"""
import hashlib, json, os
import numpy as np
from astropy.io import fits
import batman
from wotan import flatten

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "00_data"))
OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- inputs
FILES = {
    6: ("mast/hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits", "PDCSAP_FLUX"),
    32: ("mast/tess2020324010417-s0032-0000000004206066-0200-s_lc.fits", "PDCSAP_FLUX"),
    98: ("mast/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits", "SYS_RM_FLUX"),
}
# lightkurve 2.6.0 TessQualityFlags.DEFAULT_BITMASK = 17087 (attitude tweak, safe mode, coarse point, earth point,
# argabrightening, desaturation, manual exclude, impulsive outlier, bad calibration) -- printed from the installed package.
SPOC_DEFAULT = 17087
# QLP S98: documented bits (QLP data-product description Sec. 2.3.3: 1,4,8,16,32,128,1024,2048,4096) plus the
# undocumented 2^29 / 2^30 bits that coincide with the Earth scattered-light intervals at the start of each S98 orbit
# (HEASARC S98 notes) -- this reproduces the manuscript's "scattered-light intervals removed" selection.
QLP_DOC = 1 | 2 | 4 | 8 | 16 | 32 | 128 | 1024 | 2048 | 4096
QLP_SL = (1 << 29) | (1 << 30)

# ---------------------------------------------------------------- ephemerides used for window placement / masking
# A: manuscript values (only used to centre windows; the search range is +-144 min, >> any ephemeris error).
P_A0, T0_A0 = 3.1827856, 3998.7375
# B: task / manuscript values; T14_B from the manuscript range 1.6-1.7 h.
P_B, T0_B, T14_B = 11.13274, 4000.865, 1.66 / 24

# ---------------------------------------------------------------- template (fixed shape)
# Shape: k = Rp/R* = 0.0212, b = 0.35, rho* = 2.51 rho_sun (manuscript Sec. 1-2; only the shape matters for timing).
# Quadratic LD: Claret (2017, A&A 600, A30; VizieR J/A+A/600/A30 table25, ATLAS, LSM) at Teff 4000 K, [M/H]=0,
# xi=2 km/s: logg 4.5 -> (0.4090, 0.2526), logg 5.0 -> (0.2604, 0.3657); linear interpolation to logg 4.63 gives
# (0.370, 0.282). Sensitivity to LD is tested separately (t3 --ld variants).
RHO_SUN = 1.98841e33 / (4.0 / 3.0 * np.pi * 6.957e10 ** 3)
G_CGS = 6.6743e-8
SHAPE = dict(k=0.0212, b=0.35, rho=2.51, u1=0.370, u2=0.282)


def a_over_rs(rho, P):
    return (rho * RHO_SUN * G_CGS * (P * 86400.0) ** 2 / (3 * np.pi)) ** (1.0 / 3.0)


class Template:
    """Unit-depth transit profile s(tau) (tau = t - tc in days), 1 at mid-transit for the unbinned model,
    and cadence-integrated (boxcar) versions for any exposure time."""

    def __init__(self, P=P_A0, k=SHAPE["k"], b=SHAPE["b"], rho=SHAPE["rho"], u1=SHAPE["u1"], u2=SHAPE["u2"]):
        pm = batman.TransitParams()
        pm.t0, pm.per, pm.rp = 0.0, P, k
        pm.a = a_over_rs(rho, P)
        pm.inc = np.degrees(np.arccos(b / pm.a))
        pm.ecc, pm.w, pm.u, pm.limb_dark = 0.0, 90.0, [u1, u2], "quadratic"
        self.pm = pm
        self.dtau = 2.0 / 86400.0                      # 2-s fine grid
        self.half = 0.30                               # template support +-0.3 d (transit is +-0.04 d)
        self.tau = np.arange(-self.half, self.half + self.dtau / 2, self.dtau)
        fl = batman.TransitModel(pm, self.tau).light_curve(pm)
        self.s = (1.0 - fl) / (1.0 - fl[np.argmin(np.abs(self.tau))])
        inn = self.s > 0
        self.t14 = float(self.tau[inn][-1] - self.tau[inn][0])
        full = self.s >= 0.999 * self.s[inn].min() if False else None
        self._binned = {}

    def binned(self, cad):
        key = round(cad * 86400.0)
        if key not in self._binned:
            n = max(1, int(round(cad / self.dtau)))
            if n > 1:
                ker = np.ones(n) / n
                sb = np.convolve(self.s, ker, mode="same")
            else:
                sb = self.s.copy()
            self._binned[key] = sb
        return self._binned[key]

    def __call__(self, trel, cad):
        return np.interp(trel, self.tau, self.binned(cad), left=0.0, right=0.0)


# ---------------------------------------------------------------- data
def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def manifest_sha(rel):
    m = json.load(open(os.path.join(DATA, "DATA_MANIFEST.json")))
    for e in m["files"]:
        if e["relpath"] == rel:
            return e["sha256"], e.get("url"), e.get("utc")
    return None, None, None


def load_lc(sector, s98_mode="clean"):
    """Returns dict with t, f (normalised raw flux), cad, provenance.
    s98_mode: 'clean' = documented QLP bits + 2^29/2^30 (manuscript selection); 'doc' = documented bits only."""
    rel, col = FILES[sector]
    p = os.path.join(DATA, rel)
    sh = sha256(p); msh, url, utc = manifest_sha(rel)
    assert sh == msh, f"sha256 mismatch for {rel}"
    with fits.open(p) as h:
        hd0, hd1 = h[0].header, h[1].header
        d = h[1].data
        t = np.asarray(d["TIME"], float); f = np.asarray(d[col], float)
        q = np.asarray(d["QUALITY"]).astype(np.int64)
        meta = {k: (hd1.get(k) if hd1.get(k) is not None else hd0.get(k)) for k in
                ("TIMESYS", "TIMEREF", "BJDREFI", "BJDREFF", "TIMEUNIT", "TIMEDEL", "CROWDSAP", "FLFRCSAP")}
    assert meta["TIMESYS"] == "TDB" and int(meta["BJDREFI"]) == 2457000
    if sector == 98:
        mask = QLP_DOC | (QLP_SL if s98_mode == "clean" else 0)
    else:
        mask = SPOC_DEFAULT
    good = np.isfinite(t) & np.isfinite(f) & (f > 0) & ((q & mask) == 0)
    t, f = t[good], f[good]
    f = f / np.median(f)
    cad = float(np.median(np.diff(t)))
    return dict(sector=sector, t=t, f=f, cad=cad, file=rel, col=col, sha256=sh, url=url, utc=utc, meta=meta,
                bitmask=int(mask), n_good=int(good.sum()), n_total=int(len(good)))


# ---------------------------------------------------------------- epochs / masks
def a_epochs(t, P=P_A0, T0=T0_A0):
    n0 = int(np.floor((t.min() - T0) / P)) - 1; n1 = int(np.ceil((t.max() - T0) / P)) + 1
    n = np.arange(n0, n1 + 1)
    return n, T0 + n * P


def b_times(t):
    n0 = int(np.floor((t.min() - T0_B) / P_B)) - 1; n1 = int(np.ceil((t.max() - T0_B) / P_B)) + 1
    return T0_B + np.arange(n0, n1 + 1) * P_B


def near(t, centres, hw):
    if len(centres) == 0:
        return np.zeros(len(t), bool)
    d = np.abs(t[:, None] - np.asarray(centres)[None, :])
    return (d < hw).any(axis=1)


# ---------------------------------------------------------------- configuration of the per-event measurement
CFG = dict(
    W=0.35,                  # half-width of the fit window (d)
    R=144.0,                 # half-range of the mid-time search (min)
    R_alt=60.0,              # manuscript-like search half-range (min), evaluated on the same chi2 curve
    step=0.25,               # grid step (min)
    chunk=1.1,               # half-width of the local detrending chunk (d)
    win_len=0.75,            # wotan biweight window (d)
    maskA=1.5,               # A mask half-width in detrending, in units of T14_A
    maskB=1.0,               # B mask half-width in detrending, in units of T14_B
    rmB=0.5 * T14_B + 0.4 / 24,  # B cadences removed from the fit window: |t - tcB| < T14_B/2 + 24 min
    clip=4.0,                # upward clipping threshold (robust sigma) inside the window
    cov_set=0.5, cov_tim=0.7, snr_tim=2.0,
)


def detrend_chunk(t, f, centre, maskcentres_A, t14A, extra_mask=None, cfg=CFG):
    c = np.abs(t - centre) < cfg["chunk"]
    tc, fc = t[c], f[c]
    m = near(tc, maskcentres_A, cfg["maskA"] * t14A) | near(tc, b_times(tc) if len(tc) else [], cfg["maskB"] * T14_B)
    if extra_mask is not None:
        m |= np.abs(tc - extra_mask) < cfg["maskA"] * t14A
    if len(tc) < 10:
        return None
    # segments (same break tolerance as wotan); a segment with < 5 unmasked points cannot be detrended -> dropped
    seg = np.r_[0, np.where(np.diff(tc) > 0.375)[0] + 1, len(tc)]
    keep = np.ones(len(tc), bool)
    for a, b in zip(seg[:-1], seg[1:]):
        if (~m[a:b]).sum() < 5:
            keep[a:b] = False
    tc, fc, m = tc[keep], fc[keep], m[keep]
    if len(tc) < 10:
        return None
    tr = flatten(tc, fc, method="biweight", window_length=cfg["win_len"], edge_cutoff=0, break_tolerance=0.375,
                 return_trend=True, mask=m)[1]
    bad = ~np.isfinite(tr)
    if bad.any():
        if (~bad).sum() < 3:
            return None
        tr[bad] = np.interp(tc[bad], tc[~bad], tr[~bad])
    return tc, fc / tr


class Scanner:
    """Least-squares scan of the mid-time: f = c0 + c1 x - D s(x - dt), x = t - centre, on a dt grid.
    The design depends only on (t - centre), so many flux vectors can be scanned at once (prayer bead, bootstrap)."""

    def __init__(self, x, cad, tpl, cfg=CFG):
        self.x = x; self.N = len(x)
        self.dts = np.arange(-cfg["R"], cfg["R"] + cfg["step"] / 2, cfg["step"]) / 1440.0
        S = tpl(x[None, :] - self.dts[:, None], cad)                     # (ndt, N)
        self.S = S
        N = float(len(x)); sx = x.sum(); sxx = (x * x).sum()
        sS = S.sum(1); sxS = S @ x; sSS = (S * S).sum(1)
        G = np.empty((len(self.dts), 3, 3))
        G[:, 0, 0] = N; G[:, 0, 1] = G[:, 1, 0] = sx; G[:, 1, 1] = sxx
        G[:, 0, 2] = G[:, 2, 0] = -sS; G[:, 1, 2] = G[:, 2, 1] = -sxS; G[:, 2, 2] = sSS
        self.Ginv = np.linalg.pinv(G, rcond=1e-12, hermitian=True)   # pinv: D = 0 where the shifted template misses all data
        self.alt = np.abs(self.dts * 1440.0) <= cfg["R_alt"] + 1e-9

    def chi2(self, F):
        """F: (N,) or (N,K). Returns chi2 (ndt[,K]) and coefficients (ndt,3[,K])."""
        one = F.ndim == 1
        if one:
            F = F[:, None]
        b0 = F.sum(0)[None, :].repeat(len(self.dts), 0)                  # (ndt,K)
        b1 = (self.x @ F)[None, :].repeat(len(self.dts), 0)
        b2 = -(self.S @ F)
        B = np.stack([b0, b1, b2], axis=1)                                # (ndt,3,K)
        C = np.einsum("dij,djk->dik", self.Ginv, B)
        chi = (F * F).sum(0)[None, :] - np.einsum("dik,dik->dk", B, C)
        if one:
            return chi[:, 0], C[:, :, 0]
        return chi, C

    @staticmethod
    def refine(dts, chi, j):
        """parabolic refinement of the grid minimum."""
        if 0 < j < len(dts) - 1:
            y0, y1, y2 = chi[j - 1], chi[j], chi[j + 1]
            den = y0 - 2 * y1 + y2
            if den > 0:
                off = 0.5 * (y0 - y2) / den
                return dts[j] + off * (dts[1] - dts[0])
        return dts[j]

    def best(self, chi, restrict=None):
        """index of minimum (optionally within a boolean restriction), refined time (d)."""
        if restrict is None:
            j = int(np.argmin(chi))
            return j, self.refine(self.dts, chi, j)
        idx = np.where(restrict)[0]
        j = int(idx[np.argmin(chi[idx])])
        # refine only if neighbours are inside the restriction
        if restrict[max(j - 1, 0)] and restrict[min(j + 1, len(chi) - 1)]:
            return j, self.refine(self.dts, chi, j)
        return j, self.dts[j]


def white_interval(dts, chin, j):
    """contiguous Delta chi2 <= 1 interval around the minimum of a normalised chi2 curve (chi2/s2).
    Returns (half-width in d, touches_edge)."""
    lvl = chin[j] + 1.0
    L = j
    while L > 0 and chin[L - 1] <= lvl:
        L -= 1
    Rr = j
    while Rr < len(chin) - 1 and chin[Rr + 1] <= lvl:
        Rr += 1
    edge = (L == 0) or (Rr == len(chin) - 1)
    # linear interpolation of crossings
    def cross(i_in, i_out):
        y0, y1 = chin[i_in], chin[i_out]
        if y1 == y0:
            return dts[i_in]
        return dts[i_in] + (lvl - y0) / (y1 - y0) * (dts[i_out] - dts[i_in])
    lo = cross(L, L - 1) if L > 0 else dts[0]
    hi = cross(Rr, Rr + 1) if Rr < len(chin) - 1 else dts[-1]
    return 0.5 * (hi - lo), edge


def measure(t, fd, centre, cad, tpl, rng, n_boot=200, pb_max=None, cfg=CFG):
    """Measure one event in an already-detrended window (t, fd around 'centre').
    Returns dict with best dt (R and R_alt), depth, white errors, prayer-bead and bootstrap errors."""
    x = t - centre
    sc = Scanner(x, cad, tpl, cfg)
    chi, C = sc.chi2(fd)
    j, dtb = sc.best(chi)
    ja, dta = sc.best(chi, sc.alt)
    N = len(x)
    s2 = chi[j] / max(N - 4, 1)
    chin = chi / s2
    ew, edge = white_interval(sc.dts, chin, j)
    ewa, edgea = white_interval(sc.dts[sc.alt], chin[sc.alt], int(ja - np.where(sc.alt)[0][0]))
    c = C[j]
    D = float(c[2]); eD = float(np.sqrt(s2 * sc.Ginv[j, 2, 2]))
    Da = float(C[ja][2]); eDa = float(np.sqrt(s2 * sc.Ginv[ja, 2, 2]))
    model = c[0] + c[1] * x - c[2] * sc.S[j]
    r = fd - model
    # prayer bead (cyclic shifts of residuals; Winn et al. 2008 "rosary bead")
    shifts = np.arange(1, N)
    if pb_max is not None and len(shifts) > pb_max:
        shifts = shifts[np.linspace(0, len(shifts) - 1, pb_max).astype(int)]
    Fp = model[:, None] + np.stack([np.roll(r, k) for k in shifts], axis=1)
    chip, _ = sc.chi2(Fp)
    jp = np.argmin(chip, axis=0)
    dtp = sc.dts[jp]
    jpa = np.where(sc.alt)[0][np.argmin(chip[sc.alt], axis=0)]
    dtpa = sc.dts[jpa]
    # iid residual bootstrap
    Fb = model[:, None] + r[rng.integers(0, N, size=(N, n_boot))]
    chib, _ = sc.chi2(Fb)
    dtbs = sc.dts[np.argmin(chib, axis=0)]
    return dict(dt=float(dtb * 1440), dt_alt=float(dta * 1440), D=D, eD_white=eD, snr=D / eD, D_alt=Da, snr_alt=Da / eDa,
                e_white=float(ew * 1440), edge=bool(edge), e_white_alt=float(ewa * 1440), edge_alt=bool(edgea),
                e_pb=float(np.std(dtp) * 1440), e_pb_alt=float(np.std(dtpa) * 1440), e_boot=float(np.std(dtbs) * 1440),
                rms=float(np.sqrt(s2)), n=int(N), chi_curve=None)


def prepare_window(lc, centre, maskA_centres, t14A, extra_mask=None, cfg=CFG):
    """detrend locally, cut the fit window, remove B cadences, clip upward outliers.
    Returns (t, fd, coverage, n_expected) or None."""
    out = detrend_chunk(lc["t"], lc["f"], centre, maskA_centres, t14A, extra_mask=extra_mask, cfg=cfg)
    if out is None:
        return None
    tc, fd = out
    w = np.abs(tc - centre) < cfg["W"]
    tb = b_times(tc)
    w &= ~near(tc, tb, cfg["rmB"])
    t, y = tc[w], fd[w]
    if len(t) < 5:
        return None
    oot = np.abs(t - centre) > 0.75 * t14A
    base = y[oot] if oot.sum() >= 5 else y
    med = np.median(base); sig = 1.4826 * np.median(np.abs(base - med))
    keep = y < med + cfg["clip"] * sig
    t, y = t[keep], y[keep]
    inn = np.abs(t - centre) < t14A / 2
    cov = float(inn.sum() * lc["cad"] / t14A)
    return t, y, cov


def n_min(lc, cfg=CFG):
    return int(0.4 * 2 * cfg["W"] / lc["cad"])


def binned_beta(t, r, cad, taus_min):
    """time-averaging beta (Winn et al. 2008, eq. 2): sigma_N expected = sigma_1/sqrt(N) sqrt(M/(M-1))."""
    s1 = np.std(r)
    out = {}
    for tau in taus_min:
        Nb = int(round(tau / (cad * 1440)))
        if Nb < 2:
            continue
        # bin consecutive points inside continuous segments only
        seg = np.r_[0, np.where(np.diff(t) > 1.5 * cad)[0] + 1, len(t)]
        means = []
        for a, b in zip(seg[:-1], seg[1:]):
            nseg = (b - a) // Nb
            if nseg < 1:
                continue
            means.append(r[a:a + nseg * Nb].reshape(nseg, Nb).mean(1))
        if not means:
            continue
        mm = np.concatenate(means); M = len(mm)
        if M < 10:
            continue
        exp = s1 / np.sqrt(Nb) * np.sqrt(M / (M - 1.0))
        out[float(Nb * cad * 1440)] = float(np.std(mm) / exp)
    return out
