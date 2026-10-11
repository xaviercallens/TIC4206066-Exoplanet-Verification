"""Difference imaging + PRF localisation library (audit 03_pixels, written from scratch).

Method
------
* per transit event: all cadences with |t - tc| <= T14/2 + W (W = baseline half-width), other signal's
  transits masked (+/- 30 min), quality-flagged cadences removed.  For every pixel an ordinary least-squares fit
  f(t) = c0 + c1 (t-tc) [+ c2 (t-tc)^2] [+ steps] + D * box(t) is made, where box(t) is the fraction of each
  exposure inside [tc - T14/2, tc + T14/2].  D (sign flipped: positive = flux lost) is the per-event difference image.
* events are combined with inverse-variance weights per pixel.
* PRF fit: D = a * PRF(x, y) + b on a square region around the target (TESS_PRF, MAST start_s0004 files).
* reference position of the target: PRF scene fit of the out-of-transit image with all Gaia DR3 sources
  (positions propagated to the sector epoch, fluxes from TIC 8.2 Tmag) and a free global shift (dx, dy).
"""
import numpy as np
from scipy.optimize import minimize, least_squares
from common import *

PAD_OTHER = 0.5 / 24.0


class Dataset:
    def __init__(self, key, variant="default", W_hours=4.0, deg=1, fit_half=5, dump_exclude_h=0.0,
                 dump_step=False, bitmask=DEFAULT_BITMASK):
        self.key, self.variant = key, variant
        c = Cube(key)
        self.cube = c
        good = c.good(bitmask)
        self.stray_flag = (c.quality & STRAY_BIT) > 0
        # momentum dumps from the cube's own QUALITY bit 6 (value 32)
        self.dump_times = self._dump_times(c)
        if variant in ("clean", "clean_q30", "clean29"):
            # clean   : SPOC/TESScut straylight bit 2048 (+ QLP 2048 in S98)
            # clean29 : additionally QLP bit 2^29 (undocumented in the QLP docs; marks elevated-background
            #           intervals, median sky 672 vs 135 e-/s; it reproduces the manuscript's 3 clean + 1
            #           scattered-light 11.13 d transits in S98)
            # clean_q30: additionally QLP bit 2^30 ("QLP quality flag", tess.mit.edu/qlp)
            bad = self.stray_flag.copy()
            if c.sector == 98:
                bad |= self._qlp_flag(c.time, 2048)
                if variant == "clean_q30":
                    bad |= self._qlp_flag(c.time, 1 << 30)
                if variant == "clean29":
                    bad |= self._qlp_flag(c.time, 1 << 29)
            good &= ~bad
        if dump_exclude_h > 0:
            for td in self.dump_times:
                good &= np.abs(c.time - td) > dump_exclude_h / 24.0
        self.good = good
        self.t = c.time[good]
        self.F = c.flux[good]
        self.dt = c.timedel
        self.W = W_hours / 24.0
        self.deg = deg
        self.dump_step = dump_step
        self.epoch = c.mid_epoch_jyear()
        # Gaia scene at this epoch
        g = load_gaia(240.0)
        ra, dec = propagate(g, self.epoch)
        self.gaia = g
        self.g_ra, self.g_dec = ra, dec
        x, y = c.sky_to_pix(ra, dec)
        self.g_x, self.g_y = x, y
        it = np.where(g["source_id"].values == TARGET_GAIA)[0][0]
        self.it = it
        self.xt_wcs, self.yt_wcs = x[it], y[it]
        self.fit_half = fit_half
        xc, yc = int(round(self.xt_wcs)), int(round(self.yt_wcs))
        self.ys = slice(max(0, yc - fit_half), min(c.ny, yc + fit_half + 1))
        self.xs = slice(max(0, xc - fit_half), min(c.nx, xc + fit_half + 1))
        self.prf = prf_for_cube(c, self.xt_wcs, self.yt_wcs)
        self.medimg = np.nanmedian(self.F, axis=0)

    @staticmethod
    def _dump_times(c):
        m = (c.quality & 32) > 0
        if not m.any():
            return np.array([])
        t = c.time[m]
        groups = np.split(t, np.where(np.diff(t) > 0.1)[0] + 1)
        return np.array([np.nanmean(gg) for gg in groups])

    @staticmethod
    def _qlp_flag(t, bit):
        from astropy.io import fits as _f
        d = _f.getdata(FILES["QLP98"], 1)
        tq, q = np.array(d["TIME"], float), np.array(d["QUALITY"], np.int64)
        ok = np.isfinite(tq)
        tq, q = tq[ok], q[ok]
        idx = np.clip(np.searchsorted(tq, t), 1, len(tq) - 1)
        j = np.where(np.abs(tq[idx] - t) < np.abs(tq[idx - 1] - t), idx, idx - 1)
        close = np.abs(tq[j] - t) < 150.0 / 86400
        return close & ((q[j] & bit) > 0)

    # ------------------------------------------------------------------------------------------------------
    def other_mask(self, sig):
        other = "B" if sig == "A" else "A"
        return in_any_transit(self.t, self.dt, EPH[other], PAD_OTHER)

    def event_windows(self, sig, centres=None, min_in_frac=0.6, min_side_frac=0.3):
        """Return list of dicts (tc, idx of cadences used) for accepted events."""
        e = EPH[sig]
        T14 = e["T14"]
        if centres is None:
            _, centres = transit_centres(e, self.t.min() - 1, self.t.max() + 1)
        om = self.other_mask(sig)
        out = []
        for tc in centres:
            sel = (np.abs(self.t - tc) <= T14 / 2 + self.W) & ~om
            idx = np.where(sel)[0]
            if len(idx) < 5:
                continue
            frac = box_fraction(self.t[idx], self.dt, tc, T14)
            n_in_exp = T14 / self.dt
            pre = np.sum((self.t[idx] < tc - T14 / 2 - self.dt / 2))
            post = np.sum((self.t[idx] > tc + T14 / 2 + self.dt / 2))
            n_side_exp = self.W / self.dt
            ok = (frac.sum() >= min_in_frac * n_in_exp) and pre >= min_side_frac * n_side_exp \
                and post >= min_side_frac * n_side_exp
            out.append(dict(tc=tc, idx=idx, ok=bool(ok), in_sum=float(frac.sum()), n_in_exp=float(n_in_exp),
                            pre=int(pre), post=int(post)))
        return out

    def event_diff(self, idx, tc, sig, F=None, Fsel=None):
        """Per-pixel regression for one event; returns (D, var) over the full stamp.
        Fsel (optional) = replacement pixel data for the rows idx (used by the injections)."""
        if F is None:
            F = self.F
        t = self.t[idx]
        x = t - tc
        f = box_fraction(t, self.dt, tc, EPH[sig]["T14"])
        cols = [np.ones_like(x), x]
        if self.deg >= 2:
            cols.append(x ** 2)
        if self.dump_step and len(self.dump_times):
            for td in self.dump_times:
                if t.min() < td < t.max():
                    cols.append((t > td).astype(float))
        cols.append(f)
        X = np.vstack(cols).T
        Y = (F[idx] if Fsel is None else Fsel).reshape(len(idx), -1)
        XtX_inv = np.linalg.pinv(X.T @ X)
        beta = XtX_inv @ (X.T @ Y)
        resid = Y - X @ beta
        dof = max(len(idx) - X.shape[1], 1)
        s2 = (resid ** 2).sum(axis=0) / dof
        D = -beta[-1]
        var = s2 * XtX_inv[-1, -1]
        shp = F.shape[1:]
        return D.reshape(shp), var.reshape(shp)

    def combine(self, events, sig, F=None):
        Ds, Vs = [], []
        for ev in events:
            D, V = self.event_diff(ev["idx"], ev["tc"], sig, F, ev.get("Fsel"))
            Ds.append(D); Vs.append(V)
        Ds, Vs = np.array(Ds), np.array(Vs)
        # remove a per-event uniform offset (residual sky/scattered light), estimated as the median of the
        # fit-region pixels more than 3 px from the target; the PRF fit has a free constant anyway.
        yy, xx = np.indices(Ds.shape[1:])
        ring = np.zeros(Ds.shape[1:], bool)
        ring[self.ys, self.xs] = True
        ring &= np.hypot(xx - self.xt_wcs, yy - self.yt_wcs) > 3.0
        if ring.sum() >= 10:
            Ds = Ds - np.median(Ds[:, ring], axis=1)[:, None, None]
        w = 1.0 / Vs
        Dm = (w * Ds).sum(0) / w.sum(0)
        Vm = 1.0 / w.sum(0)
        # scatter of events about the mean (per pixel) -> reduced chi2
        chi2r = ((Ds - Dm) ** 2 * w).sum(0) / max(len(events) - 1, 1)
        return Dm, Vm, Ds, Vs, chi2r

    # ------------------------------------------------------------------------------------------------------
    def fit_region(self, img):
        return img[self.ys, self.xs]

    def model_region(self, x, y):
        full = self.prf.image(x, y, (self.cube.ny, self.cube.nx))
        return full[self.ys, self.xs]

    def chi2_at(self, D, sig, x, y):
        """Best (a, b) for fixed (x, y); returns chi2, a, b, sigma_a."""
        d = self.fit_region(D).ravel()
        s = sig.ravel()
        P = self.model_region(x, y).ravel()
        A = np.vstack([P, np.ones_like(P)]).T / s[:, None]
        bvec = d / s
        coef, *_ = np.linalg.lstsq(A, bvec, rcond=None)
        r = bvec - A @ coef
        cov = np.linalg.pinv(A.T @ A)
        return float(r @ r), float(coef[0]), float(coef[1]), float(np.sqrt(cov[0, 0]))

    def fit_position(self, D, V, x0=None, y0=None, search=2.5, step=0.5):
        sig = np.sqrt(self.fit_region(V))
        if x0 is None:
            x0, y0 = self.xt_wcs, self.yt_wcs
        best = (np.inf, x0, y0)
        for dx in np.arange(-search, search + 1e-9, step):
            for dy in np.arange(-search, search + 1e-9, step):
                c2 = self.chi2_at(D, sig, x0 + dx, y0 + dy)[0]
                if c2 < best[0]:
                    best = (c2, x0 + dx, y0 + dy)
        f = lambda p: self.chi2_at(D, sig, p[0], p[1])[0]
        r = minimize(f, [best[1], best[2]], method="Nelder-Mead",
                     options=dict(xatol=1e-3, fatol=1e-4, maxiter=400))
        c2, a, b, sa = self.chi2_at(D, sig, r.x[0], r.x[1])
        npix = sig.size
        return dict(x=float(r.x[0]), y=float(r.x[1]), chi2=c2, a=a, b=b, sig_a=sa, npix=int(npix), dof=int(npix - 4))

    # ------------------------------------------------------------------------------------------------------
    def scene_fit(self, img=None, half=7, maxdist_pix=12):
        """Fit the out-of-transit image with the Gaia scene: free global shift, target flux,
        neighbour flux scale, background plane.  Returns fitted target pixel position."""
        c = self.cube
        if img is None:
            img = self.medimg
        xc, yc = int(round(self.xt_wcs)), int(round(self.yt_wcs))
        ys = slice(max(0, yc - half), min(c.ny, yc + half + 1))
        xs = slice(max(0, xc - half), min(c.nx, xc + half + 1))
        sub = img[ys, xs]
        yy, xx = np.mgrid[ys, xs]
        Tt = self.gaia["Tmag_used"].values[self.it]
        rel = 10 ** (-0.4 * (self.gaia["Tmag_used"].values - Tt))
        near = np.where((np.hypot(self.g_x - self.xt_wcs, self.g_y - self.yt_wcs) < maxdist_pix)
                        & (np.arange(len(rel)) != self.it) & (rel > 1e-5))[0]
        shape = (c.ny, c.nx)

        def model(p):
            dx, dy, at, sn, b0, bx, by = p
            m = at * self.prf.image(self.xt_wcs + dx, self.yt_wcs + dy, shape)
            nb = np.zeros(shape)
            for j in near:
                nb += rel[j] * self.prf.image(self.g_x[j] + dx, self.g_y[j] + dy, shape)
            m = m + sn * at * nb
            return m[ys, xs] + b0 + bx * (xx - xc) + by * (yy - yc)

        err = np.sqrt(np.clip(sub, 1, None) + 30.0 ** 2)   # Poisson-like + floor (weights only)
        p0 = [0, 0, np.nansum(sub) - sub.size * np.nanmin(sub), 1.0, np.nanmin(sub), 0, 0]
        r = least_squares(lambda p: ((model(p) - sub) / err).ravel(), p0, x_scale="jac")
        res = dict(dx=float(r.x[0]), dy=float(r.x[1]), target_flux=float(r.x[2]), nb_scale=float(r.x[3]),
                   bkg=float(r.x[4]), n_neighbours_modelled=int(len(near)),
                   x_t=float(self.xt_wcs + r.x[0]), y_t=float(self.yt_wcs + r.x[1]),
                   frac_resid_rms=float(np.std(model(r.x) - sub) / np.max(sub)))
        # also a single-PRF (target only) fit on the central 5x5
        return res

    def sky_offset(self, x, y, xref, yref):
        J = pix_jacobian_arcsec(self.cube, xref, yref)
        return J @ np.array([x - xref, y - yref]), J


def mahal(mu, C, r):
    d = np.asarray(mu) - np.asarray(r)
    return float(np.sqrt(d @ np.linalg.solve(C, d)))
