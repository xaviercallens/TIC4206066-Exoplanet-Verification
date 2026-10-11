"""T6_star rerun 3: (a) eRASS1 counterpart separation (C028) and (b) what the eRASS1 X-ray level implies for the
rotation period under the Wright et al. (2011, ApJ 743, 48) rotation-activity relation (C025).

Inputs (read-only):
  archival_extra_tic4206066/spectra_activity/raw/vizier_5as/J_A+A_682_A34_erass1-m.csv (eRASS1 main row, _r = 4.771)
  Gaia DR3 position/PM (opus55 06_stellar/data/gaia_dr3_source.ecsv)
  L_X from audit FINAL_PARAMETERS.json star.activity.xray (flux 3.3686e-14 cgs, 0.2-2.3 keV) and L_bol.
Wright+2011 numbers (arXiv:1109.4634, fetched 2026-09-25 ~15:35 UTC): log R_X,sat = -3.13 +- 0.08, Ro_sat = 0.13 +- 0.02,
beta = -2.18 +- 0.16 (all stars) or -2.70 +- 0.13 (unbiased subsample), rms ~0.3 dex in R_X; turnover time
Eq. 11 log tau = 1.16 - 1.49 log M - 0.54 log^2 M; Table 2 bin V-Ks 3.36-3.68: log tau = 1.67 (beta=-2) / 1.49 (beta=-2.70).
Output: out/t6_03_xray_rossby.json
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table

ROOT = Path("/Users/paulsmith/space <3/tess-jev")
H = Path(__file__).resolve().parent
x = pd.read_csv(ROOT / "addenda/archival_extra_tic4206066/spectra_activity/raw/vizier_5as/J_A+A_682_A34_erass1-m.csv").iloc[0]
g = Table.read(ROOT / "addenda/opus55_final_audit_tic4206066/06_stellar/data/gaia_dr3_source.ecsv")[0]
X = SkyCoord(x.RA_ICRS, x.DE_ICRS, unit="deg")
G16 = SkyCoord(g["ra"], g["dec"], unit="deg")
dt = (x.MJD - 57388.5) / 365.25
Gx = SkyCoord(g["ra"] + g["pmra"] * dt / 3.6e6 / np.cos(np.radians(g["dec"])), g["dec"] + g["pmdec"] * dt / 3.6e6, unit="deg")
res = {"sep": {"vizier__r_arcsec": float(x._r), "sep_to_gaia_J2016_arcsec": float(X.separation(G16).arcsec),
               "sep_to_gaia_at_eRASS1_epoch_arcsec": float(X.separation(Gx).arcsec), "posErr_arcsec": float(x.posErr),
               "DetLike0": float(x.DetLike0), "counts": float(x.MLcts1), "exposure_s": float(x.MLExp1)}}

fa = json.load(open(ROOT / "addenda/tic4206066_PUBLICATION_AUDIT/FINAL_PARAMETERS.json"))["star"]["activity"]["xray"]
Lsun = 3.828e33; pc = 3.0857e18; d = 35.577
LX = fa["flux_cgs"] * 4 * np.pi * (d * pc) ** 2
out = {}
for Lb in (0.0884, 0.0901):
    for fx_lab, fx in (("eRASS1", fa["flux_cgs"]), ("eRASS1_minus1sig", fa["flux_cgs"] - fa["flux_lower_err_cgs"]),
                       ("eRASS1_plus1sig", fa["flux_cgs"] + 2.0e-14), ("eRASS3", 1.73e-14)):
        RX = np.log10(fx * 4 * np.pi * (d * pc) ** 2 / (Lb * Lsun))
        row = {"log_RX": float(RX)}
        for beta in (-2.18, -2.70):
            ro = 0.13 * 10 ** ((RX - (-3.13)) / beta)
            row[f"Ro_beta{beta}"] = float(ro)
            for tlab, tau in (("tau_eq11_M0.606", 10 ** (1.16 - 1.49 * np.log10(0.606) - 0.54 * np.log10(0.606) ** 2)),
                              ("tau_tab2", 10 ** (1.67 if beta == -2.18 else 1.49))):
                row[f"P_beta{beta}_{tlab}"] = float(ro * tau)
                # predicted log R_X for an 8.0 d rotator and deficit in units of the 0.3 dex scatter
                rp = -3.13 + beta * np.log10(max(8.0 / tau, 0.13) / 0.13)
                row[f"pred_logRX_P8d_beta{beta}_{tlab}"] = float(rp)
                row[f"deficit_sigma0.3_P8d_beta{beta}_{tlab}"] = float((rp - RX) / 0.3)
        out[f"Lbol{Lb}_{fx_lab}"] = row
res["LX_erg_s"] = float(LX)
res["rossby"] = out
res["tau_eq11_M0.606_d"] = float(10 ** (1.16 - 1.49 * np.log10(0.606) - 0.54 * np.log10(0.606) ** 2))
res["caveats"] = ["eRASS 0.2-2.3 keV band vs ROSAT 0.1-2.4 keV used by Wright+2011 (band conversion not applied; K/M coronae give "
                  "ROSAT-band fluxes typically somewhat higher, which would raise R_X by <~0.1-0.2 dex)",
                  "single marginal detection (DET_LIKE 8.0, 4.5 counts); eRASS:3 flux ~2x lower",
                  "relation scatter ~0.3 dex in R_X"]
json.dump(res, open(H / "out/t6_03_xray_rossby.json", "w"), indent=1)
print(json.dumps(res, indent=1))
