"""T6_star rerun 5: Lindegren+2021 parallax zero point for the target (gaiadr3-zeropoint 0.1.0 in opus55 audit_env,
used read-only) and the resulting distance. Inputs from opus55 06_stellar/data/gaia_dr3_source.ecsv (read-only)."""
import json
from pathlib import Path
import numpy as np
from astropy.table import Table
from zero_point import zpt
H = Path(__file__).resolve().parent
g = Table.read("/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/06_stellar/data/gaia_dr3_source.ecsv")[0]
zpt.load_tables()
z = float(zpt.get_zpt(np.array([float(g["phot_g_mean_mag"])]), np.array([float(g["nu_eff_used_in_astrometry"])]),
                      np.array([np.nan]), np.array([float(g["ecl_lat"])]), np.array([int(g["astrometric_params_solved"])]))[0])
p, e = float(g["parallax"]), float(g["parallax_error"])
pc = p - z
out = dict(parallax_raw=p, parallax_error=e, astrometric_params_solved=int(g["astrometric_params_solved"]), zp=z,
           parallax_corrected=pc, parallax_corrected_err_stat=e, parallax_corrected_err_with_0p010_floor=float(np.hypot(e, 0.010)),
           d_corrected_pc=1000 / pc, d_err_stat=1000 * e / pc**2, d_err_with_floor=1000 * np.hypot(e, 0.010) / pc**2,
           d_raw_pc=1000 / p, d_raw_err=1000 * e / p**2)
json.dump(out, open(H / "out/t6_05_parallax_zp.json", "w"), indent=1); print(json.dumps(out, indent=1))
