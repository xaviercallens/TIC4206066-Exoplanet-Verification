# Read-only: which QUALITY bits occur in the S98 QLP file, and does the frozen search mask remove every flagged cadence?
import numpy as np, json
from astropy.io import fits
f = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/00_data/mast/hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits"
with fits.open(f, memmap=True) as h:
    d = h[1].data; q = np.array(d['QUALITY'], np.int64); t = np.array(d['TIME'], float)
    sysrm = np.array(d['SYS_RM_FLUX'], float); det = np.array(d['DET_FLUX'], float); err = np.array(d['DET_FLUX_ERR'], float)
    td = h[1].header['TIMEDEL']
bits = {int(b): int(((q >> b) & 1).sum()) for b in range(32) if ((q >> b) & 1).any()}
MASK = 1610619071
masked = (q & MASK) != 0
out = dict(n_cad=int(q.size), TIMEDEL_s=td * 86400, bits_present_count={f"bit{b} (2^{b}={1<<b})": n for b, n in bits.items()},
           n_flagged_any=int((q != 0).sum()), n_removed_by_frozen_mask=int(masked.sum()),
           flagged_but_not_removed_by_frozen_mask=int(((q != 0) & ~masked).sum()),
           mask_bits=[b for b in range(32) if (MASK >> b) & 1],
           n_good_q0_finite_sysrm=int(((q == 0) & np.isfinite(sysrm)).sum()),
           frac_finite_DET_FLUX_ERR_good=float(np.isfinite(err[q == 0]).mean()),
           median_dt_s=float(np.nanmedian(np.diff(t)) * 86400))
print(json.dumps(out, indent=1)); json.dump(out, open('s98_quality_bits.json', 'w'), indent=1)
