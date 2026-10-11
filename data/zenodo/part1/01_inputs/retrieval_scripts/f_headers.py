"""(f) Key header info for every downloaded FITS product + momentum-dump (QUALITY bit 6, value 32) event list per sector."""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import numpy as np
from astropy.io import fits

KEYS = ["TELESCOP", "ORIGIN", "CREATOR", "PROCVER", "FILEVER", "DATA_REL", "DATE", "OBJECT", "TICID", "TICVER", "SECTOR", "CAMERA", "CCD",
        "RA_OBJ", "DEC_OBJ", "RA_ORIG", "DEC_ORIG", "PMRA", "PMDEC", "TESSMAG", "TEFF", "RADIUS",
        "TSTART", "TSTOP", "DATE-OBS", "DATE-END", "TIMESYS", "TIMEREF", "TASSIGN", "TIMEUNIT", "TELAPSE", "TELAPS",
        "BJDREFI", "BJDREFF", "BJDREFR", "INT_TIME", "EXPOSURE", "TIMEDEL", "FRAMETIM", "NUM_FRM", "DEADC", "TIMEPIXR", "TIERRELA",
        "CROWDSAP", "FLFRCSAP", "PDCVAR", "PDCMETHD", "CDPP0_5", "CDPP1_0", "CDPP2_0", "FFI_TYPE", "CALIB", "BESTAP", "SMALLAP", "LARGEAP",
        "DETMET", "APERTURE", "VERSION", "PIPELINE", "HLSPVER", "HLSPID", "HLSPLEAD", "LICENSE", "CRMITEN", "CRBLKSZ", "CRSPOC",
        "1CRV4P", "2CRV4P", "1CRV5P", "2CRV5P", "NAXIS2"]

def hdr_info(path):
    info = {}
    with fits.open(path, memmap=True) as hl:
        info["hdus"] = [h.name for h in hl]
        for i, h in enumerate(hl[:3]):
            for k in KEYS:
                if k in h.header and k not in info:
                    v = h.header[k]
                    info[k] = None if isinstance(v, fits.card.Undefined) else v
            # any other version-like keywords (QLP etc.)
            for k in h.header:
                if ("VER" in k or "PIPE" in k) and k not in info and not k.startswith(("TTYPE", "TFORM")):
                    v = h.header[k]
                    info[k] = None if isinstance(v, fits.card.Undefined) else v
        if len(hl) > 1 and isinstance(hl[1], fits.BinTableHDU):
            info["columns"] = hl[1].columns.names
            info["n_rows"] = len(hl[1].data)
    return info

def dump_events(path, label):
    with fits.open(path, memmap=True) as hl:
        d = hl[1].data
        cols = [c.upper() for c in hl[1].columns.names]
        tcol = "TIME" if "TIME" in cols else None
        qcol = next((c for c in ("QUALITY", "QUAL_FLAG", "DQUALITY") if c in cols), None)
        if tcol is None or qcol is None:
            return None
        t = np.array(d[tcol], float); q = np.array(d[qcol]).astype(np.int64)
        hdr = {**dict(hl[0].header), **dict(hl[1].header)}
    m = (q & 32) != 0
    ok = m & np.isfinite(t)
    tt = np.sort(t[ok])
    ev = []
    if tt.size:
        grp = np.split(tt, np.where(np.diff(tt) > 0.2)[0] + 1)
        for g in grp:
            ev.append(dict(t_first_btjd=round(float(g[0]), 6), t_last_btjd=round(float(g[-1]), 6),
                           t_mid_btjd=round(float(0.5 * (g[0] + g[-1])), 6), n_flagged_cadences=int(g.size)))
    # all bits used
    bits = {int(b): int(((q & (1 << b)) != 0).sum()) for b in range(0, 20) if ((q & (1 << b)) != 0).any()}
    return dict(product=label, sector=hdr.get("SECTOR"), quality_column=qcol, n_cadences=int(t.size),
                n_bit32_flagged=int(m.sum()), n_bit32_flagged_nan_time=int((m & ~np.isfinite(t)).sum()),
                bits_set_counts={f"bit{b+1}_value{1<<b}": n for b, n in bits.items()}, events=ev)

files = sorted(glob.glob(os.path.join(ROOT, "mast", "*.fits")) + glob.glob(os.path.join(ROOT, "tesscut", "*", "*.fits")))
H, DUMPS = {}, []
for f in files:
    rel = os.path.relpath(f, ROOT)
    H[rel] = hdr_info(f)
    if "_dvt" in f:
        continue
    r = dump_events(f, rel)
    if r:
        DUMPS.append(r)
    print(rel, H[rel].get("SECTOR"), H[rel].get("CAMERA"), H[rel].get("CCD"), (r or {}).get("n_bit32_flagged"), len((r or {}).get("events", [])), flush=True)
json.dump(H, open(os.path.join(ROOT, "logs", "fits_headers.json"), "w"), indent=1, default=str)
json.dump(DUMPS, open(os.path.join(ROOT, "logs", "momentum_dumps_from_quality.json"), "w"), indent=1, default=str)
# flat CSV of events
import csv
with open(os.path.join(ROOT, "logs", "momentum_dumps_from_quality.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["sector", "product", "event", "t_first_btjd", "t_last_btjd", "t_mid_btjd", "t_mid_bjd_tdb", "n_flagged_cadences"])
    for r in DUMPS:
        for i, e in enumerate(r["events"]):
            w.writerow([r["sector"], r["product"], i + 1, e["t_first_btjd"], e["t_last_btjd"], e["t_mid_btjd"], round(e["t_mid_btjd"] + 2457000, 6), e["n_flagged_cadences"]])
