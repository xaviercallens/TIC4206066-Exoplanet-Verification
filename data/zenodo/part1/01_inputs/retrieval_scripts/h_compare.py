"""Compare fresh downloads with earlier copies: byte identity (MAST, PRF), numeric identity (TESScut arrays), value identity (Gaia/TIC target rows)."""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astroquery.mast import Tesscut

A = "/Users/paulsmith/space <3/tess-jev/addenda"
OLD_MAST = [f"{A}/deep_followup_tic4206066/lc", f"{A}/tic4206066_PUBLICATION_AUDIT/data/mast"]
OLD_PIX = [f"{A}/deep_followup_tic4206066/pixels", f"{A}/tic4206066_PUBLICATION_AUDIT/data/tesscut"]
OLD_PRF = f"{A}/tic4206066_PUBLICATION_AUDIT/data/prf"
OLD_GAIA = f"{A}/tic4206066_PUBLICATION_AUDIT/data/gaia"
OLD_CAT = f"{A}/tic4206066_PUBLICATION_AUDIT/data/catalogs"
R = {"mast": {}, "tesscut": {}, "prf": {}, "gaia": {}, "tic": {}, "notes": []}

# --- MAST: byte identity by filename
for f in sorted(glob.glob(os.path.join(ROOT, "mast", "*.fits"))):
    s = sha256_file(f); bn = os.path.basename(f); r = {"sha256": s, "old": {}}
    for d in OLD_MAST:
        p = os.path.join(d, bn)
        r["old"][os.path.relpath(p, A)] = ("identical" if sha256_file(p) == s else "DIFFERENT") if os.path.exists(p) else "absent"
    R["mast"][bn] = r
    print(bn, r["old"], flush=True)

# --- TESScut: numeric identity of arrays (21x21 old copies); 31x31 inner 21x21 vs 21x21 fresh
def load_cut(p):
    with fits.open(p, memmap=False) as h:
        d = h[1].data
        out = {c: np.array(d[c]) for c in ("TIME", "TIMECORR", "CADENCENO", "RAW_CNTS", "FLUX", "FLUX_ERR", "FLUX_BKG", "FLUX_BKG_ERR", "QUALITY", "POS_CORR1", "POS_CORR2") if c in d.columns.names}
        hdr1 = h[1].header; hdr2 = h[2].header if len(h) > 2 else None
        wcs = {k: hdr2.get(k) for k in ("CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "1CRV4P", "2CRV4P", "1CRV5P", "2CRV5P", "CRVAL1P", "CRVAL2P")} if hdr2 is not None else {}
        wcs1 = {k: hdr1.get(k) for k in ("1CRV4P", "2CRV4P", "1CRV5P", "2CRV5P", "1CRPX4", "2CRPX4", "1CRVL4", "2CRVL4")}
        prim = {k: h[0].header.get(k) for k in ("DATE", "PROCVER", "DATA_REL", "RA_OBJ", "DEC_OBJ", "TSTART", "TSTOP")}
    return out, wcs, wcs1, prim

def cmp_arrays(a, b):
    res = {}
    for k in a:
        if k not in b:
            res[k] = "missing in old"; continue
        x, y = a[k], b[k]
        if x.shape != y.shape:
            res[k] = f"shape {x.shape} vs {y.shape}"; continue
        same_nan = np.array_equal(np.isnan(x), np.isnan(y)) if x.dtype.kind == "f" else True
        eq = np.array_equal(x, y, equal_nan=True) if x.dtype.kind == "f" else np.array_equal(x, y)
        if eq:
            res[k] = "identical"
        else:
            dd = np.nanmax(np.abs(x.astype(float) - y.astype(float)))
            res[k] = f"DIFFERENT max|diff|={dd:.3g} nan_pattern_same={same_nan}"
    return res

for f in sorted(glob.glob(os.path.join(ROOT, "tesscut", "*", "*.fits"))):
    rel = os.path.relpath(f, ROOT); size = rel.split("/")[1]; sec = int(os.path.basename(f).split("-")[1][1:])
    new, wcs, wcs1, prim = load_cut(f)
    r = {"primary": prim, "wcs_ext2": wcs, "old": {}}
    if size == "21x21":
        for d in OLD_PIX:
            for p in glob.glob(os.path.join(d, f"tess-s{sec:04d}-*_21x21_astrocut.fits")):
                old, owcs, owcs1, oprim = load_cut(p)
                r["old"][os.path.relpath(p, A)] = dict(arrays=cmp_arrays(new, old), old_primary=oprim,
                                                       wcs_same={k: (wcs.get(k) == owcs.get(k)) for k in wcs})
    else:  # 31x31: its central 21x21 should equal the fresh 21x21 cutout
        f21 = glob.glob(os.path.join(ROOT, "tesscut", "21x21", f"tess-s{sec:04d}-*_21x21_astrocut.fits"))[0]
        n21 = load_cut(f21)[0]
        inner = {k: (v[:, 5:26, 5:26] if v.ndim == 3 else v) for k, v in new.items()}
        r["inner21_vs_fresh21"] = cmp_arrays(inner, n21)
    R["tesscut"][rel] = r
    print(rel, json.dumps(r["old"] if size == "21x21" else r["inner21_vs_fresh21"], default=str)[:600], flush=True)

# sector availability as reported by TESScut now
try:
    st = Tesscut.get_sectors(coordinates=SkyCoord(RA, DEC, unit="deg"))
    R["tesscut_sectors_available"] = [dict(zip(st.colnames, [str(x) for x in row])) for row in st]
    R["tesscut_sectors_utc"] = utcnow()
except Exception as ex:
    R["tesscut_sectors_available"] = repr(ex)

# --- PRF: byte identity
for f in sorted(glob.glob(os.path.join(ROOT, "prf", "start_s0004", "*", "*.fits"))):
    cc = os.path.basename(os.path.dirname(f)); p = os.path.join(OLD_PRF, cc, os.path.basename(f))
    R["prf"][f"{cc}/{os.path.basename(f)}"] = ("identical" if sha256_file(p) == sha256_file(f) else "DIFFERENT") if os.path.exists(p) else "absent"
print("PRF", {v: list(R["prf"].values()).count(v) for v in set(R["prf"].values())})

# --- Gaia target row: value comparison
def cmp_rows(new, old, label):
    diffs, same, onlynew, near = {}, 0, [], []
    for c in new.index:
        if c not in old.index:
            onlynew.append(c); continue
        a, b = new[c], old[c]
        try:
            fa, fb = float(a), float(b)
            if (np.isnan(fa) and np.isnan(fb)) or fa == fb:
                same += 1
            elif abs(fa - fb) <= 1e-6 * max(1e-30, abs(fa)):  # float32 values written with fewer digits in the older CSV
                same += 1; near.append(c)
            else:
                diffs[c] = (a, b)
        except (TypeError, ValueError):
            if str(a).strip() == str(b).strip() or (str(a) in ("nan", "", "--", "None") and str(b) in ("nan", "", "--", "None")):
                same += 1
            else:
                diffs[c] = (str(a), str(b))
    return dict(label=label, n_same=same, n_same_within_rel1e6_only=len(near), within_rel1e6_cols=near, n_diff=len(diffs), diffs={k: [str(x) for x in v] for k, v in diffs.items()}, only_in_new=onlynew)

gn = pd.read_csv(os.path.join(ROOT, "catalogs", "gaiadr3_gaia_source_cone240as.csv"))
go = pd.read_csv(os.path.join(OLD_GAIA, "gaia_dr3_source.csv"))
tn = gn[gn.source_id == GAIA_DR3].iloc[0]; to = go[go.source_id == GAIA_DR3].iloc[0]
R["gaia"]["gaia_source_target"] = cmp_rows(tn, to, "gaia_source target row")
nb = pd.read_csv(os.path.join(OLD_GAIA, "gaia_dr3_neighbours_180as.csv"))
inner = set(gn[gn.dist_arcsec_j2016pos <= 180].source_id)
R["gaia"]["neighbours_180as"] = dict(fresh_n=len(inner), old_n=len(nb), only_fresh=sorted(map(int, inner - set(nb.source_id))), only_old=sorted(map(int, set(nb.source_id) - inner)))
for new, old in [("gaiadr3_astrophysical_parameters_target", "gaia_dr3_ap"), ("gaiadr3_nss_two_body_orbit_target", "gaia_dr3_nss_two_body"),
                 ("gaiadr3_nss_acceleration_astro_target", "gaia_dr3_nss_accel"), ("gaiadr3_nss_non_linear_spectro_target", "gaia_dr3_nss_nonlin"),
                 ("gaiadr3_nss_vim_fl_target", "gaia_dr3_nss_vim"), ("gaiadr3_vari_summary_target", "gaia_dr3_vari_summary")]:
    a = pd.read_csv(os.path.join(ROOT, "catalogs", new + ".csv")); b = pd.read_csv(os.path.join(OLD_GAIA, old + ".csv"))
    if len(a) and len(b):
        R["gaia"][new] = cmp_rows(a.iloc[0], b.iloc[0], new)
    else:
        R["gaia"][new] = dict(fresh_rows=len(a), old_rows=len(b))
print(json.dumps(R["gaia"], default=str)[:3000])

# --- TIC target row + cone membership
tn = pd.read_csv(os.path.join(ROOT, "catalogs", "tic82_cone300as.csv"))
to = pd.read_csv(os.path.join(OLD_CAT, "tic82_cone5arcmin.csv"))
R["tic"]["target_row"] = cmp_rows(tn[tn.ID == TIC].iloc[0], to[to.ID == TIC].iloc[0], "TIC target row")
R["tic"]["cone"] = dict(fresh_n=len(tn), old_n=len(to), only_fresh=sorted(map(int, set(tn.ID) - set(to.ID)))[:50], only_old=sorted(map(int, set(to.ID) - set(tn.ID)))[:50])
print(json.dumps(R["tic"], default=str)[:1500])
json.dump(R, open(os.path.join(ROOT, "logs", "comparison_with_existing.json"), "w"), indent=1, default=str)
