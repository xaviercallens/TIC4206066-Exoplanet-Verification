"""(4) S6 momentum dumps vs 3.18-d transit times (manuscript: 'every transit falls 0.4-3.5 h from a desaturation')."""
from pixloc import *
out = {}
for key in ["S6_ffi", "S6_tspoc_lc"]:
    if key == "S6_tspoc_lc":
        d = fits.getdata(FILES[key], 1); t = np.array(d["TIME"], float); q = np.array(d["QUALITY"], np.int64)
    else:
        c = Cube(key); t, q = c.time, c.quality
    m = (q & 32) > 0
    tt = t[m & np.isfinite(t)]
    groups = np.split(tt, np.where(np.diff(tt) > 0.1)[0] + 1)
    dumps = np.array([g.mean() for g in groups])
    n, tcs = transit_centres(EPH["A"], np.nanmin(t), np.nanmax(t))
    rows = []
    for tc in tcs:
        j = np.argmin(np.abs(dumps - tc))
        # is the transit covered by data?
        cov = np.sum(np.abs(t[np.isfinite(t)] - tc) < EPH["A"]["T14"] / 2)
        rows.append(dict(tc=float(tc), nearest_dump=float(dumps[j]), dt_h=float((tc - dumps[j]) * 24), n_cad_in_transit=int(cov)))
    out[key] = dict(dump_times=dumps, dump_interval_d=np.diff(dumps), transits=rows)
    print(key, "dumps", np.round(dumps, 4), "intervals", np.round(np.diff(dumps), 3))
    for r in rows:
        print(f"   A tc={r['tc']:.4f} nearest dump {r['nearest_dump']:.4f} dt={r['dt_h']:+.2f} h  n_in={r['n_cad_in_transit']}")
e = EPH["A"]; ncyc = (np.array([r["tc"] for r in out["S6_ffi"]["transits"]]) - e["T0"]) / e["P"]
out["ephemeris_sigma_at_S6_min"] = (np.sqrt(0.0023 ** 2 + (ncyc * 5.8e-6) ** 2) * 1440).tolist()
out["note"] = ("FFI dump flags mark one 30-min FFI per dump, so dump times are known to +/-15 min; S6 DRN: dumps every "
               "3.125 d (https://archive.stsci.edu/missions/tess/doc/tess_drn/tess_sector_06_drn08_v02.pdf, Sec. 1.2); "
               "P_A = 3.1828 d drifts by only 1.39 h per orbit relative to the dump cadence.")
dump(out, "s6_dumps.json")
