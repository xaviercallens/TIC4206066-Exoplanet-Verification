"""Verifier: where the post-egress baseline of S98 E-1, E1, E10, E14 and S32 E-573 ends (raw MAST LCs, read-only)."""
import numpy as np, json
from astropy.io import fits
D = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/00_data/mast/"
PA, T0A = 3.1827850, 3998.73599   # adopted 19-set ephemeris
T14 = 1.887/24; PB, T0B, T14B = 11.13274, 4000.865, 1.66/24
out = {}
for s, fn, col in ((98, "hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits", "SYS_RM_FLUX"),
                   (32, "tess2020324010417-s0032-0000000004206066-0200-s_lc.fits", "PDCSAP_FLUX")):
    d = fits.open(D+fn)[1].data
    t = np.array(d["TIME"], float); f = np.array(d[col], float); q = np.array(d["QUALITY"]).astype(np.int64)
    good = np.isfinite(t) & np.isfinite(f) & ((q == 0) if s == 98 else ((q & 17087) == 0))
    nB = np.round((t - T0B)/PB); inB = np.abs(t - (T0B + nB*PB)) < T14B/2 + 0.4/24
    for n in ((-1, 0, 1, 10, 14, 13) if s == 98 else (-573, -572, -568)):
        tc = T0A + n*PA; eg = tc + T14/2; ig = tc - T14/2
        w = good & ~inB & (np.abs(t - tc) < 0.3)
        tw = t[w]
        post = tw[tw > eg]; pre = tw[tw < ig]
        # contiguous post-egress data: until first gap > 3 cadences
        cad = np.median(np.diff(t[good]))
        if len(post):
            gaps = np.where(np.diff(np.r_[eg, post]) > 3*cad)[0]
            contig_end = (np.r_[eg, post][gaps[0]] if len(gaps) else post[-1])
        else:
            contig_end = eg
        # reason for the end of post-egress data within 0.3 d
        nxt = t[(t > contig_end) & (np.abs(t - tc) < 0.3)]
        reason = {}
        if len(nxt):
            k = nxt[:10]; qq = q[np.isin(t, k)]; fin = np.isfinite(f[np.isin(t, k)])
            bb = inB[np.isin(t, k)]
            reason = dict(next_cadences_quality=[int(x) for x in qq[:5]], next_in_B_mask=[bool(x) for x in bb[:5]],
                          next_finite=[bool(x) for x in fin[:5]])
        out[f"S{s}_E{n}"] = dict(tc=round(tc, 4), contiguous_post_egress_min=round((contig_end - eg)*1440, 1),
                                 n_post_cadences_within_0p3d=int(len(post)), last_post_min=round((post[-1]-eg)*1440, 1) if len(post) else None,
                                 pre_baseline_min=round((ig - pre[0])*1440, 1) if len(pre) else 0, **reason)
for k, v in out.items(): print(k, v)
json.dump(out, open(__file__.replace(".py", ".json"), "w"), indent=1)
