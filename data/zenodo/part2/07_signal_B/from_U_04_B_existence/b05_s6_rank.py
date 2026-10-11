"""Step 5: Sector 6.  Are the two predicted B dips the two strongest 1.6-h dips outside A?  Several photometries and
detrendings; independent boxes by greedy non-overlapping selection; rank-test p-value; proximity to dumps/gaps;
background and centroid behaviour at the two events."""
import numpy as np, json
from math import comb
import common_b as c
import archdata as ad

TB = [1473.7330, 1484.8657]          # manuscript ephemeris (P 11.13274, T0 4000.865) propagated to S6
PHOTS = ["pipe", "sap", "qlp", "eleanor", "tglc_aper", "tglc_psf", "own", "own05"]
DETS = [("biweight", 0.5), ("biweight", 0.75), ("biweight", 1.0), ("median", 0.75)]


def independent_boxes(g):
    x, s = g.snr_series()
    ok = np.isfinite(s)
    x, s = x[ok], s[ok]
    order = np.argsort(s)[::-1]
    taken = []
    blocked = np.zeros(len(x), bool)
    for i in order:
        if blocked[i]:
            continue
        taken.append((float(x[i]), float(s[i])))
        blocked |= np.abs(x - x[i]) < c.T14_B
    return taken      # list of (t, snr) for non-overlapping boxes, strongest first


res = {}
for ph in PHOTS:
    for meth, w in DETS:
        key = "%s|%s%.2f" % (ph, meth, w)
        try:
            t, f = ad.get(6, ph, window=w, method=meth)
        except Exception as e:
            res[key] = dict(error=str(e)); continue
        g = ad.BoxGrid(t, f)
        boxes = independent_boxes(g)
        N = len(boxes)
        # coverage-limited count of independent 1.6-h boxes = usable time / T14
        span_usable = np.sum(np.isfinite(g.depth)) * g.step
        ranks = []
        for tb in TB:
            # the box (from the greedy list) whose centre lies within T14/2 of the prediction; else the S/N at the prediction
            cand = [(k, b) for k, b in enumerate(boxes) if abs(b[0] - tb) < c.T14_B / 2]
            if cand:
                k, b = cand[0]
                ranks.append(dict(t_pred=tb, t_box=b[0], snr=b[1], rank=k + 1))
            else:
                d, s = g.lookup(np.array([tb]))
                ranks.append(dict(t_pred=tb, t_box=None, snr=float(d[0] / s[0]) if np.isfinite(s[0]) else None,
                                  rank=int(1 + sum(1 for b in boxes if b[1] > (d[0] / s[0]))) if np.isfinite(s[0]) else None))
        others = [b for b in boxes if all(abs(b[0] - tb) > c.T14_B for tb in TB)][:5]
        res[key] = dict(N_indep_greedy=N, usable_days=float(span_usable), N_indep_time=float(span_usable / c.T14_B), beta=float(g.beta),
                        B=ranks, top_others=others, top2_are_B=bool(sorted([r["rank"] or 999 for r in ranks]) == [1, 2]))
        print("%-26s N=%3d (time %.0f) beta %.2f  B: %s   next: %s" % (key, N, span_usable / c.T14_B, g.beta,
              [(r["rank"], r["snr"] and round(r["snr"], 1)) for r in ranks], [(round(o[0], 2), round(o[1], 1)) for o in others[:3]]), flush=True)

# rank-test p-value: probability that the two strongest of N independent boxes form an allowed pair.
# allowed pairs: first box anywhere in the scanned window (K positions) and the second at +P (spacing tolerance +-1 box -> 3)
def p_rank(N, K, spacings=3):
    return min(1.0, K * spacings / comb(N, 2))
Ks = {"W1_manuscript(+-0.004d)": 2 * 0.004 * 227 / c.T14_B, "W2(+-0.0045d)": 2 * 0.0045 * 227 / c.T14_B, "W3(+-0.0113d)": 2 * 0.0113 * 227 / c.T14_B}
Nref = res["pipe|biweight0.75"]["N_indep_time"]
prk = {k: p_rank(int(round(Nref)), v) for k, v in Ks.items()}
print("rank p-values (N=%d):" % Nref, prk)

# proximity to dumps / gaps; pixel-level: SAP_BKG and MOM_CENTR in-transit vs flanks (TESS-SPOC S6)
d = c.load_spoc_lc("s6_tessspoc")
t = d["t"]; good = np.isfinite(t) & ((d["q"] & c.SPOC_BAD) == 0) & np.isfinite(d["pdc"])
gaps = np.where(np.diff(t[good]) > 0.2)[0]
gap_edges = [(t[good][i], t[good][i + 1]) for i in gaps]
pix = []
for tb in TB:
    inn = good & (np.abs(t - tb) < c.T14_B / 2)
    fl = good & (np.abs(t - tb) > c.T14_B / 2) & (np.abs(t - tb) < c.T14_B / 2 + 3 * c.T14_B)
    ref = good & (np.abs(t - tb) > 0.5)
    row = dict(t=tb, dist_dump_h=[round((tb - x) * 24, 2) for x in c.DUMPS[6] if abs(tb - x) < 2], nearest_dump_h=float(np.min(np.abs(np.array(c.DUMPS[6]) - tb)) * 24),
               dist_gap_or_edge_d=float(min([abs(tb - a) for ab in gap_edges for a in ab] + [tb - t[good].min(), t[good].max() - tb])),
               n_in=int(inn.sum()), A_nearest_h=float(np.min(np.abs(c.T0_A + np.round((tb - c.T0_A) / c.P_A) * c.P_A - tb)) * 24))
    for col in ["bkg", "cx", "cy"]:
        x = d[col]
        dx = np.nanmean(x[inn]) - np.nanmean(x[fl])
        # empirical scatter of the same statistic at random epochs
        rr = []
        rng = np.random.default_rng(6)
        for _ in range(400):
            tt = rng.uniform(t[good].min() + 0.3, t[good].max() - 0.3)
            a = good & (np.abs(t - tt) < c.T14_B / 2); b = good & (np.abs(t - tt) > c.T14_B / 2) & (np.abs(t - tt) < c.T14_B / 2 + 3 * c.T14_B)
            if a.sum() >= 2 and b.sum() >= 4:
                rr.append(np.nanmean(x[a]) - np.nanmean(x[b]))
        rr = np.array(rr)
        row[col + "_delta"] = float(dx); row[col + "_z_emp"] = float(dx / c.robust_std(rr))
    pix.append(row)
    print("S6 event", tb, row, flush=True)
# own TESScut: background and centroid at the two events
p = np.load("out/ownphot_s6_ap12.npz")
tt, bk, cx, cy, qq = p["t"], p["bkg"], p["cx"], p["cy"], p["q"]
gd = np.isfinite(tt) & ((qq & c.SPOC_BAD) == 0)
for tb, row in zip(TB, pix):
    for nm, x in [("own_bkg", bk), ("own_cx", cx), ("own_cy", cy)]:
        inn = gd & (np.abs(tt - tb) < c.T14_B / 2); fl = gd & (np.abs(tt - tb) > c.T14_B / 2) & (np.abs(tt - tb) < c.T14_B / 2 + 3 * c.T14_B)
        rng = np.random.default_rng(7); rr = []
        for _ in range(400):
            t0 = rng.uniform(tt[gd].min() + 0.3, tt[gd].max() - 0.3)
            a = gd & (np.abs(tt - t0) < c.T14_B / 2); b = gd & (np.abs(tt - t0) > c.T14_B / 2) & (np.abs(tt - t0) < c.T14_B / 2 + 3 * c.T14_B)
            if a.sum() >= 2 and b.sum() >= 4: rr.append(np.nanmean(x[a]) - np.nanmean(x[b]))
        row[nm + "_z_emp"] = float((np.nanmean(x[inn]) - np.nanmean(x[fl])) / c.robust_std(np.array(rr)))
    print("own pixel", tb, {k: round(v, 2) for k, v in row.items() if k.startswith("own")})
json.dump(dict(per_variant=res, p_rank=prk, N_ref=Nref, K=Ks, pixel=pix, gap_edges=gap_edges), open("out/b05_s6_rank.json", "w"), indent=1, default=float)
