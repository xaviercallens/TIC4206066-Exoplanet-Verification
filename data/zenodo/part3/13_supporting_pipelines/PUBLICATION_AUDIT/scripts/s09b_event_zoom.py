"""Step 9b: close-up diagnostics of S98 events in QLP-flagged intervals (bit 2^29 / 2^30 / straylight), where the frozen
mask removed data: A epochs 3, 8, 12, -1 and B epochs 2, 3, -1, 4. For each: QLP SAP/SYS flux with its quality bits,
own TESScut photometry (r = 1.5/2/2.5 px), local background level and slope, centroid; binned to 20 min."""
import json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from astropy.io import fits
sys.path.insert(0, os.path.dirname(__file__))
import lcio
from audit_common import AUD, save_json
EPH = json.load(open(os.path.join(AUD, "out", "ephem_current.json")))
q = fits.open(os.path.join(lcio.MAST, lcio.PRODUCTS["QLP_S98_SYS"][0]))[1].data
tq = np.asarray(q["TIME"], float); fq = np.asarray(q["SYS_RM_FLUX"], float); qq = np.asarray(q["QUALITY"]).astype(np.int64)
bq = np.asarray(q["SAP_BKG"], float)
z = np.load(os.path.join(AUD, "out", "phot", "tesscut_s98.npz"))
cases = [("A", n) for n in (-1, 0, 2, 3, 6, 8, 11, 12)] + [("B", n) for n in (-1, 0, 1, 2, 3, 4)]
fig, axs = plt.subplots(len(cases), 3, figsize=(18, 2.6 * len(cases)))
info = []
def binned(t, f, w=20 / 1440):
    e = np.arange(t.min(), t.max() + w, w); i = np.digitize(t, e)
    bt = np.array([t[i == k].mean() for k in np.unique(i)]); bf = np.array([np.nanmedian(f[i == k]) for k in np.unique(i)])
    return bt, bf
for row, (sig, n) in zip(axs, cases):
    e = EPH[sig]; tc = e["T0"] + n * e["P"]
    m = np.abs(tq - tc) < 0.6
    mo = np.abs(z["time"] - tc) < 0.6
    rec = dict(signal=sig, epoch=n, tc=tc)
    # QLP flags in the window
    bits = {b: int(((qq[m] >> b) & 1).sum()) for b in (11, 12, 29, 30) if ((qq[m] >> b) & 1).sum()}
    rec["qlp_bits_in_window"] = bits; rec["n_qlp"] = int(m.sum())
    rec["tesscut_straylight_frac"] = float(z["straylight"][mo].mean()) if mo.any() else None
    ax = row[0]
    if m.sum() > 5:
        good = np.isfinite(fq[m])
        ff = fq[m] / np.nanmedian(fq[m][np.abs(tq[m] - tc) > e["T14"]])
        ax.plot((tq[m] - tc) * 24, (ff - 1) * 1e6, ".", ms=1.5, color="0.6")
        bt, bf = binned(tq[m][good], ff[good]); ax.plot((bt - tc) * 24, (bf - 1) * 1e6, "o-", ms=3, color="C0", label="QLP SYS 20-min")
        flag = ((qq[m] & (1 << 29)) > 0)
        ax.plot((tq[m][flag] - tc) * 24, np.full(flag.sum(), 2500), "|", color="C3", ms=4, label="bit29")
    ax.axvspan(-e["T14"] * 12, e["T14"] * 12, color="C1", alpha=0.15); ax.set_ylim(-4000, 3000); ax.set_title(f"{sig} n={n} QLP (bits {bits})", fontsize=8)
    ax = row[1]
    for r_, c_ in ((1.5, "C2"), (2.0, "C0"), (2.5, "C4")):
        if mo.sum() > 5:
            f = z[f"flux_r{r_}"][mo]; f = f / np.median(f[np.abs(z["time"][mo] - tc) > e["T14"]])
            bt, bf = binned(z["time"][mo], f); ax.plot((bt - tc) * 24, (bf - 1) * 1e6, "o-", ms=2, color=c_, label=f"r={r_}")
    ax.axvspan(-e["T14"] * 12, e["T14"] * 12, color="C1", alpha=0.15); ax.set_ylim(-4000, 3000); ax.set_title(f"TESScut own (straylight frac {rec['tesscut_straylight_frac']})", fontsize=8); ax.legend(fontsize=6)
    ax = row[2]
    if mo.sum() > 5:
        ax.plot((z["time"][mo] - tc) * 24, z["bkg"][mo], ".", ms=1.5, color="k", label="TESScut bkg/pix")
        ax2 = ax.twinx(); ax2.plot((z["time"][mo] - tc) * 24, z["cx"][mo] - np.median(z["cx"][mo]), ".", ms=1, color="C3"); ax2.set_ylim(-0.05, 0.05)
        rec["bkg_median"] = float(np.median(z["bkg"][mo])); rec["bkg_range"] = float(np.ptp(z["bkg"][mo]))
    ax.axvspan(-e["T14"] * 12, e["T14"] * 12, color="C1", alpha=0.15); ax.set_title("background (black) and x-centroid (red)", fontsize=8)
    info.append(rec)
plt.tight_layout(); plt.savefig(os.path.join(AUD, "plots", "s98_flagged_events_zoom.png"), dpi=80); plt.close()
save_json(info, "out/s09b_event_zoom.json")
for r in info: print(r)
