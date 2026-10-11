"""Figure 1 of the PROPOSED REVISION, VARIANT 2 (signal B reduced to a short paragraph), TIC 4206066 RNAAS note.

Three panels, adapted from make_figure_rev.py (the 5-panel variant-1 figure, archived in ../variant1_fullB/):
(a) = old (a): 3.18 d per-sector folds with the fitted template, logic of the package's a6_figure.py, but folded at the
      revised 19-transit ephemeris of the text (02_timing_A/verify/v4_results.json, drop_onesided_573_and_m1_R60)
      (a1_timing_depths imported read-only; bytecode writing disabled so nothing is written into the package).
(b) = old (d) exactly as last built: referee resolver R3 spec (referee/resolve_R3/out/r4_panel_d_spec.json): S32
      (POS_CORR-decorrelated) and S98 (TESScut) 3.18 d offsets and the S98 11.13 d offset from audit 03 (11.13 d ring centred on that S98 point), error per axis =
      max(injection SD, jackknife SE, robust bootstrap SD); 3" rings at the combined positions; zoom inset of the central
      22"; target '+', 21" pixel circle, Gaia DR3 stars able to produce the 3.18 d depth (audit-03 G4 census).
(c) = old (e) for the 3.18 d signal only: TRICERATOPS FPP at five x positions (fold 1, fold 2, fold 3, TIC params,
      M_c <= 0.38), same data, sources and run counts; log y (1e-3 to 2); lines 0.015 (dashed) and 0.5 (dotted).
Layout: (a) tall on the left; (b) top right, (c) bottom right; 7.1 x 4.35 in.
Writes ../figure.pdf, ../figure.png (dpi 200) and ./figure_values_v2.json.
"""
import sys
sys.dont_write_bytecode = True          # never write __pycache__ into the read-only package
import csv, datetime, hashlib, json, os, platform
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, Ellipse
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple

HERE = os.path.dirname(os.path.abspath(__file__))
REV = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/proposed_revision"   # FINAL_RECONCILIATION: inputs read from the original package (read-only)
SC_OUT = "/Users/paulsmith/space <3/tess-jev/addenda/tic4206066_FINAL_RECONCILIATION/submission_candidate"   # FINAL_RECONCILIATION: figure output
# FINAL_RECONCILIATION: per-sector template depths drawn in panel (a) = FINAL_FACTS.csv A.depth_S6 / A.depth_S32 / A.depth_S98 (ppm);
# the previous build drew the package values 519.3 / 484.0 / 574.3, the last including the unmeasurable S98 segment-edge epoch -1 (superseded).
FINAL_DEPTH_PPM = {6: 521.0, 32: 476.0, 98: 496.0}
FINAL_DEPTH_FACT = {6: "A.depth_S6", 32: "A.depth_S32", 98: "A.depth_S98"}
AUD = os.path.dirname(REV)                                    # opus55_final_audit_tic4206066/
ADD = os.path.dirname(AUD)                                    # addenda/
PKG = os.path.join(ADD, "rnaas_tic4206066_pubprep")
PKG_AN = os.path.join(PKG, "analysis"); PKG_OUT = os.path.join(PKG_AN, "out")
PUB = os.path.join(ADD, "tic4206066_PUBLICATION_AUDIT")
PIX = os.path.join(AUD, "03_pixels"); FPP8 = os.path.join(AUD, "08_fpp")
sys.path.insert(0, PKG_AN)
import a1_timing_depths as a1          # noqa: E402  (read-only import)

SOURCES = {}
def src(path, what):
    """register an input file (relative to addenda/) with its sha256"""
    rel = os.path.relpath(path, ADD)
    if rel not in SOURCES:
        SOURCES[rel] = dict(sha256=hashlib.sha256(open(path, "rb").read()).hexdigest(), used_for=[])
    SOURCES[rel]["used_for"].append(what)
    return rel

def jload(path, what):
    src(path, what); return json.load(open(path))

A1_P = os.path.join(PKG_OUT, "a1_results.json")
A1 = jload(A1_P, "panel a ephemeris and template depths")
ephA_pkg = A1["signals"]["A"]["A__pipe"]["ephemeris"]["e_dt_inj_min"]          # package 20-transit ephemeris (variant 1 fold; kept for reference)
# panel (a) fold: the revised 19-transit ephemeris quoted in the text (audit 02 verifier, key drop_onesided_573_and_m1_R60)
V4_P = os.path.join(AUD, "02_timing_A", "verify", "v4_results.json"); V4 = jload(V4_P, "panel a fold ephemeris (19 transits)")
E19 = V4["ephem"]["drop_onesided_573_and_m1_R60"]
EPH = {"A": dict(P=E19["P"], T0=E19["T0"] - 2457000.0, T14=A1["signals"]["A"]["A__pipe"]["template_T14_h"] / 24)}
W, H = 7.1, 4.35
vals = dict(meta=dict(script=os.path.relpath(os.path.abspath(__file__), ADD), created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
                      python=platform.python_version(), numpy=np.__version__, matplotlib=matplotlib.__version__, variant="2 (three panels)",
                      based_on="figure_src/make_figure_rev.py (variant 1; archived in variant1_fullB/) <- rnaas_tic4206066_pubprep/analysis/a6_figure.py",
                      panel_map={"a": "old a", "b": "old d", "c": "old e (3.18 d only)"}, size_in=[W, H], png_dpi=200),
            ephemerides_used=dict(values=EPH, source=src(V4_P, "panel a fold"), key="ephem.drop_onesided_573_and_m1_R60 (T0 converted BJD_TDB -> BTJD)",
                                  n_transits=E19["n"], eP=E19["eP"], eT0_d=E19["eT0"], chi2=E19["chi2"], dof=E19["dof"],
                                  template_T14_source=f"{src(A1_P, 'template T14')}: signals.A.A__pipe.template_T14_h",
                                  package_ephemeris_not_used=dict(P=ephA_pkg["P"], T0_BTJD=ephA_pkg["T0"], source="rnaas_tic4206066_pubprep/analysis/out/a1_results.json"),
                                  note="FINAL_RECONCILIATION: template depths per sector = FINAL_FACTS A.depth_S6/S32/S98; S98 epoch -1 omitted from the fold; 11.13 d masking unchanged from the package"))

plt.rcParams.update({"font.size": 7.5, "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True, "font.family": "serif", "mathtext.fontset": "dejavuserif", "pdf.fonttype": 42})
fig = plt.figure(figsize=(W, H))
def box(x0, y0, w, h): return [x0 / W, y0 / H, w / W, h / H]
HB = 2.2; WB = HB * 170 / 118                         # (b): equal aspect, x range 170", y range 118"
XB0 = W - 0.05 - WB; YB0 = H - 0.08 - HB
axa = fig.add_axes(box(0.55, 0.42, XB0 - 0.62 - 0.55, H - 0.08 - 0.42))
axb = fig.add_axes(box(XB0, YB0, WB, HB))
axc = fig.add_axes(box(XB0, 0.36, WB, 1.2))
COL = {6: "#1f4e9c", 32: "#b2182b", 98: "#1b7837"}; OFF = {6: 0, 32: -1500, 98: -3000}; LAB = {6: "S6, TESS-SPOC 30 min", 32: "S32, SPOC 2 min", 98: "S98, QLP 200 s"}

# ------------------------------------------------------------------ (a): as in a6_figure.py (3.18 d folds)
for k, fn in (("S6", "hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits"), ("S32", "tess2020324010417-s0032-0000000004206066-0200-s_lc.fits"),
              ("S98", "hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits")):
    src(os.path.join(a1.LC, fn), f"panel a light curve {k}")
lcs = {s: a1.detrend(a1.load_lc(k)) for s, k in ((6, "S6"), (32, "S32"), (98, "S98"))}
def binned(x, y, w, sig):
    edges = np.arange(-6, 6 + 1e-9, w); c = 0.5 * (edges[1:] + edges[:-1]); m, e = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (x >= lo) & (x < hi); m.append(np.mean(y[s]) if s.sum() >= 2 else np.nan); e.append(sig / np.sqrt(s.sum()) if s.sum() >= 2 else np.nan)
    return c, np.array(m), np.array(e)
def r(a, n=2): return [None if not np.isfinite(v) else round(float(v), n) for v in np.asarray(a, float)]
e = EPH["A"]; tpl = a1.Template("A", "A", e["P"]); pv = {}
for s in (6, 32, 98):
    L = lcs[s]; t, fd = L["t"], L["fd"]; keep = ~a1.in_tr(t, "B", 1.0)
    if s == 98:   # FINAL_RECONCILIATION: omit the S98 segment-edge epoch -1 (BTJD 3995.555), as in the FINAL S98 depth (A.depth_S98)
        keep &= np.round((t - e["T0"]) / e["P"]) != -1
    x = (((t - e["T0"]) / e["P"] + 0.5) % 1 - 0.5) * e["P"] * 24; m = keep & (np.abs(x) < 6); y = (fd[m] - 1) * 1e6; sg = 1.4826 * np.median(np.abs(y - np.median(y)))
    bw = 0.5 if s == 6 else 1 / 3
    c, mb, eb = binned(x[m], y, bw, sg); axa.errorbar(c, mb + OFF[s], eb, fmt="o", ms=2.2, lw=0.6, color=COL[s], capsize=0)
    D = FINAL_DEPTH_PPM[s] * 1e-6   # FINAL_RECONCILIATION (was: A1 ... per_sector[s].e_D_inj[0])
    tt = np.linspace(-0.25, 0.25, 600); axa.plot(tt * 24, -D * tpl.unit(tt, L["cad"]) * 1e6 + OFF[s], color="k", lw=0.8)
    axa.text(5.5, OFF[s] + 330, LAB[s], ha="right", va="bottom", fontsize=6.3, color=COL[s])
    pv[str(s)] = dict(label=LAB[s], offset_ppm=OFF[s], bin_width_h=bw, template_depth_ppm=round(D * 1e6, 1), point_scatter_ppm=round(float(sg), 1),
                      bin_centre_h=r(c, 4), bin_mean_ppm_without_offset=r(mb, 1), bin_err_ppm=r(eb, 1),
                      depth_source=f"FINAL_FACTS.csv {FINAL_DEPTH_FACT[s]} (package value {A1['signals']['A']['A__pipe']['per_sector'][str(s)]['e_D_inj'][0]*1e6:.1f} ppm not used)",
                      epoch_minus1_omitted=(s == 98))
axa.set_xlim(-6, 6); axa.set_ylim(-4600, 1150); axa.set_xlabel("hours from mid-transit"); axa.set_ylabel("relative flux (ppm) + offset")
axa.text(0.03, 0.985, f"(a)  P = {e['P']:.2f} d", transform=axa.transAxes, ha="left", va="top", fontweight="bold")
vals["panel_a"] = dict(signal="A", P_d=e["P"], T0_BTJD=e["T0"], other_signal_removed_within_T14=1.0, template="A shape of a1_timing_depths.SHAPE", per_sector=pv)

# ------------------------------------------------------------------ (b) [old panel d]: localisation (audit 03 values, referee resolver R3 selection and errors)
# Plotted per-sector points and their errors come from referee/resolve_R3/out/r4_panel_d_spec.json ("plot"): positions from audit 03,
# error per axis = max(injection SD, jackknife SE, robust bootstrap SD). Sector 6 (both signals) and the Sector 32 11.13 d point are
# not plotted (spec "do_not_plot"). Rings: radius 3" (the text's rounded +-3") centred on the combined positions of the spec.
SPEC_P = os.path.join(AUD, "referee", "resolve_R3", "out", "r4_panel_d_spec.json"); SPEC = jload(SPEC_P, "panel b points, errors, ring centres")
SUM_P = os.path.join(PIX, "summary.json"); PS = jload(SUM_P, "panel b cross-check of the spec positions")
NB_P = os.path.join(PIX, "neighbours_vs_fitted_position.csv"); src(NB_P, "panel b capable Gaia DR3 stars (G4 flux basis, three light-curve apertures)")
axb.set_aspect("equal"); axb.add_patch(plt.Circle((0, 0), 21, fill=False, ls=":", lw=0.6, color="k")); axb.text(0, 22.5, "21$''$ (1 pixel)", ha="center", va="bottom", fontsize=6)
MK = {32: "s", 98: "^"}
SETKEY = {("A", 32): "S32_tpf_poscorr_A", ("A", 98): "S98_clean29_A", ("B", 98): "S98_clean29_B"}
PTS = [dict(sig="A" if q["signal"] == "3.18 d" else "B", **q) for q in SPEC["plot"]]
RING = {"A": dict(mean=SPEC["text_numbers"]["A"]["exact"]["mean"], key="text_numbers.A.exact (= ellipses.A_S32_S98): S32 decorrelated + S98, event-scatter errors"),
        "B": dict(mean=SPEC["ellipses"]["B_S98"]["mean"], key="ellipses.B_S98: the plotted S98 point itself (S32 11.13 d point not used)")}
R_DRAW = 3.0                                      # drawn 1-sigma radius: the text's rounded +-3"
pb_ = dict(points=[], combined={}, neighbours=[], pixel_circle_radius_arcsec=21.0, target_marker=[0.0, 0.0],
           not_plotted=SPEC["do_not_plot"], spec_source=src(SPEC_P, "b"))
def draw_points(ax, ms=4):
    for q in PTS:
        s = q["sector"]; col = COL[s]; fill = q["sig"] == "A"
        ax.errorbar(q["off"][0], q["off"][1], xerr=q["err"][0], yerr=q["err"][1], fmt=MK[s], ms=ms, mfc=col if fill else "white", mec=col, color=col,
                    lw=0.7, capsize=0, zorder=3)
def draw_comb(ax, sig, lw=1.0):
    m_ = RING[sig]["mean"]
    ax.add_patch(Ellipse((m_[0], m_[1]), 2 * R_DRAW, 2 * R_DRAW, fill=False, lw=lw, ls="-" if sig == "A" else "--", color="k", zorder=4))
# layout: the full-width slot is used (x range widened to the west) so that a zoom of the central 22" fits inside the panel
XL, YL = (62, -108), (-66, 52)                       # east left; x range fills the gridspec slot at equal aspect
ZX, ZY = (11.0, -11.0), (-11.0, 11.0)                # zoom window: central 22" (caption: "inset: central 22''")
axz = axb.inset_axes([0.505, 0.075, 0.40, 0.58]); axz.set_anchor("W")   # lower right of panel b, clear of all crosses and axis labels
axz.set_aspect("equal"); axz.set_xlim(*ZX); axz.set_ylim(*ZY)
axz.tick_params(labelsize=5.8, length=2, pad=1.5, labelleft=False, labelright=True); axz.set_xticks([10, 5, 0, -5, -10]); axz.set_yticks([-10, -5, 0, 5, 10])
axz.plot(0, 0, "+", color="k", ms=6, mew=0.8, zorder=1)
draw_points(axz, ms=3.6)
for sig in ("A", "B"): draw_comb(axz, sig, lw=0.9)
kh = [Line2D([], [], marker="s", ls="none", mfc=COL[32], mec=COL[32], ms=3.6), Line2D([], [], marker="^", ls="none", mfc=COL[98], mec=COL[98], ms=3.6)]
axz.legend(kh, ["S32", "S98"], loc="lower center", ncol=2, frameon=False, fontsize=5.8, handletextpad=0.1, columnspacing=0.9, borderaxespad=0.2, handlelength=1.0)
axb.add_patch(plt.Rectangle((ZX[1], ZY[0]), ZX[0] - ZX[1], ZY[1] - ZY[0], fill=False, lw=0.5, color="k", zorder=1))
for yy in ZY:                                        # connectors: west corners of the zoom window -> east corners of the inset
    fig.add_artist(ConnectionPatch(xyA=(ZX[1], yy), coordsA=axb.transData, xyB=(ZX[0], yy), coordsB=axz.transData, lw=0.5, color="k"))
for q in PTS:
    k = SETKEY[(q["sig"], q["sector"])]; a = PS["sets"][k]
    pb_["points"].append(dict(signal=q["sig"], sector=q["sector"], colour=COL[q["sector"]], marker=MK[q["sector"]], filled=q["sig"] == "A", product=q["product"],
                              offset_E_N_arcsec=q["off"], sigma_E_N_arcsec=q["err"], error_rule="max(injection SD, jackknife SE, robust bootstrap SD) per axis (spec)",
                              injection_sd_audit03=a["inj_std"], n_events_audit03=a["n_events"], audit03_set=k,
                              max_abs_diff_spec_vs_audit03_offset=float(np.max(np.abs(np.array(q["off"]) - np.array(a["real_off"])))),
                              source=f"{src(SPEC_P, 'b')}: plot[{SPEC['plot'].index({kk: vv for kk, vv in q.items() if kk != 'sig'})}]"))
for sig in ("A", "B"):
    pb_["combined"][sig] = dict(mean_E_N_arcsec=RING[sig]["mean"], drawn_radius_arcsec=R_DRAW, drawn_style="solid" if sig == "A" else "dashed",
                                spec_key=RING[sig]["key"], text=SPEC["text_numbers"][sig]["value"],
                                exact=SPEC["text_numbers"][sig]["exact"] if sig == "A" else SPEC["text_numbers"]["B"]["exact_S98"],
                                alternative_not_drawn=None if sig == "A" else dict(key="text_numbers.B.exact_S32PC_S98 (S32 decorrelated + S98; drawn in the previous build)",
                                                                                     mean=SPEC["text_numbers"]["B"]["exact_S32PC_S98"]["mean"]),
                                source=f"{src(SPEC_P, 'b')}")
nb = []
with open(NB_P) as fh:
    for row in csv.DictReader(fh):
        if float(row["dreqA_min"]) <= 1.0:
            nb.append(dict(source_id=row["source_id"], TIC=int(float(row["TIC"])), Tmag=float(row["Tmag"]), sep_arcsec=float(row["sep"]),
                           dE=float(row["dE"]), dN=float(row["dN"]), required_eclipse_depth_min=float(row["dreqA_min"]),
                           mahalanobis_from_audit_position=float(row["z_A_S32tpfPC+S98clean29_b0"])))
inr = [n for n in nb if min(XL) < n["dE"] < max(XL) and min(YL) < n["dN"] < max(YL)]
for n in nb: n["plotted"] = n in inr
axb.plot([n["dE"] for n in inr], [n["dN"] for n in inr], "x", color="k", ms=4.5, mew=0.9); axb.plot(0, 0, "+", color="k", ms=8, mew=1.0, zorder=5)
pb_["neighbours"] = nb
pb_["neighbour_rule"] = ("audit-03 G4 census, delta = 500 ppm: required eclipse depth <= 1 in at least one of the three light-curve apertures used "
                         "(S6 TESS-SPOC, S32 SPOC, S98 QLP 3x3); positions = mean over those products' epochs; column dreqA_min")
pb_["axis_limits"] = dict(x=list(XL), y=list(YL), east="left", zoom_inset_x=list(ZX), zoom_inset_y=list(ZY),
                          note="offsets, error bars and rings are drawn in the zoom inset (central 22\"); the full field shows the target, 21\" circle, "
                               "capable stars and the zoom window")
axb.set_xlim(*XL); axb.set_ylim(*YL); axb.set_xlabel(r"$\Delta\alpha$ (arcsec, east left)"); axb.set_ylabel(r"$\Delta\delta$ (arcsec)")
hand = [Line2D([], [], color="k", lw=0.9), Line2D([], [], color="k", lw=0.9, ls="--"), Line2D([], [], ls="none"), Line2D([], [], marker="x", ls="none", color="k", ms=4.5)]
labs = ["3.18 d (filled)", "11.13 d (open)", "1$\\sigma$ rings: S32+S98 (3.18 d),\nS98 (11.13 d)", "Gaia DR3 stars"]   # line handles = ring style
axb.legend(hand, labs, loc="upper right", frameon=False, handletextpad=0.4, borderaxespad=0.2, borderpad=0.2, labelspacing=0.2, handlelength=2.0)
axb.text(0.03, 0.97, "(b)", transform=axb.transAxes, ha="left", va="top", fontweight="bold")
vals["panel_b"] = pb_


# ------------------------------------------------------------------ (c) [old panel e, 3.18 d only]: TRICERATOPS FPP
def stats(arr):
    a = np.asarray(arr, float); return dict(mean=float(a.mean()), sd=float(a.std(ddof=1)), min=float(a.min()), max=float(a.max()), n_runs=int(len(a)))
def own(sig):
    p = os.path.join(PKG_OUT, "fpp", f"fpp_ownfold_{sig}.json"); j = jload(p, f"panel e fold 1 {sig}")
    return [x["FPP"] for x in j["runs"]], src(p, "e"), "package's own S32+S98 fold (a4_fpp_ownfold.py), stellar M 0.606 R 0.623 Teff 4030"
def pub(tag):
    p = os.path.join(PUB, "out", "fpp", "runs", tag, "summary.json"); j = jload(p, f"panel e {tag}")
    return [x["FPP"] for x in j["runs"]], src(p, "e"), j["log"]
def a08(tag):
    p = os.path.join(FPP8, "runs", tag, "runs.json"); j = jload(p, f"panel e {tag}")
    return [x["FPP"] for x in j["runs"]], src(p, "e"), j
def v08(tag):
    p = os.path.join(FPP8, "verify", "runs_v", tag + ".json"); j = jload(p, f"panel e {tag}")
    return [x["FPP"] for x in j["runs"]], src(p, "e"), j
XS = ["fold 1", "fold 2", "fold 3", "TIC params", r"$M_{\rm c}\leq0.38$" + "\n(post hoc)"]   # FINAL_RECONCILIATION pass-1 fix R1-05: label the post hoc cut in the panel
XKEY = ["fold1", "fold2", "fold3", "tic", "cap038"]
pc_ = dict(signal="A (3.18 d) only", x_labels=XS, markers="mean; thick bar = SD (ddof=1) of the runs; thin bar = run-to-run min-max",
           lines=dict(dashed=0.015, dotted=0.5), ylim=[1e-3, 2], points={})
f1, s1, d1 = own("A")
f2, s2, l2 = pub("A_base_S32+S98")
f3, s3, j3 = a08("A_base")
f4, s4, l4 = pub("A_tic_S32+S98")
f5, s5, j5 = v08("trunc0.38_A_v_A_S32+S98_5min")
pts = {"fold1": (f1, s1, d1 + " (6-min bins)"),
       "fold2": (f2, s2, "publication-audit baseline S32+S98 fold (tic4206066_PUBLICATION_AUDIT, 6-min bins), stellar M 0.606 R 0.623 Teff 4030"),
       "fold3": (f3, s3, f"this audit's fold (08_fpp {j3['fold']}, 6-min bins, empirical sigma), stellar M 0.606 R 0.623 Teff 4030"),
       "tic": (f4, s4, f"publication-audit fold 2 with TIC 8.2 stellar parameters M {l4['star_params']['mass']} R {l4['star_params']['rad']} Teff {l4['star_params']['Teff']}"),
       "cap038": (f5, s5, f"audit-08 VERIFIER direct TRICERATOPS runs with the companion-mass prior truncated at 0.38 Msun ({j5['mode']}), verifier's own fold "
                          f"{j5['fold']} (5-min bins, Teff {j5['Teff']})")}
for i, k in enumerate(XKEY):
    arr, sp, desc = pts[k]; st = stats(arr)
    axc.vlines(i, st["min"], st["max"], color="k", lw=0.6)
    axc.errorbar(i, st["mean"], st["sd"], fmt="o", ms=3.5, color="k", mfc="k", lw=1.2, capsize=0)
    pc_["points"][k] = dict(x=i, label=XS[i], **st, draws_per_run=1000000, description=desc, source=sp, per_run_FPP=[float(v) for v in arr])
axc.set_yscale("log"); axc.set_ylim(1e-3, 2); axc.axhline(0.015, color="k", lw=0.6, ls="--"); axc.text(2.5, 0.0135, "0.015 (validation)", ha="center", va="top", fontsize=6.3)   # pass-1 fix R1-05
axc.axhline(0.5, color="k", lw=0.6, ls=":"); axc.text(2.5, 0.45, "0.5 (likely planet)", ha="center", va="top", fontsize=6.3)   # pass-2 fix P2-R1-03
axc.set_xticks(range(len(XS))); axc.set_xticklabels(XS, fontsize=6.5); axc.set_xlim(-0.5, len(XS) - 0.5); axc.set_ylabel("TRICERATOPS FPP")
axc.text(0.015, 0.96, "(c)  3.18 d", transform=axc.transAxes, ha="left", va="top", fontweight="bold")
vals["panel_c"] = pc_

# ------------------------------------------------------------------ plotted values vs the numbers manuscript.tex (variant 2) quotes
cA, cB = vals["panel_b"]["combined"]["A"], vals["panel_b"]["combined"]["B"]
for n in nb:
    n["sigma_distance_drawn_3arcsec"] = float(np.hypot(n["dE"] - cA["mean_E_N_arcsec"][0], n["dN"] - cA["mean_E_N_arcsec"][1]) / R_DRAW)
P_ = pc_["points"]; dflt = ("fold1", "fold2", "fold3", "tic")
vals["manuscript_consistency"] = dict(
    note="numbers quoted in FINAL_RECONCILIATION/submission_candidate/manuscript.tex vs the plotted values (pass-1 fix R2-09: texts refer to the submission candidate)",
    a_fold_ephemeris=dict(plotted=dict(P=e["P"], T0_BTJD=e["T0"]), text=dict(P="3.182785 +- 0.000006", T0_BJD="2460998.736 +- 0.003 (19 transits)"),
                          P_rounds_to_text=bool(round(e["P"], 6) == 3.182785), T0_rounds_to_text=bool(round(e["T0"] + 2457000, 3) == 2460998.736),
                          eP_eT0=[E19["eP"], E19["eT0"]], dT0_vs_package_min=float((e["T0"] - ephA_pkg["T0"]) * 1440)),
    b_A_combined=dict(text="(+5 +- 3, 0 +- 3)''", plotted_centre=cA["mean_E_N_arcsec"], drawn_sigma=R_DRAW, exact_sigma=cA["exact"]["sig"],
                      sigma_distance_from_target_drawn=float(np.hypot(*cA["mean_E_N_arcsec"]) / R_DRAW), exact_p_vs_target=cA["exact"]["p"]),
    b_B_combined=dict(text="(+1 +- 3, -1 +- 3)''", plotted_centre=cB["mean_E_N_arcsec"], drawn_sigma=R_DRAW, exact_sigma=cB["exact"]["sig"]),
    b_capable_stars=dict(text="\\gtrsim43'' (>10 sigma)", nearest_sep=min(n["sep_arcsec"] for n in nb),
                         min_sigma_distance_drawn_3arcsec=min(n["sigma_distance_drawn_3arcsec"] for n in nb)),
    c_default_means=dict(text="0.03-0.04 (4-20 runs, 5-10 min bins)", plotted=[round(P_[k]["mean"], 4) for k in dflt],
                         runs=[P_[k]["n_runs"] for k in dflt], run_range=[min(P_[k]["min"] for k in dflt), max(P_[k]["max"] for k in dflt)]),
    c_capped=dict(text="excluding companions > 0.38 Msun (post hoc prior cut): ~0.002-0.007", plotted=[P_["cap038"]["mean"], P_["cap038"]["sd"], P_["cap038"]["min"], P_["cap038"]["max"]],
                  runs=P_["cap038"]["n_runs"]))

fig.savefig(os.path.join(SC_OUT, "figure.pdf")); fig.savefig(os.path.join(SC_OUT, "figure.png"), dpi=200)
vals["sources"] = SOURCES
json.dump(vals, open(os.path.join(HERE, "figure_values_final.json"), "w"), indent=1, default=float)
print(json.dumps(vals["manuscript_consistency"], default=float))
print("figure written")
