"""Named analysis configurations for the signal-B depth audit (one consistent method + one-switch variants)."""
import os, json, copy
from b_run import *

PBASE = {6: "S6_tessspoc_pdc", 32: "S32_spoc2m_pdc", 98: "S98_qlp_sysrm"}

# fallback shapes (replaced by the s02 joint-fit best-fit shapes when those files exist)
_FALLBACK = {"prior": dict(k=0.0328, b=0.877, rho=2.51, u1=0.32, u2=0.31),
             "free": dict(k=0.0297, b=0.45, rho=12.0, u1=0.32, u2=0.31)}


def shape(name):
    fn = os.path.join(OUT, f"s02_joint_{name}_default.json")
    if name == "free" and os.path.exists(fn):
        # free-density posterior is a long b-rho ridge; use its central shape (posterior-median b and T14)
        d = json.load(open(fn))
        return shape_from_T14(d["b"]["med"], d["T14_h"]["med"], k=d["k"]["med"])
    if os.path.exists(fn):
        d = json.load(open(fn)); b = d["best"]
        return dict(k=b["k"], b=b["b"], rho=10 ** b["log10rho"], u1=b["u1"], u2=b["u2"])
    return dict(_FALLBACK[name])


_TPLS = {}


def template_for(cfg):
    nm = cfg.get("template", "prior")
    if nm not in _TPLS:
        _TPLS[nm] = Template(shape(nm))
    return _TPLS[nm]


BASE = dict(products=PBASE, spoc_mask=17087, s98_mask=MASK_S98_BASE, A_mode="fixed", A_depth_ppm=522.0,
            detrend=dict(method="biweight", window=0.75), poly=1, hw=0.3, dt="fixed", supersample=True, template="prior")


def var(**kw):
    c = copy.deepcopy(BASE)
    for k, v in kw.items():
        if k == "products":
            p = dict(c["products"]); p.update(v); c["products"] = p
        else:
            c[k] = v
    return c


CONFIGS = {
    "BASE": BASE,
    # manuscript-like choices (A cadences within +-T14_A dropped, per-event free timing, lightkurve default mask)
    "MSLIKE": var(A_mode="mask1.0", dt="event"),
    "MSLIKE_fixeph": var(A_mode="mask1.0"),
    # quality mask
    "KEEP512": var(spoc_mask=16575),
    # A treatment
    "A_MODEL": var(A_mode="model"),
    "A_MASK06": var(A_mode="mask0.6"),
    "A_MASK10": var(A_mode="mask1.0"),
    # detrending
    "BW050": var(detrend=dict(method="biweight", window=0.5)),
    "BW100": var(detrend=dict(method="biweight", window=1.0)),
    "BW150": var(detrend=dict(method="biweight", window=1.5)),
    "COSINE": var(detrend=dict(method="cosine", window=1.0, robust=True)),
    "RSPLINE": var(detrend=dict(method="rspline", window=0.75)),
    "HSPLINE": var(detrend=dict(method="hspline", window=0.75)),
    "LOCALPOLY2": var(detrend=dict(method="none"), poly=2, hw=0.3),
    "LOCALPOLY1": var(detrend=dict(method="none"), poly=1, hw=0.25),
    "POLY0": var(poly=0),
    # timing
    "DT_SECTOR": var(dt="sector"),
    "DT_EVENT": var(dt="event"),
    # template / cadence
    "FREE_TPL": var(template="free"),
    "NO_SUPERSAMPLE": var(supersample=False),
    # coverage / event selection in S98
    "NO_PARTIAL": var(drop_partial=True),
    "S98_KEEP_B30": var(s98_mask=QLP_SPACECRAFT | QLP_STRAY | QLP_B29 | 4096),      # partial event 4 becomes full
    "S98_WITH_EP2": var(s98_mask=QLP_SPACECRAFT | QLP_STRAY | QLP_B30 | 4096),       # adds epoch 2 (b29-only interval)
    # products / flux scale
    "S32_TESSSPOC10": var(products={32: "S32_tessspoc10_pdc"}),
    "S32_SAPCROWD": var(products={32: "S32_spoc2m_sap"}),
    "S6_SAPCROWD": var(products={6: "S6_tessspoc_sap"}),
    "S6_QLP": var(products={6: "S6_qlp_ksp"}),
    "S6_TGLC": var(products={6: "S6_tglc_aper"}),
    "S6_ELEANOR": var(products={6: "S6_eleanor_corr"}),
    "QLP_ALL": var(products={6: "S6_qlp_ksp", 32: "S32_qlp_ksp"}),
    "QLP_ALL_SAP": var(products={6: "S6_qlp_sap", 32: "S32_qlp_sap", 98: "S98_qlp_sap"}),
    "S98_DET": var(products={98: "S98_qlp_det"}),
    "S98_SAP": var(products={98: "S98_qlp_sap"}),
    "S98_OWN_r2.0": var(products={98: "own98_r2.0"}),
    # own TESScut photometry for all three sectors, several apertures
    "OWN_r1.5": var(products={6: "own6_r1.5", 32: "own32_r1.5", 98: "own98_r1.5"}),
    "OWN_r2.0": var(products={6: "own6_r2.0", 32: "own32_r2.0", 98: "own98_r2.0"}),
    "OWN_r2.5": var(products={6: "own6_r2.5", 32: "own32_r2.5", 98: "own98_r2.5"}),
    "OWN_r3.0": var(products={6: "own6_r3.0", 32: "own32_r3.0", 98: "own98_r3.0"}),
}
# own photometry with the 512 flag kept is identical (TESScut cubes carry no PDC outlier flags)
