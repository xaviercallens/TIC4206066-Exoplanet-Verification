"""Fetch Claret (2017, A&A 600, A30) TESS quadratic limb-darkening coefficients from VizieR (J/A+A/600/A30)
near the target (Teff ~4012 K, logg ~4.63, [M/H]=0) and interpolate bilinearly in (Teff, logg).
Output: out/ldc_claret2017.json (+ raw rows csv)."""
import json, os, datetime
import numpy as np
from astroquery.vizier import Vizier
D = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(D, "out")
TEFF, LOGG = 4012.0, 4.631
v = Vizier(row_limit=-1, columns=["**"])
res = {}
for tab, cols in (("J/A+A/600/A30/table25", ("aLSM", "bLSM")), ("J/A+A/600/A30/tableab", ("aLSM", "bLSM"))):
    t = v.query_constraints(catalog=tab, Teff="3700..4300", logg="4.0..5.0", Z="0")[0]
    df = t.to_pandas()
    df.to_csv(os.path.join(OUT, f"ldc_raw_{tab.split('/')[-1]}.csv"), index=False)
    if "xi" in df.columns:
        df = df[np.isclose(df["xi"], 2.0)]
    if "Mod" in df.columns:
        mods = df["Mod"].unique().tolist()
    else:
        mods = []
    # PHOENIX table has 'Type'; keep all rows with the same method, prefer r-method if flagged
    out = {}
    for mod in (mods or [None]):
        dd = df if mod is None else df[df["Mod"] == mod]
        if "Type" in dd.columns:
            types = dd["Type"].unique().tolist()
        else:
            types = [None]
        for ty in types:
            d2 = dd if ty is None else dd[dd["Type"] == ty]
            if "L/HP" in d2.columns:
                d2 = d2[np.isclose(d2["L/HP"], d2["L/HP"].iloc[0])]
            Ts = np.sort(d2["Teff"].unique()); Gs = np.sort(d2["logg"].unique())
            try:
                T0 = Ts[Ts <= TEFF].max(); T1 = Ts[Ts >= TEFF].min(); G0 = Gs[Gs <= LOGG].max(); G1 = Gs[Gs >= LOGG].min()
            except ValueError:
                continue
            def val(T, G, c):
                r = d2[(d2["Teff"] == T) & (d2["logg"] == G)]
                return float(r[c].iloc[0])
            wt = 0 if T1 == T0 else (TEFF - T0) / (T1 - T0); wg = 0 if G1 == G0 else (LOGG - G0) / (G1 - G0)
            ab = []
            for c in cols:
                v00, v01, v10, v11 = val(T0, G0, c), val(T0, G1, c), val(T1, G0, c), val(T1, G1, c)
                ab.append((1 - wt) * (1 - wg) * v00 + (1 - wt) * wg * v01 + wt * (1 - wg) * v10 + wt * wg * v11)
            out[f"{mod}|{ty}"] = dict(u1=float(ab[0]), u2=float(ab[1]), grid=dict(T=[float(T0), float(T1)], logg=[float(G0), float(G1)]))
    res[tab] = out
res["_meta"] = dict(Teff=TEFF, logg=LOGG, MH=0.0, xi_ATLAS=2.0, source="VizieR J/A+A/600/A30 (Claret 2017)",
                    utc=datetime.datetime.now(datetime.UTC).isoformat())
json.dump(res, open(os.path.join(OUT, "ldc_claret2017.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
