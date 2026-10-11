"""Verifier: R5 paired seeds, Teff 4030 (default) vs 4012 (default_teff4012), same fold/target object (read-only)."""
import json
U = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/referee/resolve_R5/runs/"
a = {r["seed"]: r["FPP"] for r in json.load(open(U + "default_A_adoptedfold.json"))["runs"]}
b = {r["seed"]: r["FPP"] for r in json.load(open(U + "default_teff4012_A_adoptedfold.json"))["runs"]}
out = {str(s): dict(T4030=a[s], T4012=b[s], ratio=b[s] / a[s]) for s in sorted(set(a) & set(b))}
json.dump(out, open("v03_r5_teff_pairs.json", "w"), indent=1); print(out)
