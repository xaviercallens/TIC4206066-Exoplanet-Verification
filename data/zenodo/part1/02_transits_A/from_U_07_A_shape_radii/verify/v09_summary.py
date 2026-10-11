"""Verifier: tabulate own joint fits -> out/summary_verify_fits.csv"""
import glob, json, csv, os
rows = []
for fn in sorted(glob.glob("out/fit_*.json")):
    d = json.load(open(fn))
    k, b, T, R = d["k"], d["b"], d["T14_h"], d["Rp"]
    rows.append(dict(tag=d["tag"], n_epochs=d["n_epochs"], n_over_tau=round(d["n_over_tau"], 1),
                     k=f"{k[2]:.5f} -{k[2]-k[1]:.5f}/+{k[3]-k[2]:.5f}", b=f"{b[2]:.2f} -{b[2]-b[1]:.2f}/+{b[3]-b[2]:.2f}",
                     T14_h=f"{T[2]:.3f} +-{(T[3]-T[1])/2:.3f}", Rp=f"{R[2]:.3f} +-{(R[3]-R[1])/2:.3f}", rho_q977=round(d["rho_q977"], 2),
                     u1=round(d["u1"][2], 3), u2=round(d["u2"][2], 3)))
with open("out/summary_verify_fits.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); [w.writerow(r) for r in rows]
for r in rows:
    print(r)
