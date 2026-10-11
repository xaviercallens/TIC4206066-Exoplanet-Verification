"""Write out_numeric_rows.csv: every manuscript number in scope vs recomputed value."""
import json, csv
S = json.load(open('out_summary.json'))
em = json.load(open(S['emcee_file']))
f = lambda v, d=7: f'{v:.{d}f}'
rows = []
def add(q, ms, rc, st, method, note=''):
    rows.append(dict(section='2 (3.18 d signal)', quantity=q, manuscript_value=ms, recomputed_value=rc, status=st, method=method, note=note))
P1 = S['P_route1_T1']; P3 = S['P_route1_T3']; Pl = S['P_lsq']; Pe = S['P_emcee']
add('P_A (d)', '3.1827856 +- 0.0000058',
    f"route1 T1 {f(P1[0])} +- {f(P1[1])}; T3 {f(P3[0])} +- {f(P3[1])}; LSQ {f(Pl[0])} +- {f(Pl[1])}; emcee {f(Pe[0])} +- {f(Pe[1])}; BLS {S['P_bls']:.6f}",
    'MATCH', 'own template timing + injection errors; global batman LSQ/emcee; BLS', 'all within 0.9 sigma')
T1 = S['T0_route1_T1']; T3 = S['T0_route1_T3']; Tl = S['T0_lsq']; Te = S['T0_emcee']
add('T0_A (BJD_TDB)', '2460998.7375 +- 0.0023',
    f"route1 T1 {2457000+T1[0]:.4f} +- {T1[1]:.4f}; T3 {2457000+T3[0]:.4f}; Huber {2457000+S['T0_route1_huber'][0]:.4f}; LSQ {2457000+Tl[0]:.5f} +- {Tl[1]:.5f}; emcee {2457000+Te[0]:.5f} +{Te[2]:.4f}/-{Te[1]:.4f}",
    'MINOR_DIFF', 'same', 'recomputed T0 earlier by 1.2-5.3 min (<=1.5 sigma); method dependent')
json.dump(rows, open('/dev/null', 'w'))
with open('out_numeric_rows_core.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); [w.writerow(r) for r in rows]
print(rows)
