"""Where did the pool pass-2 B ephemeris (P 11.132339, T0 4003.43952) fall in S6/S32, relative to B (adopted) and A transits?"""
import json, numpy as np
Pp, Tp = 11.132339360439559, 4003.439520363492
PB, TB = 11.13274, 4000.865
PA, TA = 3.182785, 2460998.736 - 2457000
out = []
for lo, hi, s in ((1468.0, 1491.0, 6), (2174.0, 2201.0, 32)):
    for n in range(int((lo - Tp) / Pp) - 1, int((hi - Tp) / Pp) + 2):
        t = Tp + n * Pp
        if not (lo <= t <= hi): continue
        dA = ((t - TA) / PA + 0.5) % 1 - 0.5
        dB = ((t - TB) / PB + 0.5) % 1 - 0.5
        out.append(dict(sector=s, n=n, parity="even" if n % 2 == 0 else "odd", t=round(t, 4), dt_to_A_h=round(dA * PA * 24, 2), dt_to_B_h=round(dB * PB * 24, 2)))
json.dump(out, open("/Users/paulsmith/space <3/tess-jev/addenda/tic4206066_FINAL_RECONCILIATION/_work/T5_signal_B/verify/v_pool_epoch.json", "w"), indent=1)
for o in out: print(o)
