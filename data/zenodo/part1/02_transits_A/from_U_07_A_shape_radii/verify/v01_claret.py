"""Verifier: independent query of Claret (2017) TESS quadratic LDCs (VizieR J/A+A/600/A30), PHOENIX quasi-spherical
(tableab, Type q) and ATLAS (table25, Z=0, xi=2). Records URL, UTC time, sha256 of the saved table."""
import json, hashlib, datetime
import numpy as np
from astroquery.vizier import Vizier
out = {}
TEFF, LOGG = 4012.0, 4.438 + np.log10(0.606) - 2 * np.log10(0.623)
for tab, cons in [("J/A+A/600/A30/tableab", {"Teff": "4000..4100", "logg": "4.5..5.0"}),
                  ("J/A+A/600/A30/table25", {"Teff": "4000..4250", "logg": "4.5..5.0", "Z": "0", "xi": "2"})]:
    t = Vizier(columns=["**"], row_limit=-1).query_constraints(catalog=tab, **cons)[0]
    fn = "claret_" + tab.split("/")[-1] + ".ecsv"
    t.write(fn, format="ascii.ecsv", overwrite=True)
    print(tab); print(t)
    out[tab] = dict(file=fn, n=len(t), sha256=hashlib.sha256(open(fn, "rb").read()).hexdigest(),
                    utc=datetime.datetime.utcnow().isoformat() + "Z", constraints=cons)
json.dump(dict(provenance=out, teff=TEFF, logg=float(LOGG)), open("claret_provenance_verify.json", "w"), indent=1)
