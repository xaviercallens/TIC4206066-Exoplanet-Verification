"""Download Claret (2017, A&A 600, A30) TESS quadratic LD tables from VizieR and record provenance."""
import json, hashlib, datetime
from astroquery.vizier import Vizier
out = {}
for tab, cons in [("J/A+A/600/A30/tableab", {"Teff": "3800..4200", "logg": "4.0..5.0"}),
                  ("J/A+A/600/A30/table25", {"Teff": "3750..4250", "logg": "4.0..5.0", "Z": "-0.5..0.5"})]:
    v = Vizier(columns=["**"], row_limit=-1)
    t = v.query_constraints(catalog=tab, **cons)[0]
    fn = "claret/" + tab.split("/")[-1] + ".ecsv"
    t.write(fn, format="ascii.ecsv", overwrite=True)
    sha = hashlib.sha256(open(fn, "rb").read()).hexdigest()
    out[tab] = dict(file=fn, constraints=cons, n=len(t), sha256=sha, utc=datetime.datetime.utcnow().isoformat() + "Z",
                    url="https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=" + tab)
    print(tab, len(t))
    print(t)
json.dump(out, open("claret/claret_provenance.json", "w"), indent=1)
