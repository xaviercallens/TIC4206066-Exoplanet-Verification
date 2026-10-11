"""Audit 9: RUWE/RV comparison sample and IR neighbours via VizieR TAP (Gaia archive TAP was unresponsive)."""
import hashlib, json, datetime, sys
from pathlib import Path
import pyvo
HERE = Path(__file__).resolve().parent; D = HERE / "data"
svc = pyvo.dal.TAPService("https://tapvizier.cds.unistra.fr/TAPVizieR/tap")
def run(name, q):
    t = svc.search(q, maxrec=200000).to_table()
    p = D / f"{name}.ecsv"; t.write(p, format="ascii.ecsv", overwrite=True)
    rec = dict(name=name, url="https://tapvizier.cds.unistra.fr/TAPVizieR/tap", query=q, file=p.name,
               sha256=hashlib.sha256(p.read_bytes()).hexdigest(), nrows=len(t),
               utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"))
    open(D / "PROVENANCE.jsonl", "a").write(json.dumps(rec) + "\n"); print(name, len(t))
RA, DEC = 80.85983, -1.32696
for what in sys.argv[1:]:
    if what == "ruwe":
        run("ruwe_rv_comparison_vizier", 'SELECT "Source", Gmag, "BP-RP", Plx, RUWE, epsi, sepsi, IPDfmp, amax, NSS, '
            '"e_RV", "o_RV", "RVamp", "RVchi2", "RVgof", GRVSmag FROM "I/355/gaiadr3" '
            'WHERE Gmag BETWEEN 10.2 AND 10.8 AND "BP-RP" BETWEEN 1.60 AND 1.90 AND Plx > 10')
    if what == "ir":
        for cat, nm in [("II/246/out", "2mass_10as"), ("II/328/allwise", "allwise_10as")]:
            run(nm, f'SELECT * FROM "{cat}" WHERE 1=CONTAINS(POINT(\'ICRS\', RAJ2000, DEJ2000), CIRCLE(\'ICRS\', {RA}, {DEC}, {10/3600}))')
