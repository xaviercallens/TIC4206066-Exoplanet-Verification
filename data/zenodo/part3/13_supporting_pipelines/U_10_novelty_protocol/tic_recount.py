"""Audit 11: independent recount of the QLP Sector 98 stars with TIC 8.2 Tmag <= 11.5.
Inputs fetched fresh from public archives (not from the project's index stage):
  * QLP S98 target list  https://archive.stsci.edu/hlsps/qlp/target_lists/s0098.csv
  * TIC 8.2 Tmag for every listed TIC ID via MAST Catalogs (astroquery.mast.Catalogs.query_criteria, catalog="Tic")
Checkpointed per batch (tic_batches/*.csv) so it can resume; log in tic_recount.log; result tic_recount.json.
"""
import os, sys, json, time, hashlib, datetime
import numpy as np, requests

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw"); os.makedirs(RAW, exist_ok=True)
BD = os.path.join(HERE, "tic_batches"); os.makedirs(BD, exist_ok=True)
URL = "https://archive.stsci.edu/hlsps/qlp/target_lists/s0098.csv"
LOG = os.path.join(HERE, "queries.log")


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def qlog(source, url, st, nb, sha, result):
    with open(LOG, "a") as f:
        f.write("\t".join([now(), source, url, str(st), str(nb), sha, result]) + "\n")


fn = os.path.join(RAW, "qlp_s0098_target_list.csv")
if not os.path.exists(fn):
    r = requests.get(URL, headers={"User-Agent": "novelty-audit/1.0"}, timeout=600)
    open(fn, "wb").write(r.content)
    qlog("MAST QLP S98 target list", URL, r.status_code, len(r.content), hashlib.sha256(r.content).hexdigest(),
         f"lines={r.content.count(b'\\n')}")
ids = []
for line in open(fn):
    if line.startswith("#") or not line.strip():
        continue
    ids.append(int(line.split(",")[0]))
ids = np.array(sorted(set(ids)))
print("unique TIC in QLP S98 list:", len(ids), flush=True)

from astroquery.mast import Catalogs
B = 5000
nb = int(np.ceil(len(ids) / B))
for i in range(nb):
    out = os.path.join(BD, f"b{i:04d}.csv")
    if os.path.exists(out):
        continue
    chunk = [str(x) for x in ids[i * B:(i + 1) * B]]
    for att in range(5):
        try:
            t = Catalogs.query_criteria(catalog="Tic", ID=chunk)
            t = t[["ID", "Tmag"]]
            t.write(out + ".tmp", format="csv", overwrite=True)
            os.replace(out + ".tmp", out)
            break
        except Exception as e:
            print("batch", i, "attempt", att, repr(e)[:200], flush=True)
            time.sleep(20 * (att + 1))
    print(now(), "batch", i + 1, "/", nb, flush=True)

import pandas as pd
tabs = [pd.read_csv(os.path.join(BD, f"b{i:04d}.csv")) for i in range(nb) if os.path.exists(os.path.join(BD, f"b{i:04d}.csv"))]
t = pd.concat(tabs).drop_duplicates("ID")
res = dict(qlp_list_url=URL, qlp_list_sha256=hashlib.sha256(open(fn, "rb").read()).hexdigest(), n_unique_tic_in_list=int(len(ids)),
           n_tic_rows_returned=int(len(t)), n_missing_in_tic=int(len(set(ids) - set(t.ID.astype(int)))),
           n_tmag_nan=int(t.Tmag.isna().sum()), n_tmag_le_11p5=int((t.Tmag <= 11.5).sum()),
           target_4206066_Tmag=float(t.loc[t.ID == 4206066, "Tmag"].iloc[0]) if (t.ID == 4206066).any() else None,
           n_batches=nb, n_batches_done=len(tabs), finished_utc=now())
json.dump(res, open(os.path.join(HERE, "tic_recount.json"), "w"), indent=1)
qlog("MAST TIC 8.2 (astroquery Catalogs.query_criteria catalog=Tic, ID batches of 5000)", "MAST Mast.Catalogs.Filtered.Tic", 200, 0, "-",
     f"TIC rows {len(t)}; Tmag<=11.5: {res['n_tmag_le_11p5']}")
print(json.dumps(res, indent=1))
