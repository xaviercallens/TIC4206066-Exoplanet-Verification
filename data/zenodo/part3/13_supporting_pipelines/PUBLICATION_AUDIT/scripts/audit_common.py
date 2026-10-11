"""Shared helpers for the TIC 4206066 publication audit.

Everything written by this audit goes under AUD (addenda/tic4206066_PUBLICATION_AUDIT/).
Every network download is logged (URL, UTC time, HTTP status, bytes, sha256, Last-Modified/ETag)
to data/provenance.jsonl. Nothing outside AUD is ever written.
"""
import datetime as _dt
import hashlib
import json
import os
import time
import urllib.request

AUD = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PROJ = os.path.abspath(os.path.join(AUD, "..", ".."))           # tess-jev/
OLD_DEEP = os.path.join(PROJ, "addenda", "deep_followup_tic4206066")
OLD_COMP = os.path.join(PROJ, "addenda", "companion_check_tic4206066")
OLD_PRIO = os.path.join(PROJ, "addenda", "priority_recheck_tic4206066_SEALED")
OLD_RNAAS = os.path.join(PROJ, "addenda", "rnaas_tic4206066_draft")
PROV = os.path.join(AUD, "data", "provenance.jsonl")

TIC = 4206066
GAIA_DR3 = 3220388198192519424


def utcnow():
    return _dt.datetime.now(_dt.timezone.utc).isoformat()


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def log_prov(rec):
    os.makedirs(os.path.dirname(PROV), exist_ok=True)
    with open(PROV, "a") as f:
        f.write(json.dumps(rec, sort_keys=True) + "\n")


def fetch(url, dest, note="", retries=3, timeout=300, data=None, headers=None):
    """Download url -> dest (absolute or relative to AUD). Logs provenance. Returns record."""
    if not os.path.isabs(dest):
        dest = os.path.join(AUD, dest)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    last_err = None
    for attempt in range(retries):
        try:
            t0 = utcnow()
            req = urllib.request.Request(url, data=data, headers=headers or {"User-Agent": "tic4206066-audit/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                status = r.status
                hdr = {k: r.headers.get(k) for k in ("Last-Modified", "ETag", "Content-Length", "Content-Type")}
                tmp = dest + ".part"
                with open(tmp, "wb") as f:
                    while True:
                        b = r.read(1 << 20)
                        if not b:
                            break
                        f.write(b)
            os.replace(tmp, dest)
            rec = dict(url=url, dest=os.path.relpath(dest, AUD), utc_start=t0, utc_end=utcnow(), http_status=status,
                       bytes=os.path.getsize(dest), sha256=sha256_file(dest), headers=hdr, note=note,
                       post_data=(data.decode()[:2000] if isinstance(data, (bytes, bytearray)) else None))
            log_prov(rec)
            return rec
        except Exception as e:  # noqa
            last_err = repr(e)
            time.sleep(3 * (attempt + 1))
    rec = dict(url=url, dest=os.path.relpath(dest, AUD), utc=utcnow(), error=last_err, note=note)
    log_prov(rec)
    raise RuntimeError(f"fetch failed {url}: {last_err}")


def save_json(obj, rel):
    p = rel if os.path.isabs(rel) else os.path.join(AUD, rel)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as f:
        json.dump(obj, f, indent=1, sort_keys=False, default=_json_default)
    return p


def _json_default(o):
    import numpy as np
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)
