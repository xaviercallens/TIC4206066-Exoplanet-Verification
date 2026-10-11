"""Shared helpers for the 00_data provisioning step (own code, written for this audit)."""
import hashlib, json, os, time, datetime
import requests

ROOT = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/00_data"
PROV = os.path.join(ROOT, "logs", "provenance.jsonl")
UA = {"User-Agent": "opus55-audit-tic4206066/1.0 (python-requests)"}

TIC = 4206066
GAIA_DR3 = 3220388198192519424
RA, DEC = 80.85983, -1.32696


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def log_prov(entry):
    os.makedirs(os.path.dirname(PROV), exist_ok=True)
    with open(PROV, "a") as f:
        f.write(json.dumps(entry, default=str) + "\n")


def http_get(url, dest, params=None, retries=4, timeout=600, note=None, stream=True):
    """GET url -> dest, with retries. Logs url, utc, sha256, bytes, selected HTTP headers."""
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    last = None
    for k in range(retries):
        t0 = utcnow()
        try:
            with requests.get(url, params=params, headers=UA, stream=stream, timeout=timeout) as r:
                final_url = r.url
                status = r.status_code
                if status != 200:
                    last = f"HTTP {status}"
                    time.sleep(3 * (k + 1))
                    continue
                tmp = dest + ".part"
                with open(tmp, "wb") as f:
                    for b in r.iter_content(1 << 20):
                        f.write(b)
                os.replace(tmp, dest)
                hdrs = {h: r.headers.get(h) for h in ("Content-Length", "Content-Type", "ETag", "Last-Modified", "Content-Disposition") if r.headers.get(h)}
            e = dict(relpath=os.path.relpath(dest, ROOT), url=final_url, utc_start=t0, utc=utcnow(),
                     http_status=status, bytes=os.path.getsize(dest), sha256=sha256_file(dest),
                     http_headers=hdrs, note=note)
            log_prov(e)
            return e
        except Exception as ex:  # network hiccup
            last = repr(ex)
            time.sleep(3 * (k + 1))
    e = dict(relpath=os.path.relpath(dest, ROOT), url=url, utc=utcnow(), error=last, note=note)
    log_prov(e)
    return e
