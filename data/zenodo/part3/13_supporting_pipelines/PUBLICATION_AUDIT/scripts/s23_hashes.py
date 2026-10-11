"""Step 23: outputs_sha256.json — sha256 of every file in the audit directory (deliverables, scripts, tests, outputs,
plots, logs, downloaded data), excluding only the hash file itself. Run last."""
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from audit_common import AUD, utcnow  # noqa

OUT = os.path.join(AUD, "outputs_sha256.json")
files = {}
total = 0
for root, dirs, names in os.walk(AUD):
    dirs[:] = sorted(d for d in dirs if d != "__pycache__")
    for n in sorted(names):
        p = os.path.join(root, n)
        if os.path.abspath(p) == OUT or n == ".DS_Store":
            continue
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        rel = os.path.relpath(p, AUD)
        files[rel] = dict(sha256=h.hexdigest(), bytes=os.path.getsize(p))
        total += os.path.getsize(p)
json.dump(dict(generated_utc=utcnow(), audit_dir="addenda/tic4206066_PUBLICATION_AUDIT", n_files=len(files), total_bytes=total,
               note="every file under the audit directory except this one; __pycache__ and .DS_Store skipped", files=files),
          open(OUT, "w"), indent=1)
print(f"{len(files)} files, {total / 1e6:.1f} MB -> outputs_sha256.json")
