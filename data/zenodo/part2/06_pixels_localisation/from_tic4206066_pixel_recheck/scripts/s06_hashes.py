"""SHA256 of every file in the package (except SHA256SUMS itself)."""
import os
from common import PKG, sha256
lines = []
for dp, dn, fn in os.walk(PKG):
    dn[:] = sorted(d for d in dn if d != "__pycache__")
    for f in sorted(fn):
        rp = os.path.relpath(os.path.join(dp, f), PKG)
        if rp == "SHA256SUMS" or f == ".DS_Store": continue
        lines.append(f"{sha256(os.path.join(dp, f))}  {rp}")
open(os.path.join(PKG, "SHA256SUMS"), "w").write("\n".join(lines) + "\n")
print(len(lines), "files hashed")
