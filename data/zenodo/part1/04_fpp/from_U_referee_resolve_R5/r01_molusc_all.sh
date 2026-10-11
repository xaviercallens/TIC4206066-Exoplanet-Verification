#!/bin/zsh
# R5: fresh MOLUSC run (same settings as audit 08: age 5 Gyr, RUWE + Gaia contrast, 1e6 companions, M=0.606)
# with --all so that every generated companion and its Gaia/RUWE rejection flags are written.
# Used only to build the Gaia survival function S(m, a); output kept in the session scratchpad (large file).
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
T=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/trice_env/bin/python
C=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/r5_cache
date -u
cd /private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/a08_cache/MOLUSC/code
git log -1 --format=%H 2>/dev/null
$T BinaryStarGUI.py -v -a --age 5 --ruwe --gaia $C/r5 05h23m26.36s -- -01d19m38.2s 1000000 0.606 2>&1 | grep -v -i "warn\|interp\|Rect\|bisp\|Regular\|Linear\|Clough\|transition\|new_row\|^ *$"
date -u
shasum -a 256 $C/r5_kept.csv $C/r5_all.csv
echo R5_MOLUSC_DONE
