#!/bin/zsh
# waits for queue B, then paired Teff 4012 vs 4030 runs (seeds 1000-1004, same as the package's adopted-fold runs)
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/referee/resolve_R5"
until grep -q Q_B_DONE logs/q_B.log; do sleep 20; done
PY=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/trice_env/bin/python
for S in A B; do $PY r05_direct_runs.py default_teff4012 $S 5 2>&1 | grep -v "pkg_resources\|Warn"; done
echo Q_C_DONE
