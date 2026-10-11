#!/bin/zsh
# waits for queue B, then r07 (RUWE robustness, short) and paired Teff 4012 runs (seeds 1000-1004)
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/referee/resolve_R5"
until grep -q Q_B_DONE logs/q_B.log; do sleep 20; done
PY=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/trice_env/bin/python
$PY r07_ruwe_robustness.py > logs/r07.log 2>&1
for S in B A; do $PY r05_direct_runs.py default_teff4012 $S 5 2>&1 | grep -v "pkg_resources\|Warn"; done
echo Q_C_DONE
