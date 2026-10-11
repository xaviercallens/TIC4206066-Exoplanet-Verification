#!/bin/zsh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/referee/resolve_R5"
PY=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/trice_env/bin/python
MODES=(gaia_ol0.38_bl0.2 ol0.38_bl0.2 ol0.38_all gaia_only gaia_ol0.45_bl0.2)
for m in $MODES; do $PY r05_direct_runs.py $m A 10 2>&1 | grep -v "pkg_resources\|Warn"; done
echo Q_A_DONE
