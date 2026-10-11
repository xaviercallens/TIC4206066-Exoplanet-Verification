#!/bin/zsh
# signal A, reduced run counts (machine shared with other agents; ~80 s per run)
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
cd "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/referee/resolve_R5"
PY=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/trice_env/bin/python
for spec in "gaia_ol0.38_bl0.2 8" "ol0.38_bl0.2 5" "gaia_only 5" "ol0.38_all 5"; do
  set -- ${=spec}
  $PY r05_direct_runs.py $1 A $2 2>&1 | grep -v "pkg_resources\|Warn"
done
echo Q_A_DONE
