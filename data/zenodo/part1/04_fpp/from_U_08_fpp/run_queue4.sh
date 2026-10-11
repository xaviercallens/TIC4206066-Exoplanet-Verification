#!/bin/zsh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
H="/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp"
until grep -q QUEUE2_DONE "$H/logs/queue2.log"; do sleep 20; done
T=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/trice_env/bin/python
MR=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/a08_cache/molusc_run
cd "$H"
until [ -f $MR/fresh5_kept.csv ]; do sleep 20; done
$T f02_trice_run.py A_moluscFresh5_N5e5 A fold_A_S32+S98_6min 3 --N 500000 --molusc $MR/fresh5_kept.csv --seed0 3200
$T f02_trice_run.py A_moluscFresh5_N2e6 A fold_A_S32+S98_6min 2 --N 2000000 --molusc $MR/fresh5_kept.csv --seed0 3300
echo QUEUE4_DONE
