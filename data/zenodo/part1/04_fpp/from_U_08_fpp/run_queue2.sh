#!/bin/zsh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
PY=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/trice_env/bin/python
M=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/a08_cache/molusc_prev/molusc_age5_kept.csv
cd "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp"
$PY f02_trice_run.py A_moluscPrev5 A fold_A_S32+S98_6min 10 --molusc $M --seed0 3000
$PY f02_trice_run.py B_moluscPrev5 B fold_B_S32+S98_6min 10 --molusc $M --seed0 3000
$PY f02_trice_run.py A_ap5x5 A fold_A_S32+S98_6min 3 --ap 5x5 --seed0 2300
$PY f02_trice_run.py B_ap5x5 B fold_B_S32+S98_6min 3 --ap 5x5 --seed0 2300
$PY f02_trice_run.py A_T0p A fold_A_S32+S98_6min_T0p 3 --seed0 2400
$PY f02_trice_run.py A_T0m A fold_A_S32+S98_6min_T0m 3 --seed0 2400
$PY f02_trice_run.py B_T0p B fold_B_S32+S98_6min_T0p 3 --seed0 2400
$PY f02_trice_run.py B_T0m B fold_B_S32+S98_6min_T0m 3 --seed0 2400
$PY f02_trice_run.py A_S32only A fold_A_S32_6min 3 --seed0 2500
$PY f02_trice_run.py A_S98only A fold_A_S98_6min 3 --seed0 2500
$PY f02_trice_run.py B_S98only B fold_B_S98_6min 3 --seed0 2500
$PY f02_trice_run.py A_sigformal A fold_A_S32+S98_6min 3 --sigma formal --seed0 2600
$PY f02_trice_run.py B_sigformal B fold_B_S32+S98_6min 3 --sigma formal --seed0 2600
echo QUEUE2_DONE
