#!/bin/zsh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
PY=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/trice_env/bin/python
cd "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp"
$PY f02_trice_run.py A_base A fold_A_S32+S98_6min 20 --seed0 1000
$PY f02_trice_run.py B_base B fold_B_S32+S98_6min 20 --seed0 1000
$PY f02_trice_run.py A_tic A fold_A_S32+S98_6min 3 --M 0.63 --R 0.648133 --Teff 4010 --plx 28.0487 --seed0 2000
$PY f02_trice_run.py B_tic B fold_B_S32+S98_6min 3 --M 0.63 --R 0.648133 --Teff 4010 --plx 28.0487 --seed0 2000
$PY f02_trice_run.py A_teff4012 A fold_A_S32+S98_6min 3 --Teff 4012 --seed0 2100
$PY f02_trice_run.py B_teff4012 B fold_B_S32+S98_6min 3 --Teff 4012 --seed0 2100
$PY f02_trice_run.py A_10min A fold_A_S32+S98_10min 3 --exptime_min 10 --seed0 2200
$PY f02_trice_run.py B_10min B fold_B_S32+S98_10min 3 --exptime_min 10 --seed0 2200
echo QUEUE1_DONE
