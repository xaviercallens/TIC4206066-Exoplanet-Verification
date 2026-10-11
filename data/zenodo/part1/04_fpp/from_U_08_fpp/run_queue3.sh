#!/bin/zsh
# waits for queue1 or queue2 to finish (max 2 workers), then fresh MOLUSC (age 5 Gyr, RUWE + Gaia contrast) and TRICERATOPS MOLUSC-mode runs
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
H="/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066/08_fpp"
until grep -q QUEUE1_DONE "$H/logs/queue1.log" || grep -q QUEUE2_DONE "$H/logs/queue2.log"; do sleep 20; done
T=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/trice_env/bin/python
MR=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/a08_cache/molusc_run
date -u; cd /private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/a08_cache/MOLUSC/code
$T BinaryStarGUI.py -v --age 5 --ruwe --gaia $MR/fresh5 05h23m26.36s -- -01d19m38.2s 1000000 0.606 2>&1 | grep -v -i "warn\|interp\|Rect\|bisp\|Regular\|Linear\|Clough\|transition\|new_row\|^ *$"
date -u; shasum -a 256 $MR/fresh5_kept.csv
cd "$H"
$T f02_trice_run.py A_moluscFresh5 A fold_A_S32+S98_6min 10 --molusc $MR/fresh5_kept.csv --seed0 3100
$T f02_trice_run.py B_moluscFresh5 B fold_B_S32+S98_6min 10 --molusc $MR/fresh5_kept.csv --seed0 3100
MP=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/a08_cache/molusc_prev/molusc_age5_kept.csv
$T f02_trice_run.py A_moluscPrev5 A fold_A_S32+S98_6min 10 --molusc $MP --seed0 3000
$T f02_trice_run.py B_moluscPrev5 B fold_B_S32+S98_6min 10 --molusc $MP --seed0 3000
echo QUEUE3_DONE
