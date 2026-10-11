#!/bin/bash
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/audit_env/bin/python
until [ -f logs/fits2_done.txt ]; do sleep 5; done
q1() {
$PY v03_fit.py o-1_hw0.45_SYS_RM_FLUX mann claret 98:-1 15000 > logs/f16.log 2>&1
$PY v05_vshape.py o2_hw0.45_SYS_RM_FLUX none > logs/v1.log 2>&1
$PY v05_vshape.py o-1_hw0.45_SYS_RM_FLUX 98:-1 > logs/v3.log 2>&1
}
q2() {
$PY v03_fit.py o-1_hw0.45_SYS_RM_FLUX mann claret none 15000 > logs/f17.log 2>&1
$PY v07_secscan.py > logs/scan.log 2>&1
$PY v05_vshape.py o2_hw0.45_SYS_RM_FLUX 98:-1 > logs/v2.log 2>&1
}
q1 & q2 & wait
echo ALLDONE > logs/fits3_done.txt
