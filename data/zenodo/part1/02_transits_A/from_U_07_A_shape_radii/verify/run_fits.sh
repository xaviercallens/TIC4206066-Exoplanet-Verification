#!/bin/bash
# two parallel queues of verifier fits
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/audit_env/bin/python
N=15000
q1() {
$PY v03_fit.py o2_hw0.45_SYS_RM_FLUX mann claret none $N > logs/f1.log 2>&1
$PY v03_fit.py o2_hw0.45_SYS_RM_FLUX uni claret 98:-1 $N > logs/f7.log 2>&1
$PY v03_fit.py o2_hw0.45_SYS_RM_FLUX mann free 98:-1 $N > logs/f11.log 2>&1
$PY v03_fit.py o2_hw0.45_DET_FLUX mann claret none $N > logs/f5.log 2>&1
$PY v03_fit.py o2_hw0.45_SYS_RM_FLUX uni claret none $N > logs/f9.log 2>&1
$PY v03_fit.py o2_hw0.45_DET_FLUX uni claret none $N > logs/f13.log 2>&1
$PY v03_fit.py o1_hw0.45_SYS_RM_FLUX mann claret none $N > logs/f14.log 2>&1
}
q2() {
$PY v03_fit.py o2_hw0.45_SYS_RM_FLUX mann claret 98:-1 $N > logs/f2.log 2>&1
$PY v03_fit.py o2_hw0.45_SYS_RM_FLUX logu claret 98:-1 $N > logs/f8.log 2>&1
$PY v03_fit.py o2_hw0.45_SYS_RM_FLUX uni free 98:-1 $N > logs/f12.log 2>&1
$PY v03_fit.py o2_hw0.45_DET_FLUX mann claret 98:-1 $N > logs/f6.log 2>&1
$PY v03_fit.py o2_hw0.45_SYS_RM_FLUX logu claret none $N > logs/f10.log 2>&1
$PY v03_fit.py o2_hw0.45_DET_FLUX logu claret none $N > logs/f15.log 2>&1
$PY v03_fit.py o1_hw0.45_SYS_RM_FLUX mann claret 98:-1 $N > logs/f4.log 2>&1
$PY v03_fit.py o2_hw0.45_SYS_RM_FLUX mann claret 98:-1,6:-795 $N > logs/f3.log 2>&1
}
q1 & q2 & wait
echo ALLDONE > logs/fits_done.txt
