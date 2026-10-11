#!/bin/bash
# reduced second queue (after the two running fits finish); max 2 concurrent processes
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/audit_env/bin/python
N=15000
while pgrep -f "v03_fit.py o2_hw0.45_SYS_RM_FLUX (uni|logu) claret 98:-1" > /dev/null; do sleep 5; done
q1() {
$PY v04_tests.py 0.0212 0.35 2.51 0.316 0.313 SYS_RM_FLUX ms_template > logs/t1.log 2>&1
$PY v05_vshape.py o2_hw0.45_SYS_RM_FLUX none > logs/v1.log 2>&1
$PY v05_vshape.py o2_hw0.45_SYS_RM_FLUX 98:-1 > logs/v2.log 2>&1
$PY v03_fit.py o2_hw0.45_SYS_RM_FLUX mann free 98:-1 $N > logs/f11.log 2>&1
$PY v03_fit.py o2_hw0.45_DET_FLUX mann claret 98:-1 $N > logs/f6.log 2>&1
$PY v03_fit.py o1_hw0.45_SYS_RM_FLUX mann claret 98:-1 $N > logs/f4.log 2>&1
}
q2() {
$PY v03_fit.py o2_hw0.45_SYS_RM_FLUX uni free 98:-1 $N > logs/f12.log 2>&1
$PY v03_fit.py o2_hw0.45_DET_FLUX mann claret none $N > logs/f5.log 2>&1
$PY v03_fit.py o2_hw0.45_DET_FLUX uni claret none $N > logs/f13.log 2>&1
}
q1 & q2 & wait
echo ALLDONE > logs/fits2_done.txt
