until [ -f out/fit_mann_atlas_jit_base.json ]; do sleep 10; done
OMP_NUM_THREADS=1 /private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/audit_env/bin/python s04_tests.py > logs/tests.log 2>&1
/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/audit_env/bin/python s03_jointfit.py mann claret jit raw 30000 > logs/fit_mann_claret_jit_raw.log 2>&1
/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/audit_env/bin/python s03_jointfit.py mann claret jit det98 30000 > logs/fit_mann_claret_jit_det98.log 2>&1
/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/audit_env/bin/python s03_jointfit.py uni claret jit raw 60000 > logs/fit_uni_claret_jit_raw.log 2>&1
/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/audit_env/bin/python s03_jointfit.py mann free jit raw 40000 > logs/fit_mann_free_jit_raw.log 2>&1
OMP_NUM_THREADS=1 /private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/audit_env/bin/python s06_depths_localbaseline.py mann_claret_jit_raw > logs/depths_localbaseline.log 2>&1
