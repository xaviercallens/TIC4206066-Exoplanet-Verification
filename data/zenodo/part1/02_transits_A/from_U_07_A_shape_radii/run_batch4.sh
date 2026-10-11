until [ -f out/fit_mann_free_jit_base.json ]; do sleep 10; done
/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/audit_env/bin/python s03_jointfit.py mann claret beta base 30000 > logs/fit_mann_claret_beta_base.log 2>&1
/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/audit_env/bin/python s03_jointfit.py mann claret jit b29 30000 > logs/fit_mann_claret_jit_b29.log 2>&1
/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/audit_env/bin/python s03_jointfit.py logu claret jit raw 60000 > logs/fit_logu_claret_jit_raw.log 2>&1
/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/audit_env/bin/python s03_jointfit.py uni free jit raw 60000 > logs/fit_uni_free_jit_raw.log 2>&1
EXCL_S98_EPOCHS=-1,8 /private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/audit_env/bin/python s03_jointfit.py mann claret jit base 30000 > logs/fit_mann_claret_jit_base_excl.log 2>&1
EXCL_S98_EPOCHS=-1,8 /private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/audit_env/bin/python s03_jointfit.py mann claret jit raw 30000 > logs/fit_mann_claret_jit_raw_excl.log 2>&1
