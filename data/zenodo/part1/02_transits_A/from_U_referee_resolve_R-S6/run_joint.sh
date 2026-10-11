#!/bin/bash
# Queue of joint fits for R-S6; run with at most 2 parallel processes.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
PY=/private/tmp/claude-501/-Users-paulsmith-space--3/327b8aec-95f9-4559-ac96-4dbbd1cede72/scratchpad/audit_env/bin/python
cat > logs/joint_queue.txt <<'Q'
6,32,98 raw jit
32,98 raw jit
6,32,98 raw jit SPLIT6=1
6,32,98 raw jit SS6=1
6,32,98 raw jit SPLIT6=1 SS6=1
32,98 raw jit EXCL=98:-1
6,32,98 raw jit EXCL=98:-1
6,32,98 bw jit
32,98 bw jit
6,32,98 bw beta
32,98 bw beta
6,32,98 bw jit SPLIT6=1
Q
tr '\n' '\0' < logs/joint_queue.txt | xargs -0 -P 2 -I{} bash -c 'args="{}"; tag=$(echo "$args" | tr " ,=:" "____"); $0 s2_joint.py $args > logs/joint_$tag.log 2>&1; echo "done $args"' "$PY"
echo ALL_DONE
