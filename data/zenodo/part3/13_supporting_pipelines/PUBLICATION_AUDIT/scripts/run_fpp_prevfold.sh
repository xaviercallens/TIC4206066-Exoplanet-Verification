#!/bin/zsh
# queued after run_fpp_all.sh: earlier-fold FPP tests for B (s13c)
PY=/private/tmp/claude-501/-Users-paulsmith-space--3/21ae3b23-e2c8-49b9-927a-27d6e2d5b6a3/scratchpad/trice_env/bin/python
cd "/Users/paulsmith/space <3/tess-jev/addenda/tic4206066_PUBLICATION_AUDIT/scripts"
until [ -f ../logs/fpp_alldone.flag ]; do sleep 20; done
for v in prevfold_repro prevfold_base; do
  $PY s13c_fpp_prevfold.py $v 20 1000000 > ../logs/fpp_B_${v}_S32+S98.log 2>&1
done
echo DONE > ../logs/fpp_prevfold_done.flag
