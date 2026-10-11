#!/bin/zsh
# sequential TRICERATOPS runs (each uses all cores via parallel=True)
PY=/private/tmp/claude-501/-Users-paulsmith-space--3/21ae3b23-e2c8-49b9-927a-27d6e2d5b6a3/scratchpad/trice_env/bin/python
cd "/Users/paulsmith/space <3/tess-jev/addenda/tic4206066_PUBLICATION_AUDIT/scripts"
for spec in "A repro S32+S98" "B repro S32+S98" "A base S32+S98" "B base S32+S98" "A loc S32+S98" "B loc S32+S98" \
            "A molusc5 S32+S98" "B molusc5 S32+S98" "A molusc1 S32+S98" "B molusc1 S32+S98" "A tic S32+S98" "B tic S32+S98" \
            "A base S32only" "A base S98only" "B base S98only" "A base S6only"; do
  set -- ${=spec}
  $PY s13_fpp_run.py $1 $2 $3 20 1000000 > ../logs/fpp_$1_$2_$3.log 2>&1
done
echo ALLDONE > ../logs/fpp_alldone.flag
