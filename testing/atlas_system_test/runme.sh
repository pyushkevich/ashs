#!/bin/bash
echo "ASHS_ROOT:   ${ASHS_ROOT?}"
echo "Output Dir:  ${1?}"

OUTDIR=${1?}

$ASHS_ROOT/bin/ashs_train \
  -D config/manifest.txt \
  -L config/snaplabels.txt \
  -w $OUTDIR \
  -d \
  -r config/rules.txt \
  -x config/xval.txt \
  -C config/ashs_config_test.txt | tee $OUTDIR/ashs_train_stdout.txt

