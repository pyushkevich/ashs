#!/bin/bash

# Common code executed at startup of all 'master' scripts, i.e., ashs_main,
# ashs_train

# The ROOT variable must be set
if [[ ! $ASHS_ROOT ]]; then
  echo "ASHS_ROOT is not set. Please set this variable to point to the root ASHS directory"
  exit -1;
fi

# Get the architecture and check ability to run binaries
ARCH=$(uname);
ASHS_BIN=$ASHS_ROOT/ext/$ARCH/bin
#ASHS_ANTS=$ASHS_BIN/ants
ASHS_ANTS=$ASHS_BIN/ants_1042
ASHS_FSL=$ASHS_BIN/fsl
if [[ ! $($ASHS_BIN/c3d -version | grep 'Version') ]]; then
  echo "Can not execute command \'$ASHS_BIN/c3d -version\'. Wrong architecture?"
  exit -1
fi

# Set the path for the ASHS programs, to ensure that we don't call some other
# version of ants or c3d. This is nicer than having to prefix every call to
# c3d by the path. 
PATH=$ASHS_BIN:$ASHS_ANTS:$ASHS_FSL:$PATH

# Check that the data directory exists
if [[ ! -f $ASHS_ROOT/data/train/train21/tse_native.nii.gz ]]; then
  echo "Data files appear to be missing. Can't locate $ASHS_ROOT/data/train/train21/tse_native.nii.gz"
  exit -1
fi
