#!/bin/bash
#$ -S /bin/bash
set -x -e

source ashs_lib.sh


# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_atlas_pairwise.sh
	Root: ${ROOT?}
	Working directory: ${WORK?}
	PATH: ${PATH?}
	Subject: ${id?}
  Train Set: ${TRAIN?}
  Output: ${FNOUT?}
	Side: ${side?}
BLOCK1

# Just call the label fusion library command
ashs_label_fusion $id $FNOUT
