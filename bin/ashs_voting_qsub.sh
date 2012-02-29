#!/bin/bash
#$ -S /bin/bash
set -x -e

source $ROOT/bin/ashs_lib.sh

# Determine side based on the TASK ID
SIDES=(left right)
side=${SIDES[$(((SGE_TASK_ID - 1) % 2))]}

# Bootstrap mode or not
BOOTSTRAP=${1?}

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_voting_qsub.sh
	Root: ${ROOT?}
	Working directory: ${WORK?}
	Subjob ID: ${SGE_TASK_ID}
	Side: ${side?}
	PATH: ${PATH?}
	BOOTSTRAP: ${BOOTSTRAP}
BLOCK1

# Just run label fusion
ashs_label_fusion_apply $BOOTSTRAP
