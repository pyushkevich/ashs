#!/bin/bash
#$ -S /bin/bash
set -x -e

# Library
source $ROOT/bin/ashs_lib.sh

# Determine training case and side based on the TASK ID
SIDES=(left right)
side=${SIDES[$(((SGE_TASK_ID - 1) % 2))]}
tid=$(printf %03d $(((SGE_TASK_ID - 1)/ 2)))

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_multiatlas_qsub.sh
	Root: ${ROOT?}
	Working directory: ${WORK?}
	PATH: ${PATH?}
	Subjob ID: ${SGE_TASK_ID}
	Side: ${side?}
	Training subject: ${tid?}
BLOCK1

# Existing directory structure
WFSL=$WORK/flirt_t2_to_t1
WANT=$WORK/ants_t1_to_temp

# Training directory and training data
TDIR=$ASHS_ATLAS/train/train${tid}
TSEG=$TDIR/tse_native_chunk_${side}_seg.nii.gz

# Create directory for this registration
WREG=$WORK/multiatlas/tseg_${side}_train${tid}
mkdir -p $WREG

# Go to the work directory
cd $WORK

# Run ANTS with current image as fixed, training image as moving
ashs_ants_pairwise 0
