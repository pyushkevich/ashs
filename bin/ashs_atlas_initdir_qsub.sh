#!/bin/bash
#$ -S /bin/bash
set -x -e

# Include the common file
source ashs_lib.sh

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_atlas_initdir_qsub.sh
	Root: ${ROOT?}
	Working directory: ${WORK?}
	PATH: ${PATH?}
BLOCK1

# Create the working directory for this ID
id=${1?}
MYWORK=$WORK/atlas/$id
WFSL=$MYWORK/flirt_t2_to_t1
mkdir -p $MYWORK $WFSL

# Copy the images into the working directory
$ASHS_BIN/c3d -type ushort \
  $2 -o $MYWORK/mprage.nii.gz \
  $3 -o $MYWORK/tse.nii.gz \
  $4 -o $MYWORK/seg_left.nii.gz \
  $5 -o $MYWORK/seg_right.nii.gz

# Peform registration between the two modalities
ashs_align_t1t2 $MYWORK $WFSL
