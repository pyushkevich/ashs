#!/bin/bash
#$ -S /bin/bash

#######################################################################
#
#  Program:   ASHS (Automatic Segmentation of Hippocampal Subfields)
#  Module:    $Id$
#  Language:  BASH Shell Script
#  Copyright (c) 2012 Paul A. Yushkevich, University of Pennsylvania
#  
#  This file is part of ASHS
#
#  ASHS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details. 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

set -x -e

# Include the common file
source ${ASHS_ROOT?}/bin/ashs_lib.sh

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_atlas_initdir_qsub.sh
	Root: ${ASHS_ROOT?}
	Working directory: ${ASHS_WORK?}
	PATH: ${PATH?}
BLOCK1

# Create the working directory for this ID
id=${1?}
MYWORK=$ASHS_WORK/atlas/$id
WFSL=$MYWORK/flirt_t2_to_t1
mkdir -p $MYWORK $WFSL

# Parse the manifest file for this ID
FILES=($(cat $ASHS_TRAIN_MANIFEST | awk -v id=$id '$1 == id { print $2,$3,$4,$5 }'))

# Copy the images into the working directory and set the transforms
# of the segmentations to equal the transforms of the input images.
if [[ $ASHS_SKIP && \
      -f $MYWORK/mprage.nii.gz && -f $MYWORK/tse.nii.gz && \
      -f $MYWORK/seg_left.nii.gz && -f $MYWORK/seg_right.nii.gz && \
      -f $MYWORK/mprage_lr.nii.gz ]]; 
then
  echo "Skipping initial image copy for subject $id"
else 
  $ASHS_BIN/c3d -type ushort \
    ${FILES[0]} -as MPR -o $MYWORK/mprage.nii.gz \
    ${FILES[1]} -o $MYWORK/tse.nii.gz -popas REF \
    -push REF ${FILES[2]} -copy-transform -o $MYWORK/seg_left.nii.gz \
    -push REF ${FILES[3]} -copy-transform -o $MYWORK/seg_right.nii.gz \
    -push MPR -smooth-fast 4vox -resample 12.5% -o $MYWORK/mprage_lr.nii.gz
fi

# Peform registration between the two modalities
# Check if the is an override for the affine matrix
if [[ -f $ASHS_TRAIN_TRANSFORM_MANIFEST ]]; then
  OVERRIDE_MAT=$(cat $ASHS_TRAIN_TRANSFORM_MANIFEST | awk -v id=$id '$1==id {print $2}')
  OVERRIDE_MODE=$(cat $ASHS_TRAIN_TRANSFORM_MANIFEST | awk -v id=$id '$1==id {print $3}')
  ashs_align_t1t2 $MYWORK $OVERRIDE_MAT $OVERRIDE_MODE
else
  ashs_align_t1t2 $MYWORK
fi

