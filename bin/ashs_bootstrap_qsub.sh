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

# Library
source ${ASHS_ROOT?}/bin/ashs_lib.sh

# Determine training case and side based on the TASK ID
SIDES=(left right)
side=${SIDES[$(((SGE_TASK_ID - 1) % 2))]}
tid=$(printf %03d $(((SGE_TASK_ID - 1)/ 2)))

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_multiatlas_qsub.sh
	Root: ${ASHS_ROOT?}
	Working directory: ${ASHS_WORK?}
	PATH: ${PATH?}
	Subjob ID: ${SGE_TASK_ID}
	Side: ${side?}
	Training subject: ${tid?}
BLOCK1

# Existing directory structure
WFSL=$ASHS_WORK/flirt_t2_to_t1
WANT=$ASHS_WORK/ants_t1_to_temp

# Training directory and training data
TDIR=$ASHS_ATLAS/train/train${tid}
TSEG=$TDIR/tse_native_chunk_${side}_seg.nii.gz

# Create directory for this registration
WREG=$ASHS_WORK/bootstrap/tseg_${side}_train${tid}
mkdir -p $WREG

# TODO: organize this better!!!

# Get the fusion directory
FDIR=$ASHS_WORK/multiatlas/fusion

# Perform fit atlas segmentation to the MALF segmentation

ml_affine $FDIR/lfseg_heur_${side}.nii.gz $TSEG $WREG/bs_affine.mat
c3d_affine_tool $WREG/bs_affine.mat -oitk $WREG/bs_affine_itk.txt

# Halfway space registration

# Split transform into halves
c3d_affine_tool $WREG/bs_affine.mat -sqrt \
	-o $WREG/sqrt_fwd.mat -oitk $WREG/sqrt_fwd_itk.txt \
	-inv -o $WREG/sqrt_inv.mat -oitk $WREG/sqrt_inv_itk.txt

# Convert segmentation into a mesh
c3d $FDIR/lfseg_heur_${side}.nii.gz -thresh 1 inf 1 0 -o $TMPDIR/mybin.nii.gz
vtklevelset $TMPDIR/mybin.nii.gz $TMPDIR/mymesh.vtk 0.5

# Transform one of the meshes forward
warpmesh $TMPDIR/mymesh.vtk $TMPDIR/mymeshhw.vtk $WREG/sqrt_fwd.mat

# Generate reference space
RES=$(echo $ASHS_TEMPLATE_TARGET_RESOLUTION | sed -e "s/x/ /g" -e "s/mm//g")
mesh2img -vtk $TMPDIR/mymeshhw.vtk -f -a $RES 5 $WREG/refspace.nii.gz

# Warp both images into reference space
c3d $WREG/refspace.nii.gz $TDIR/tse.nii.gz -reslice-matrix $WREG/sqrt_fwd.mat -o $TMPDIR/moving_hw.nii.gz
c3d $WREG/refspace.nii.gz $ASHS_WORK/tse.nii.gz -reslice-matrix $WREG/sqrt_inv.mat -o $TMPDIR/fixed_hw.nii.gz

# Create mask in reference space
c3d $WREG/refspace.nii.gz  $FDIR/lfseg_heur_${side}.nii.gz -thresh 1 inf 1 0 \
	-int 0 -reslice-matrix $WREG/sqrt_inv.mat -thresh 0.5 inf 1 0 -dilate 1 10x10x2vox \
	-o $TMPDIR/mask.nii.gz

# Run ANTS in this space
ANTS 3 \
  -x $TMPDIR/mask.nii.gz \
	-m PR[$TMPDIR/fixed_hw.nii.gz,$TMPDIR/moving_hw.nii.gz,1,4] \
	-o $WREG/antsreg.nii.gz \
	-i $ASHS_PAIRWISE_ANTS_ITER -t SyN[$ASHS_PAIRWISE_ANTS_STEPSIZE] -v \
	--continue-affine false

# Warp the moving ASHS_TSE image into the space of the native ASHS_TSE image using one interpolation.
# Since we only care about the region around the segmentation, we use tse_native_chunk
WarpImageMultiTransform 3 $TDIR/tse.nii.gz \
	$WREG/atlas_to_native.nii.gz \
	-R $ASHS_WORK/tse_native_chunk_${side}.nii.gz \
	-i $WREG/sqrt_inv_itk.txt \
	$WREG/antsregWarp.nii.gz \
	$WREG/sqrt_fwd_itk.txt

# Warp the segmentation labels the same way. This should work with WarpImageMultiTransform --use-ML
# but for some reason that is still broken. Let's use the old way
ATLAS_SEG=$TSEG
LSET=($(c3d $ATLAS_SEG -dup -lstat | awk 'NR > 1 {print $1}'))

for ((i=0; i < ${#LSET[*]}; i++)); do

	LID=$(printf '%03d' $i)
	c3d $ATLAS_SEG -thresh ${LSET[i]} ${LSET[i]} 1 0 -smooth 0.24mm -o $TMPDIR/label_${LID}.nii.gz

	WarpImageMultiTransform 3 $TMPDIR/label_${LID}.nii.gz \
		$TMPDIR/label_${LID}_warp.nii.gz \
		-R $ASHS_WORK/tse_native_chunk_${side}.nii.gz \
		-i $WREG/sqrt_inv_itk.txt \
		$WREG/antsregWarp.nii.gz \
		$WREG/sqrt_fwd_itk.txt

done

# Perform voting using replacement rules
RULES=$(for ((i=0; i < ${#LSET[*]}; i++)); do echo $i ${LSET[i]}; done)
c3d $TMPDIR/label_*_warp.nii.gz -vote -replace $RULES -o $WREG/atlas_to_native_segvote.nii.gz

# In tidy mode, we can clean up after this step
if [[ $ASHS_TIDY ]]; then
	rm -rf $WREG/antsreg*
fi

