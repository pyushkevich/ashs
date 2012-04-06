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

source ashs_lib.sh

cd $WORK/atlas/$id

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_atlas_pairwise.sh
	Root: ${ROOT?}
	Working directory: ${WORK?}
	PATH: ${PATH?}
	Subject: ${id?}
	Side: ${side?}
	Training subject: ${tid?}
BLOCK1

# Existing directory structure
WFSL=flirt_t2_to_t1
WANT=ants_t1_to_temp

# Training directory and training data
TDIR=$WORK/atlas/${tid}
TSEG=$TDIR/seg_${side}.nii.gz

# Create directory for this registration
WREG=tseg_${side}_train${tid}
mkdir -p $WREG

# Run ANTS with current image as fixed, training image as moving
ashs_atlas_pairwise 1

<<'TRASH'
# Warp the moving TSE image into the space of the native TSE image using one interpolation.
# Since we only care about the region around the segmentation, we use tse_native_chunk
WarpImageMultiTransform 3 $TDIR/tse.nii.gz \
  $WREG/atlas_to_native.nii.gz \
  -R tse_native_chunk_${side}.nii.gz \
  -i flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt \
  -i ants_t1_to_temp/ants_t1_to_tempAffine.txt \
  ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii.gz \
	$WREG/antsregWarp.nii.gz \
	$WREG/antsregAffine.txt \
	$TDIR/ants_t1_to_temp/ants_t1_to_tempWarp.nii.gz \
	$TDIR/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
	$TDIR/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt

# Warp the segmentation labels the same way. This should work with WarpImageMultiTransform --use-ML
# but for some reason that is still broken. Let's use the old way
LSET=($(c3d $TDIR/seg_${side}.nii.gz -dup -lstat | awk 'NR > 1 {print $1}'))

for ((i=0; i < ${#LSET[*]}; i++)); do

  LID=$(printf '%03d' $i)
  c3d $TDIR/seg_${side}.nii.gz -thresh ${LSET[i]} ${LSET[i]} 1 0 -smooth 0.24mm -o $TMPDIR/label_${LID}.nii.gz

  WarpImageMultiTransform 3 $TMPDIR/label_${LID}.nii.gz \
    $TMPDIR/label_${LID}_warp.nii.gz \
    -R tse_native_chunk_${side}.nii.gz \
    -i flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt \
    -i ants_t1_to_temp/ants_t1_to_tempAffine.txt \
    ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii.gz \
    $WREG/antsregWarp.nii.gz \
    $WREG/antsregAffine.txt \
    $TDIR/ants_t1_to_temp/ants_t1_to_tempWarp.nii.gz \
    $TDIR/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
    $TDIR/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt

done

# Perform voting using replacement rules
RULES=$(for ((i=0; i < ${#LSET[*]}; i++)); do echo $i ${LSET[i]}; done)
c3d $TMPDIR/label_*_warp.nii.gz -vote -replace $RULES -o $WREG/atlas_to_native_segvote.nii.gz

TRASH

# THIS SEEMS BROKEN!
#WarpImageMultiTransform 3 $TDIR/seg_${side}.nii.gz \
#  $WREG/atlas_to_native_seg.nii.gz \
#  -R tse.nii.gz \
#  -i flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt \
#  -i ants_t1_to_temp/ants_t1_to_tempAffine.txt \
#  ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii.gz \
#	$WREG/antsregWarp.nii.gz \
#	$WREG/antsregAffine.txt \
#	$TDIR/ants_t1_to_temp/ants_t1_to_tempWarp.nii.gz \
#	$TDIR/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
#	$TDIR/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt \
#  --use-ML 0.24mm

<<'SKIP'
mkdir -p $WREG/sf
for((LAB=0;LAB<=10;LAB++)); do

	SFID=$(printf "sf%02i" $LAB)

	$BIN/c3d $TSEG -thresh $LAB $LAB 1 0 -smooth 0.24mm -o $WREG/sf/${SFID}_src.nii.gz

	# TODO: replace training full-warps with already combined warps residing in template space
	$BIN_ANTS/WarpImageMultiTransform 3 $WREG/sf/${SFID}_src.nii.gz \
    $WREG/sf/${SFID}_to_native.nii.gz -R $WORK/tse.nii.gz \
		-i $WORK/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt \
		-i $WORK/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
		$WORK/ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii \
		$WREG/antsregWarp.nii.gz \
		$TDIR/ants_t1_to_tempWarp.nii.gz \
		$TDIR/ants_t1_to_tempAffine.txt \
		$TDIR/flirt_t2_to_t1_ITK.txt

	$BIN_ANTS/WarpImageMultiTransform 3 $WREG/sf/${SFID}_src.nii.gz \
		$WREG/sf/${SFID}_to_chunk.nii.gz -R $WORK/tse_to_chunktemp_${side}.nii.gz \
		$WREG/antsregWarp.nii.gz \
		$TDIR/ants_t1_to_tempWarp.nii.gz \
		$TDIR/ants_t1_to_tempAffine.txt \
		$TDIR/flirt_t2_to_t1_ITK.txt

done

$BIN/c3d $WORK/tse_to_chunktemp_${side}.nii.gz $WREG/tse_reslice_to_chunk_oneinterp.nii.gz \
	-ncc 4x4x4 $WORK/tse_to_chunktemp_${side}_regmask.nii.gz -times \
	-o $WREG/tse_nccmap.nii.gz

SKIP
