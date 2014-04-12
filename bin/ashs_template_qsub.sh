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

# Read the library
source ${ASHS_ROOT?}/bin/ashs_lib.sh

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_template_qsub.sh
	Root: ${ASHS_ROOT?}
	Working directory: ${ASHS_WORK?}
	PATH: ${PATH?}
  Skip Rigid: ${ASHS_SKIP_RIGID}
  Skip Ants: ${ASHS_SKIP_ANTS}
BLOCK1

# Ensure directory structure
WFSL=$ASHS_WORK/flirt_t2_to_t1
WANT=$ASHS_WORK/ants_t1_to_temp
mkdir -p $WANT $WFSL

# Set some variables
TEMP_T1_FULL=$ASHS_ATLAS/template/template.nii.gz
TEMP_T1_MASK=$ASHS_ATLAS/template/template_bet_mask.nii.gz

# Copy the images into the working directory
if [[ $ASHS_MPRAGE -nt $ASHS_WORK/mprage.nii.gz ]]; then
  c3d $ASHS_MPRAGE -type ushort -o $ASHS_WORK/mprage.nii.gz
fi

if [[ $ASHS_TSE -nt $ASHS_WORK/tse.nii.gz ]]; then
  c3d $ASHS_TSE -type ushort -o $ASHS_WORK/tse.nii.gz 
fi

# Histogram match the images to a reference image (used later down the road, but better to do it now)
for kind in tse mprage; do
  c3d $ASHS_ATLAS/ref_hm/ref_${kind}.nii.gz $ASHS_WORK/${kind}.nii.gz \
    -histmatch 5 -o $ASHS_WORK/${kind}_histmatch.nii.gz
done

# --- RIGID ALIGNMENT T1/T2 ---
ashs_align_t1t2 $ASHS_WORK $WFSL

# Use FLIRT to register T1 to template
WAFF=$ASHS_WORK/affine_t1_to_template
mkdir -p $WAFF

if [[ -f $WAFF/t1_to_template_ITK.txt && $ASHS_SKIP_RIGID ]]; then

  echo "Skipping Affine Registration"

else

	# Try using FLIRT
  export FSLOUTPUTTYPE=NIFTI_GZ

  # Run flirt with template as reference
  flirt -v -anglerep quaternion -ref $TEMP_T1_FULL -in $ASHS_WORK/mprage.nii.gz \
    -datatype short -o $WAFF/test_flirt_affine.nii.gz \
    -omat $WAFF/flirt_intermediate_affine.mat -cost corratio -searchcost corratio

  # Convert the transform to ITK
  c3d_affine_tool $WAFF/flirt_intermediate_affine.mat -ref $TEMP_T1_FULL -src $ASHS_WORK/mprage.nii.gz \
    -fsl2ras -oitk $WAFF/flirt_t1_to_template_ITK.txt

  # Try using ANTS
	ANTS 3 -m MI[$TEMP_T1_FULL,$ASHS_WORK/mprage.nii.gz,1,32] -o $WAFF/antsaffineonly.nii.gz -i 0 
	WarpImageMultiTransform 3 $ASHS_WORK/mprage.nii.gz $WAFF/test_ants_affine.nii.gz \
		-R $TEMP_T1_FULL $WAFF/antsaffineonlyAffine.txt

	# Compare the two images
	SIM_FLIRT=$(c3d $TEMP_T1_FULL $WAFF/test_flirt_affine.nii.gz -nmi | awk '{print int(1000*$3)}');
	SIM_ANTS=$(c3d $TEMP_T1_FULL $WAFF/test_ants_affine.nii.gz -nmi | awk '{print int(1000*$3)}');
	if [[ $SIM_FLIRT -gt $SIM_ANTS ]]; then
		cp -a $WAFF/flirt_t1_to_template_ITK.txt $WAFF/t1_to_template_ITK.txt
	else
		cp -a $WAFF/antsaffineonlyAffine.txt $WAFF/t1_to_template_ITK.txt
	fi

fi

# Use ANTS to warp the ASHS_MPRAGE image to the template
if [[ $ASHS_SKIP_ANTS \
  && -f $WANT/ants_t1_to_tempAffine.txt \
  && -f $WANT/ants_t1_to_tempWarp.nii \
  && -f $WANT/ants_t1_to_tempInverseWarp.nii ]]; then

    echo "SKIPPING ANTS"
  
else

    ANTS 3 -m PR[$TEMP_T1_FULL,$ASHS_WORK/mprage.nii.gz,1,4] \
      -x $TEMP_T1_MASK \
      -o $WANT/ants_t1_to_temp.nii \
      -a $WAFF/t1_to_template_ITK.txt \
      -i ${ASHS_TEMPLATE_ANTS_ITER} -v | tee $WANT/ants_output.txt
  
    shrink_warp 3 $WANT/ants_t1_to_tempWarp.nii.gz $WANT/ants_t1_to_tempWarp.nii.gz
    shrink_warp 3 $WANT/ants_t1_to_tempInverseWarp.nii.gz $WANT/ants_t1_to_tempInverseWarp.nii.gz

fi

# Apply the transformation to the T1 image and to the T1 segmentation
WarpImageMultiTransform 3 $ASHS_WORK/mprage.nii.gz \
  $WANT/reslice_mprage_to_template.nii.gz -R $TEMP_T1_FULL \
  $WANT/ants_t1_to_tempWarp.nii $WANT/ants_t1_to_tempAffine.txt

# Apply the transformation to the T2 image 
WarpImageMultiTransform 3 $ASHS_WORK/tse.nii.gz $WANT/reslice_tse_to_template.nii.gz -R $TEMP_T1_FULL \
  $WANT/ants_t1_to_tempWarp.nii $WANT/ants_t1_to_tempAffine.txt $WFSL/flirt_t2_to_t1_ITK.txt

# Transform all of the images into the chunk template space
ashs_reslice_to_template $ASHS_WORK $ASHS_ATLAS

# If there is a reference segmentation image, process it too
if [[ -f $ASHS_REFSEG_RIGHT && -f $ASHS_REFSEG_LEFT ]]; then
  mkdir -p $ASHS_WORK/refseg
  c3d $ASHS_REFSEG_LEFT -type ushort -o $ASHS_WORK/refseg/refseg_left.nii.gz
  c3d $ASHS_REFSEG_RIGHT -type ushort -o $ASHS_WORK/refseg/refseg_right.nii.gz

  if [[ $ASHS_HEURISTICS ]]; then
    for side in left right; do

      c3d $ASHS_WORK/tse_native_chunk_${side}.nii.gz $ASHS_WORK/refseg/refseg_${side}.nii.gz \
        -int 0 -reslice-identity -type short -o $ASHS_WORK/refseg/refseg_native_chunk_${side}.nii.gz
    
      mkdir -p $ASHS_WORK/refseg/heurex
      subfield_slice_rules $ASHS_WORK/refseg/refseg_native_chunk_${side}.nii.gz \
        $ASHS_HEURISTICS $ASHS_WORK/refseg/heurex/heurex_${side}_%04d.nii.gz

    done
  fi
fi

# We don't need to keep the histmatch images around
if [[ $ASHS_TIDY ]]; then
	rm -rf $ASHS_WORK/tse_histmatch.nii.gz $ASHS_WORK/mprage_histmatch.nii.gz
fi
