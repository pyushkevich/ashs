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
WAFF=$ASHS_WORK/affine_t1_to_template
mkdir -p $WANT $WFSL $WAFF

# Set some variables
TEMP_T1_FULL=$ASHS_ATLAS/template/template.nii.gz
TEMP_T1_MASK=$ASHS_ATLAS/template/template_bet_mask.nii.gz

# Subject files generated in this script
ashs_subj_vars $ASHS_WORK

# Additional variables for local use
SUBJ_AFF_T1TEMP_RESLICE=$WAFF/t1_to_template_affine.nii.gz

# Copy the images into the working directory
if [[ $ASHS_MPRAGE -nt $SUBJ_MPRAGE ]]; then
  ### TODO: I took this out because it was messing up CL with old atlases!
  ### c3d $ASHS_MPRAGE -stretch 1% 99% 0 4096 -clip 0 4096 -type short -o $SUBJ_MPRAGE
  c3d $ASHS_MPRAGE -type short -o $SUBJ_MPRAGE
fi

if [[ $ASHS_TSE -nt $SUBJ_TSE ]]; then
  ### TODO: I took this out because it was messing up CL with old atlases!
  ### c3d $ASHS_TSE -stretch 1% 99% 0 4096 -clip 0 4096 -type short -o $SUBJ_TSE 
  c3d $ASHS_TSE -type short -o $SUBJ_TSE 
fi

# --- RIGID ALIGNMENT T1/T2 ---
ashs_align_t1t2 $ASHS_WORK $WFSL

# Report some progress
job_progress 0.25

# Use FLIRT to register T1 to template
if [[ -f $SUBJ_AFF_T1TEMP_MAT && $ASHS_SKIP_RIGID ]]; then

  echo "Skipping Affine Registration"

else


  # Start with a very fast rigid transform - hope is it won't hurt things and
  # will actually prevent affine from doing insane scaling
  time greedy -d 3 -a -dof 6 -m NCC 2x2x2 \
    -i $TEMP_T1_FULL $SUBJ_MPRAGE \
    -o $WAFF/greedy_t1_to_template_init_rigid.mat -n 400x0x0x0 \
    -ia-image-centers -search 400 5 5
    
  # Use greedy
  time greedy -d 3 -a -m NCC 2x2x2 \
    -i $TEMP_T1_FULL $SUBJ_MPRAGE \
    -o $WAFF/greedy_t1_to_template.mat -n 400x80x40x0 \
    -ia $WAFF/greedy_t1_to_template_init_rigid.mat

  greedy -d 3 -r $WAFF/greedy_t1_to_template.mat -rf $TEMP_T1_FULL \
    -rm $ASHS_WORK/mprage.nii.gz $WAFF/test_greedy_affine.nii.gz

  # Store the transform
  cp -a $WAFF/greedy_t1_to_template.mat $SUBJ_AFF_T1TEMP_MAT
  ln -sf $WAFF/test_greedy_affine.nii.gz $SUBJ_AFF_T1TEMP_RESLICE

  # Compute the inverse matrix
  c3d_affine_tool $SUBJ_AFF_T1TEMP_MAT -inv -o $SUBJ_AFF_T1TEMP_INVMAT

fi

# Report some more progress
job_progress 0.5

# Use ANTS to warp the ASHS_MPRAGE image to the template
if [[ $ASHS_SKIP_ANTS && -f $SUBJ_T1TEMP_WARP ]]; then

    echo "SKIPPING Deformable registration"
  
else
    
    time greedy -d 3 -m NCC 2x2x2 -e 0.5 -n ${ASHS_TEMPLATE_ANTS_ITER} \
      -i $TEMP_T1_FULL $SUBJ_AFF_T1TEMP_RESLICE \
      -o $SUBJ_T1TEMP_WARP \
      -oinv $SUBJ_T1TEMP_INVWARP \
      -gm $TEMP_T1_MASK

fi

# Report some more progress
job_progress 0.9

# Apply the transformation to the T1 image and to the T1 segmentation
greedy -d 3 -rm $SUBJ_MPRAGE $WANT/reslice_mprage_to_template.nii.gz \
  -rf $TEMP_T1_FULL -r $SUBJ_T1TEMP_WARP $SUBJ_AFF_T1TEMP_MAT

# Apply the transformation to the T2 image 
greedy -d 3 -rm $SUBJ_TSE $WANT/reslice_tse_to_template.nii.gz \
  -rf $TEMP_T1_FULL -r $SUBJ_T1TEMP_WARP $SUBJ_AFF_T1TEMP_MAT $SUBJ_AFF_T2T1_MAT

# Transform all of the images into the chunk template space
ashs_reslice_to_template $ASHS_WORK $ASHS_ATLAS

# If there is a reference segmentation image, process it too
if [[ -f $ASHS_REFSEG_RIGHT && -f $ASHS_REFSEG_LEFT ]]; then
  mkdir -p $ASHS_WORK/refseg
  c3d $ASHS_REFSEG_LEFT -type short -o $ASHS_WORK/refseg/refseg_left.nii.gz
  c3d $ASHS_REFSEG_RIGHT -type short -o $ASHS_WORK/refseg/refseg_right.nii.gz

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

# Generate a QC image for the registration
for side in $SIDES; do
  ashs_registration_qc $side
done

# Report final progress
job_progress 1

