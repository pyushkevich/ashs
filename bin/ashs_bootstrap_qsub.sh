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
side=${1?}
tid=${2?}

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_bootstrap_qsub.sh
	Root: ${ASHS_ROOT?}
	Working directory: ${ASHS_WORK?}
	PATH: ${PATH?}
	Side: ${side?}
	Training subject: ${tid?}
BLOCK1

# Do we skip bootstrapping
if [[ $ASHS_NO_BOOTSTRAP -ne 1 ]]; then

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

  # Get the subject and atlas data
  ashs_subj_side_vars $ASHS_WORK $side
  ashs_atlas_side_vars $tid $side 0

  # Perform fit atlas segmentation to the MALF segmentation
  ml_affine $FDIR/lfseg_heur_${side}.nii.gz $ATLAS_SEG $WREG/bs_affine.mat

  # Halfway space registration

  # Split transform into halves
  c3d_affine_tool $WREG/bs_affine.mat -sqrt \
    -o $WREG/sqrt_fwd.mat -inv -o $WREG/sqrt_inv.mat 

  # Convert segmentation into a mesh
  c3d $FDIR/lfseg_heur_${side}.nii.gz -thresh 1 inf 1 0 -o $TMPDIR/mybin.nii.gz
  vtklevelset $TMPDIR/mybin.nii.gz $TMPDIR/mymesh.vtk 0.5

  # Transform one of the meshes forward
  warpmesh $TMPDIR/mymesh.vtk $TMPDIR/mymeshhw.vtk $WREG/sqrt_fwd.mat

  # Generate reference space
  RES=$(echo $ASHS_TEMPLATE_TARGET_RESOLUTION | sed -e "s/x/ /g" -e "s/mm//g")

  # TODO: this is a big margin - to do with greedy limitations!
  mesh2img -vtk $TMPDIR/mymeshhw.vtk -f -a $RES 5 $WREG/refspace.nii.gz

  # Warp both images into reference space
  c3d $WREG/refspace.nii.gz $ATLAS_TSE -reslice-matrix $WREG/sqrt_fwd.mat -o $TMPDIR/moving_hw.nii.gz
  c3d $WREG/refspace.nii.gz $SUBJ_TSE -reslice-matrix $WREG/sqrt_inv.mat -o $TMPDIR/fixed_hw.nii.gz

  # Create mask in reference space
  c3d $WREG/refspace.nii.gz  $FDIR/lfseg_heur_${side}.nii.gz -thresh 1 inf 1 0 \
    -int 0 -reslice-matrix $WREG/sqrt_inv.mat -thresh 0.5 inf 1 0 -dilate 1 10x10x2vox \
    -o $TMPDIR/mask.nii.gz

  # Bootstrap transform file
  BOOT_WARP=$WREG/greedy_warp.nii.gz

  # Report intermittent progress
  job_progress 0.5

  # Run ANTS in this space

  if [[ -f $BOOT_WARP && $ASHS_SKIP ]]; then
    echo "Skipping ANTS registration"
  else
    # TODO: use the right parameters!
    time greedy -d 3 $ASHS_GREEDY_THREADS \
      -gm $TMPDIR/mask.nii.gz \
      -m NCC 2x2x2 \
      -i $TMPDIR/fixed_hw.nii.gz $TMPDIR/moving_hw.nii.gz \
      -o $BOOT_WARP \
      -e 0.5 -n 60x60x20
  fi

  # Warp the moving ASHS_TSE image into the space of the native ASHS_TSE image using one interpolation.
  # Since we only care about the region around the segmentation, we use tse_native_chunk

  # Define the resliced images
  ATLAS_RESLICE=$WREG/atlas_to_native.nii.gz
  ATLAS_RESLICE_SEG=$WREG/atlas_to_native_segvote.nii.gz

  greedy -d 3 $ASHS_GREEDY_THREADS \
    -rm $ATLAS_TSE $ATLAS_RESLICE \
    -ri LABEL ${ASHS_LABEL_SMOOTHING} -rm $ATLAS_SEG $ATLAS_RESLICE_SEG \
    -rf $SUBJ_SIDE_TSE_NATCHUNK \
    -r $WREG/sqrt_inv.mat,-1 $BOOT_WARP $WREG/sqrt_fwd.mat

  # In tidy mode, we can clean up after this step
  if [[ $ASHS_TIDY ]]; then
    rm -rf $BOOT_WARP
  fi

fi

# Report progress
job_progress 1
