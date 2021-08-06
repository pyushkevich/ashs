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
  Skip Ants: ${ASHS_SKIP_REGN}
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
if [[ $ASHS_MPRAGE -nt $SUBJ_RAWMPRAGE ]]; then
  ### TODO: I took this out because it was messing up CL with old atlases!
  ### c3d $ASHS_MPRAGE -stretch 1% 99% 0 4096 -clip 0 4096 -type short -o $SUBJ_MPRAGE
  c3d $ASHS_MPRAGE -type short -o $SUBJ_RAWMPRAGE

  # Perform preprocessing - denoising
  SUBJ_MPRAGE_DENOISE=$ASHS_WORK/mprage_d.nii.gz
  echo "ASHS_MPRAGE_DENOISE = $ASHS_MPRAGE_DENOISE"
  if [[ $ASHS_MPRAGE_DENOISE == 1 ]]; then
    NLMDenoise \
      -i $SUBJ_RAWMPRAGE \
      -o $SUBJ_MPRAGE_DENOISE
  else
    cp $SUBJ_RAWMPRAGE $SUBJ_MPRAGE_DENOISE
  fi

  # Perform preprocessing - SR upsample
  if [[ $ASHS_MPRAGE_SRUPSAMPLE == 1 ]]; then

    # get orientation code
    orient_code=$(c3d $ASHS_MPRAGE -info | cut -d ';' -f 5 | cut -d ' ' -f 5)
    if [[ $orient_code == "Oblique," ]]; then
      orient_code=$(c3d $ASHS_MPRAGE -info | cut -d ';' -f 5 | cut -d ' ' -f 8)
    fi

    # Swap the dimention of the image to RPI
    c3d $SUBJ_MPRAGE_DENOISE \
      -clip 0 inf -type short \
      -swapdim RPI \
      -o $SUBJ_MPRAGE_DENOISE

    # Perform upsample
    NLMUpsample \
      -i $SUBJ_MPRAGE_DENOISE \
      -o $SUBJ_MPRAGE \
      -lf $ASHS_MPRAGE_SRUPSAMPLE_FACTOR

    # Swap the dimention back to the original orientation code
    c3d $SUBJ_MPRAGE -swapdim $orient_code \
      -o $SUBJ_MPRAGE

  else
    cp $SUBJ_MPRAGE_DENOISE $SUBJ_MPRAGE
  fi

  # remove intermediate file
  rm -f $SUBJ_MPRAGE_DENOISE

fi

if [[ $ASHS_TSE -nt $SUBJ_RAWTSE ]]; then
  ### TODO: I took this out because it was messing up CL with old atlases!
  ### c3d $ASHS_TSE -stretch 1% 99% 0 4096 -clip 0 4096 -type short -o $SUBJ_TSE 
  c3d $ASHS_TSE -type short -o $SUBJ_RAWTSE 

  # Perform preprocessing - denoising
  SUBJ_TSE_DENOISE=$ASHS_WORK/tse_d.nii.gz
  if [[ $ASHS_TSE_DENOISE == 1 ]]; then
    NLMDenoise \
      -i $SUBJ_RAWTSE \
      -o $SUBJ_TSE_DENOISE
  else
    cp $SUBJ_RAWTSE $SUBJ_TSE_DENOISE
  fi

  # Perform preprocessing - SR upsample
  if [[ $ASHS_TSE_SRUPSAMPLE == 1 ]]; then

    # get orientation code
    orient_code=$(c3d $ASHS_TSE -info | cut -d ';' -f 5 | cut -d ' ' -f 5)
    if [[ $orient_code == "Oblique," ]]; then
      orient_code=$(c3d $ASHS_TSE -info | cut -d ';' -f 5 | cut -d ' ' -f 8)
    fi

    # Swap the dimention of the image to RPI
    c3d $SUBJ_TSE_DENOISE -swapdim RPI \
      -o $SUBJ_TSE_DENOISE

    # Perform upsample
    NLMUpsample \
      -i $SUBJ_TSE_DENOISE \
      -o $SUBJ_TSE \
      -lf $ASHS_TSE_SRUPSAMPLE_FACTOR

    # Swap the dimention back to the original orientation code
    c3d $SUBJ_TSE -clip 0 inf -type short -swapdim $orient_code \
      -o $SUBJ_TSE

  else
    cp $SUBJ_TSE_DENOISE $SUBJ_TSE
  fi

  # remove intermediate file
  rm -f $SUBJ_TSE_DENOISE


fi

# --- RIGID ALIGNMENT T1/T2 ---
ashs_align_t1t2 $ASHS_WORK $ASHS_INPUT_T2T1_MAT $ASHS_INPUT_T2T1_MODE

# Report some progress
if [[ $ASHS_MPRAGE_SRUPSAMPLE == 1 && $ASHS_TSE_SRUPSAMPLE == 1 ]]; then
  job_progress 0.8
elif [[ $ASHS_MPRAGE_SRUPSAMPLE == 1 || $ASHS_TSE_SRUPSAMPLE == 1 ]]; then
  job_progress 0.6
else
  job_progress 0.25
fi


# Use FLIRT to register T1 to template
if [[ -f $SUBJ_AFF_T1TEMP_MAT && $ASHS_SKIP_RIGID ]]; then

  echo "Skipping Affine Registration"

else

  if [[ $ASHS_T1TEMP_RIGID_MASK ]]; then
    BMASK=" -gm $TEMP_T1_MASK "
  else
    BMASK=""
  fi

  # Start with a very fast rigid transform - hope is it won't hurt things and
  # will actually prevent affine from doing insane scaling
  time greedy -d 3 $ASHS_GREEDY_THREADS -a -dof 6 -m NCC 2x2x2 \
    -i $TEMP_T1_FULL $SUBJ_MPRAGE \
    -o $WAFF/greedy_t1_to_template_init_rigid.mat -n 400x0x0x0 \
    -ia-image-centers -search 400 5 5 \
    $BMASK
    
  # Use greedy
  time greedy -d 3 $ASHS_GREEDY_THREADS -a -m NCC 2x2x2 \
    -i $TEMP_T1_FULL $SUBJ_MPRAGE \
    -o $WAFF/greedy_t1_to_template.mat -n 400x80x40x0 \
    -ia $WAFF/greedy_t1_to_template_init_rigid.mat \
    $BMASK

fi

greedy -d 3 $ASHS_GREEDY_THREADS -r $WAFF/greedy_t1_to_template.mat -rf $TEMP_T1_FULL \
  -rm $ASHS_WORK/mprage.nii.gz $WAFF/test_greedy_affine.nii.gz

# Store the transform
cp -a $WAFF/greedy_t1_to_template.mat $SUBJ_AFF_T1TEMP_MAT
ln -sf $WAFF/test_greedy_affine.nii.gz $SUBJ_AFF_T1TEMP_RESLICE

# Compute the inverse matrix
c3d_affine_tool $SUBJ_AFF_T1TEMP_MAT -inv -o $SUBJ_AFF_T1TEMP_INVMAT

# Report some more progress
if [[ $ASHS_MPRAGE_SRUPSAMPLE == 1 && $ASHS_TSE_SRUPSAMPLE == 1 ]]; then
  job_progress 0.9
elif [[ $ASHS_MPRAGE_SRUPSAMPLE == 1 || $ASHS_TSE_SRUPSAMPLE == 1 ]]; then
  job_progress 0.8
else
  job_progress 0.5
fi

# Use ANTS to warp the ASHS_MPRAGE image to the template
if [[ $ASHS_SKIP_REGN && -f $SUBJ_T1TEMP_WARP ]]; then

    echo "SKIPPING Deformable registration"
  
else
    
    time greedy -d 3 $ASHS_GREEDY_THREADS -m NCC 2x2x2 -e 0.5 -n ${ASHS_TEMPLATE_ITER} \
      -i $TEMP_T1_FULL $SUBJ_AFF_T1TEMP_RESLICE \
      -o $SUBJ_T1TEMP_WARP \
      -oinv $SUBJ_T1TEMP_INVWARP \
      -gm $TEMP_T1_MASK

fi

# Report some more progress
if [[ $ASHS_MPRAGE_SRUPSAMPLE == 1 && $ASHS_TSE_SRUPSAMPLE == 1 ]]; then
  job_progress 0.98
elif [[ $ASHS_MPRAGE_SRUPSAMPLE == 1 || $ASHS_TSE_SRUPSAMPLE == 1 ]]; then
  job_progress 0.95
else
  job_progress 0.9
fi

# Apply the transformation to the T1 image and to the T1 segmentation
greedy -d 3 $ASHS_GREEDY_THREADS -rm $SUBJ_MPRAGE $WANT/reslice_mprage_to_template.nii.gz \
  -rf $TEMP_T1_FULL -r $SUBJ_T1TEMP_WARP $SUBJ_AFF_T1TEMP_MAT

# Apply the transformation to the T2 image 
greedy -d 3 $ASHS_GREEDY_THREADS -rm $SUBJ_TSE $WANT/reslice_tse_to_template.nii.gz \
  -rf $TEMP_T1_FULL -r $SUBJ_T1TEMP_WARP $SUBJ_AFF_T1TEMP_MAT $SUBJ_AFF_T2T1_INVMAT

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

