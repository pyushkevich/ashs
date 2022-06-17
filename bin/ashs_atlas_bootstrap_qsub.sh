#!/bin/bash
#$ -S /bin/bash

#######################################################################
#
#  Program:   ASHS (Automatic Segmentation of Hippocampal Subfields)
#  Module:    $Id: ashs_atlas_loo_qsub.sh 86 2012-04-24 17:23:00Z yushkevich $
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

# Read parameters
id=${1?}
side=${2?}
XID=${3?}
tid=${4?}

source ${ASHS_ROOT?}/bin/ashs_lib.sh

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_atlas_bootstrap.sh
	Root: ${ASHS_ROOT?}
	Working directory: ${ASHS_WORK?}
	PATH: ${PATH?}
	Subject: ${id?}
  Training case: ${tid?}
  Crossval Id: ${XID?}
	Side: ${side?}
BLOCK1

# Create output directory
WBS=$ASHS_WORK/xval/${XID}_${side}/bootstrap/test_${id}/tseg_${side}_train${tid}
mkdir -p $WBS

# The uncorrected segmentation for this iteration
SEGUNC=$WXVAL/${id}_${side}_xval_seg_uncorr.nii.gz
SEGTRUTH=$ASHS_WORK/atlas/$id/tse_native_chunk_${side}_seg.nii.gz

# Create the posterior maps
POSTUNC=$WXVAL/${id}_${side}_xval_posterior_%03d_uncorr.nii.gz

# Call the label fusion library command
ashs_label_fusion $id $SEGUNC $POSTUNC

# Apply the correction
pushd $ASHS_WORK/train/${XID}_${side}

# If there are heuristics, make sure they are supplied to the LF program
if [[ $ASHS_HEURISTICS ]]; then
  pushd $ASHS_WORK/atlas/$id
  EXCLCMD=$(for fn in $(ls $ASHS_WORK/atlas/$id/heurex/heurex_${side}_*.nii.gz); do \
    echo "-x $(echo $fn | sed -e "s/.*_//g" | awk '{print 1*$1}') $fn"; \
    done)
  popd
fi

sa $ASHS_WORK/atlas/$id/tse_native_chunk_${side}.nii.gz $SEGUNC adaboost $SEGCOR $EXCLCMD

popd

# Compute overlaps
LSET=($(c3d $SEGTRUTH -dup -lstat | awk 'NR > 1 {print $1}'))

# Compute overlaps 
for type in corr uncorr; do
  for lab in ${LSET[*]}; do
    OVL=$(c3d $WXVAL/${id}_${side}_xval_seg_${type}.nii.gz $SEGTRUTH -overlap $lab | awk '$1=="OVL:" {print 1.0 * $6}')
    echo $idtest $side $ixval $lab $OVL
  done > $WXVAL/${id}_${side}_xval_ovl_${type}.txt

  # For both corrected and uncorrected, resample the segmentation to native
  # space, full size.
  c3d $ASHS_WORK/atlas/$id/tse.nii.gz $WXVAL/${id}_${side}_xval_seg_${type}.nii.gz \
    -int 0 -reslice-identity  -o $WXVAL/${id}_${side}_xval_seg_${type}_fullsize.nii.gz

done
