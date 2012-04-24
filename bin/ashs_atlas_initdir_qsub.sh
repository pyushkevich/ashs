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

# Copy the images into the working directory
$ASHS_BIN/c3d -type ushort \
  $2 -o $MYWORK/mprage.nii.gz \
  $3 -o $MYWORK/tse.nii.gz \
  $4 -o $MYWORK/seg_left.nii.gz \
  $5 -o $MYWORK/seg_right.nii.gz

# Peform registration between the two modalities
ashs_align_t1t2 $MYWORK $WFSL
