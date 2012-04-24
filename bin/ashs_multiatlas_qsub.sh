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
WREG=$ASHS_WORK/multiatlas/tseg_${side}_train${tid}
mkdir -p $WREG

# Go to the work directory
cd $ASHS_WORK

# Run ANTS with current image as fixed, training image as moving
ashs_ants_pairwise 0
