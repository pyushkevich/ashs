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

# Read parameters
id=${1?}
tid=${2?}
side=${3?}

source ${ASHS_ROOT?}/bin/ashs_lib.sh

cd $ASHS_WORK/atlas/$id

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_atlas_pairwise.sh
	Root: ${ASHS_ROOT?}
	Working directory: ${ASHS_WORK?}
	PATH: ${PATH?}
	Subject: ${id?}
	Side: ${side?}
	Training subject: ${tid?}
BLOCK1

# Existing directory structure
WFSL=flirt_t2_to_t1
WANT=ants_t1_to_temp

# Training directory and training data
TDIR=$ASHS_WORK/atlas/${tid}
TSEG=$TDIR/seg_${side}.nii.gz

# Create directory for this registration
WREG=pairwise/tseg_${side}_train${tid}
mkdir -p $WREG

# Run ANTS with current image as fixed, training image as moving
ashs_ants_pairwise 1
