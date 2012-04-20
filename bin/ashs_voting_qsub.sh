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

source $ASHS_ROOT/bin/ashs_lib.sh

# Determine side based on the TASK ID
SIDES=(left right)
side=${SIDES[$(((SGE_TASK_ID - 1) % 2))]}

# Bootstrap mode or not
BOOTSTRAP=${1?}

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_voting_qsub.sh
	Root: ${ASHS_ROOT?}
	Working directory: ${ASHS_WORK?}
	Subjob ID: ${SGE_TASK_ID}
	Side: ${side?}
	PATH: ${PATH?}
	BOOTSTRAP: ${BOOTSTRAP}
BLOCK1

# Just run label fusion
ashs_label_fusion_apply $BOOTSTRAP
