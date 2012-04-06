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

# Determine side based on the TASK ID
sid=$(((SGE_TASK_ID - 1) % 2))
if [[ $sid == 0 ]]; then
	side=left
	SL=$SLL
else
	side=right
	SL=$SLR
fi

# If no slice labeling, use automatic result
if [[ ! -f $SL ]]; then
  SL=$WORK/subfields/consensus_wgtavg_${side}_native.nii.gz
fi

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_biascorr_qsub.sh
	Root: ${ROOT?}
	Working directory: ${WORK?}
	Subjob ID: ${SGE_TASK_ID}
	Side: ${side?}
	Slice Labeling: ${SL}
	PATH: ${PATH?}
BLOCK1

# directory for the subfields (separate for different parameter values)
WSUB=$WORK/subfields

# Quit if slice labeling does not exist
if [[ ! -f $WSUB/consensus_heuristic_wgtavg_${side}_native.nii.gz ]]; then
	echo "Segmentation from MASV not available. Quitting"
	exit 0;
fi

# Call Hongzhi's BC code
$BIN/bc $WORK/tse.nii.gz $WSUB/consensus_heuristic_wgtavg_${side}_native.nii.gz \
  $sid $ROOT/data/adaboost $WSUB/bcfh_wgtavg_${side}_native.nii.gz

# Now apply subfield remapping (maintain heuristics)
$BIN/subfield_leveler $SL $WSUB/bcfh_wgtavg_${side}_native.nii.gz $WSUB/bcfh_heuristic_wgtavg_${side}_native.nii.gz

# Copy files into the 'final' directory
$BIN/c3d $WSUB/bcfh_heuristic_wgtavg_${side}_native.nii.gz -type ushort \
  -o $WORK/final/${SUBJID}_${side}_subfields_final.nii.gz

