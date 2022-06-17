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

source ${ASHS_ROOT?}/bin/ashs_lib.sh

# Bootstrap mode or not
BOOTSTRAP=${1?}
side=${2?}

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_voting_qsub.sh
	Root: ${ASHS_ROOT?}
	Working directory: ${ASHS_WORK?}
	Side: ${side?}
	PATH: ${PATH?}
	BOOTSTRAP: ${BOOTSTRAP}
BLOCK1

if [[ $BOOTSTRAP -ne 1 || $ASHS_NO_BOOTSTRAP -ne 1 ]]; then

  # Just run label fusion
  ashs_label_fusion_apply $BOOTSTRAP

  # Generate the QC images
  MALFMODE=$(if [[ $BOOTSTRAP -eq 1 ]]; then echo bootstrap; else echo multiatlas; fi)
  ashs_segmentation_qc $side $MALFMODE heur
  ashs_segmentation_qc $side $MALFMODE corr_usegray
  ashs_segmentation_qc $side $MALFMODE corr_nogray

fi

# Report progress
job_progress 1
