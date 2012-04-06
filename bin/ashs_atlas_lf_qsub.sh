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

source ashs_lib.sh


# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_atlas_pairwise.sh
	Root: ${ROOT?}
	Working directory: ${WORK?}
	PATH: ${PATH?}
	Subject: ${id?}
  Train Set: ${TRAIN?}
  Output: ${FNOUT?}
	Side: ${side?}
BLOCK1

# Just call the label fusion library command
ashs_label_fusion $id $FNOUT
