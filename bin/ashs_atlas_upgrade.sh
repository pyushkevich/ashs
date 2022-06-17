#!/bin/bash - 

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

# This script updates the ASHS atlas to the current version

# Check the root dir
if [[ ! $ASHS_ROOT ]]; then
  echo "Please set ASHS_ROOT to the ASHS root directory before running $0"
  exit -2
fi

# Source the master script in order to get the PATH variable right
source ${ASHS_ROOT?}/bin/ashs_common_master.sh

# The atlas directory must be provided as the first input
ASHS_ATLAS=${1?  "Provide the path to the atlas to upgrade"}

# Check for the vars file
if [[ ! -f $ASHS_ATLAS/ashs_atlas_vars.sh ]]; then
  echo "Missing ashs_atlas_vars.sh in the atlas directory"
  exit -1
fi

# Get the number of atlases, other information
source $ASHS_ATLAS/ashs_atlas_vars.sh
TRIDS=$(for((i = 0; i < $ASHS_ATLAS_N; i++)); do echo $(printf "%03i" $i); done)

# Check for existence of ANTS-style warps, and rename them
for warp in $(find $ASHS_ATLAS -name 'ants_t1_to_chunktemp_*Warp.nii.gz'); do

  # Warp files are just copied to a new name (greedy_t1_to_template_${side}_warp.nii.gz)
  echo "Converting $warp"
  cp $warp $(dirname $warp)/$(echo $(basename $warp) | sed -e "s/ants_t1_to_chunktemp_/greedy_t1_to_template_/" -e "s/Warp.nii.gz/_warp.nii.gz/")

done

# Affine files need to be translated into the RAS format
for aff in $(find $ASHS_ATLAS -name 'ants_t1_to_tempAffine.txt'); do
  echo "Converting $aff"
  c3d_affine_tool -itk $aff -o $(dirname $aff)/t1_to_template_affine.mat
done

for aff in $(find $ASHS_ATLAS -name 'flirt_t2_to_t1_ITK.txt'); do
  echo "Converting $aff"
  c3d_affine_tool -itk $aff -o $(dirname $aff)/flirt_t2_to_t1.mat
done

# The version information in the atlas needs to be brought up to date
source ${ASHS_ROOT}/bin/ashs_version.sh
  cat > $ASHS_ATLAS/ashs_atlas_vars.sh <<-CONTENT
		ASHS_ATLAS_VERSION_FULL=${ASHS_VERSION_FULL}
		ASHS_ATLAS_VERSION_DATE=${ASHS_VERSION_DATE}
		ASHS_ATLAS_N=${ASHS_ATLAS_N}
	CONTENT

