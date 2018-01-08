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

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_extractstats_qsub.sh
	Root: ${ASHS_ROOT?}
	Working directory: ${ASHS_WORK?}
	PATH: ${PATH?}
  SUBJECT: ${ASHS_SUBJID?}
BLOCK1

function voxel_size()
{
	echo $(c3d $1 -info-full | grep Spacing | sed -e "s/[^0-9\.]/ /g" | awk '{print $1*$2*$3}')
}

# directory for the subfields (separate for different parameter values)
WSTAT=$ASHS_WORK/final
mkdir -p $WSTAT

# Generate a list of subfield indices and names. This awk command skips comment lines,
# ignores the zero label, extracts the label name from quotes, and replaces whitespace
# in the label name with underlines
cat $ASHS_ATLAS/snap/snaplabels.txt | \
  awk '$1 > 0 {split($0,arr,"\""); sub(/[ \t]+/,"_",arr[2]); print $1,arr[2]}' \
  > $TMPDIR/labels.txt

LABIDS=($(cat $TMPDIR/labels.txt | awk '{print $1}'))
LABNAMES=($(cat $TMPDIR/labels.txt | awk '{print $2}'))

# Names of segmentations
for segtype in raw heur corr_usegray corr_nogray manual; do
  for side in $SIDES; do

    if [[ $segtype == "manual" ]]; then
      SBC=$ASHS_WORK/refseg/refseg_${side}.nii.gz
    else
      SBC=$ASHS_WORK/bootstrap/fusion/lfseg_${segtype}_${side}.nii.gz
    fi

    if [[ -f $SBC ]]; then

      # Generate the voxel and extent statistics
      STATS=$TMPDIR/rawvols_${segtype}_${side}.txt
      c3d $SBC -dup -lstat | tee $STATS

      # Create an output file
      FNBODYVOL=$WSTAT/${ASHS_SUBJID}_${side}_${segtype}_volumes.txt 
      rm -rf $FNBODYVOL

      # Dump volumes into that file
      for ((ilab = 0; ilab < ${#LABIDS[*]}; ilab++)); do

        # The id of the label
        i=${LABIDS[ilab]};
        SUB=${LABNAMES[ilab]};

        # Get the extent along z axis
        NBODY=$(cat $STATS | awk -v id=$i '$1 == id {print $10}')

        # Get the volume of this subfield
        VSUB=$(cat $STATS | awk -v id=$i '$1 == id {print $7}')

        # Write the volume information to output file
        if [[ $NBODY ]]; then
          echo $ASHS_SUBJID $side $SUB $NBODY $VSUB >> $FNBODYVOL
        fi

      done

      # If there is a reference segmentation, generate overlap statistics
      REFSEG=$ASHS_WORK/refseg/refseg_${side}.nii.gz
      if [[ -f $REFSEG && $segtype != "manual" ]]; then

        # Get the overlap statistics
        OVLFILE=$TMPDIR/ovl_${segtype}_${side}.txt
        c3d $REFSEG -int 0 -dup $SBC -reslice-identity -label-overlap > $OVLFILE

        # Extract the statistics for each label
        OUTOVL=$WSTAT/${ASHS_SUBJID}_${side}_${segtype}_overlap.txt
        rm -rf $OUTOVL

        # Dump volumes into that file
        for ((ilab = 0; ilab < ${#LABIDS[*]}; ilab++)); do

          # The id of the label
          i=${LABIDS[ilab]};
          SUB=${LABNAMES[ilab]};

          # Extract the overlap
          OVL=$(cat $OVLFILE | awk -v k=$i '$1==k && NR>4 {print $4}')

          if [[ $OVL ]]; then
            echo $ASHS_SUBJID $side $SUB $OVL >> $OUTOVL
          fi

        done

      fi

    fi
  done
done

# Last thing: compute ICV (for now using BET mask)
if [[ -f $ASHS_ATLAS/template/template_bet_mask.nii.gz ]]; then

  # Warp the BET mask
  greedy -d 3 \
    -rm $ASHS_ATLAS/template/template_bet_mask.nii.gz $TMPDIR/icv.nii.gz \
    -rf $ASHS_WORK/mprage.nii.gz \
    -r $ASHS_WORK/affine_t1_to_template/t1_to_template_affine_inv.mat \
       $ASHS_WORK/ants_t1_to_temp/greedy_t1_to_template_invwarp.nii.gz

  # Get T1 voxel volume
  VVOX=$(voxel_size $TMPDIR/icv.nii.gz)

  # Get ICV
  VOLUME=$(c3d $TMPDIR/icv.nii.gz -thresh 0.5 inf 1 0 -dup -lstat | tail -n 1 | awk '{print $6}')
  ICV=$(echo "$VVOX $VOLUME" | awk '{print $1*$2}')

  # Write the volume information to output file
  echo $ASHS_SUBJID $ICV > $WSTAT/${ASHS_SUBJID}_icv.txt

fi

# Report final progress
job_progress 1
