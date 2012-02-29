#!/bin/bash
#$ -S /bin/bash
set -x -e

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_extractstats_qsub.sh
	Root: ${ROOT?}
	Working directory: ${WORK?}
	PATH: ${PATH?}
  SUBJECT: ${SUBJID?}
BLOCK1

function voxel_size()
{
	echo $(c3d $1 -info-full | grep Spacing | sed -e "s/[^0-9\.]/ /g" | awk '{print $1*$2*$3}')
}

# directory for the subfields (separate for different parameter values)
WSTAT=$WORK/final
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
for segtype in corr heur; do
  for side in left right; do

    SBC=$WORK/bootstrap/fusion/lfseg_${segtype}_${side}.nii.gz
    if [[ -f $SBC ]]; then

      # Get voxel volume
      VVOX=$(voxel_size $SBC)

      # Create an output file
      FNBODYVOL=$WSTAT/${SUBJID}_${side}_${segtype}_volumes.txt 
      rm -rf $FNBODYVOL

      # Dump volumes into that file
      SUB=("bkg" "CA1" "CA2" "DG" "CA3" "HEAD" "TAIL" "misc" "SUB" "ERC" "PHG")
      for ((ilab = 0; ilab < ${#LABIDS[*]}; ilab++)); do

        # The id of the label
        i=${LABIDS[ilab]};
        SUB=${LABNAMES[ilab]};

        # Get the number of slices in this subfield 
        NBODY=$(c3d $SBC -thresh $i $i 1 0 -trim 0mm -info \
          | grep 'dim = ' | cut -f 1 -d ';' | awk '{print $7;}' | sed -e "s/]//")

        # Get the volume of this subfield
        VOLUME=$(c3d $SBC -thresh $i $i 1 0 -voxel-sum | awk {'print $3'})

        # Get the volume of this subfield
        VSUB=$(echo "$VVOX $VOLUME" | awk '{print $1*$2}')

        # Write the volume information to output file
        echo $SUBJID $side $SUB $NBODY $VSUB >> $FNBODYVOL

      done

    fi
  done
done

# Last thing: compute ICV (for now using BET mask)
if [[ -f $ASHS_ATLAS/template/template_bet_mask.nii.gz ]]; then

  # Warp the BET mask
  WarpImageMultiTransform 3 $ASHS_ATLAS/template/template_bet_mask.nii.gz \
    $TMPDIR/icv.nii.gz -R $WORK/mprage.nii.gz \
    -i $WORK/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
    $WORK/ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii

  # Get T1 voxel volume
  VVOX=$(voxel_size $TMPDIR/icv.nii.gz)

  # Get ICV
  VOLUME=$(c3d $TMPDIR/icv.nii.gz -thresh 0.5 inf 1 0 -dup -lstat | tail -n 1 | awk '{print $6}')
  ICV=$(echo "$VVOX * $VOLUME" | awk '{print $1*$2}')

  # Write the volume information to output file
  echo $SUBJID $ICV > $WSTAT/${SUBJID}_icv.txt

fi
