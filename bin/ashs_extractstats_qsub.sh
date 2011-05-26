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

# directory for the subfields (separate for different parameter values)
WSUB=$WORK/subfields
WSTAT=$WORK/final
mkdir -p $WSTAT

# Names of segmentations
for side in left right; do

	SBC=$WSUB/bcfh_heuristic_wgtavg_${side}_native.nii.gz
	if [[ -f $SBC ]]; then

    # Get voxel volume
    DVOX=($($BIN/c3d $SBC -info | cut -f 3 -d ';' | sed -e "s/.*\[//" -e "s/\].*//" -e 's/, / /g'))
    VVOX=$(echo "${DVOX[0]} * ${DVOX[1]} * ${DVOX[2]}" | bc -l);

    # Create an output file
    FNBODYVOL=$WSTAT/${SUBJID}_${side}_volumes.txt 
    rm -rf $FNBODYVOL

    # Dump volumes into that file
    SUB=("bkg" "CA1" "CA2" "DG" "CA3" "HEAD" "TAIL" "misc" "SUB" "ERC" "PHG")
    for i in 1 2 3 4 5 6 8 9 10; do

      # Get the number of slices in this subfield 
      NBODY=$($BIN/c3d $SBC -thresh $i $i 1 0 -trim 0mm -info \
        | grep 'dim = ' | cut -f 1 -d ';' | awk '{print $7;}' | sed -e "s/]//")

      # Get the volume of this subfield
      VOLUME=$($BIN/c3d $SBC -thresh $i $i 1 0 -voxel-sum | awk {'print $3'})

      # Get the volume of this subfield
      VSUB=$(echo "$VVOX * $VOLUME" | /usr/bin/bc -l)

      # Write the volume information to output file
      echo $SUBJID $side ${SUB[i]} $NBODY $VSUB >> $FNBODYVOL

    done

	fi
done

# Last thing: compute ICV (for now using BET mask)
$BIN/ants/WarpImageMultiTransform 3 $ROOT/data/template/template_bet_mask.nii.gz \
  $TMPDIR/icv.nii.gz -R $WORK/mprage.nii.gz \
  -i $WORK/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
  $WORK/ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii

# Get T1 voxel volume
DVOX=($($BIN/c3d $TMPDIR/icv.nii.gz -info | cut -f 3 -d ';' | sed -e "s/.*\[//" -e "s/\].*//" -e 's/, / /g'))
VVOX=$(echo "${DVOX[0]} * ${DVOX[1]} * ${DVOX[2]}" | bc -l);

# Get ICV
VOLUME=$($BIN/c3d $TMPDIR/icv.nii.gz -thresh 0.5 inf 1 0 -dup -lstat | tail -n 1 | awk '{print $6}')
ICV=$(echo "$VVOX * $VOLUME" | /usr/bin/bc -l)

# Write the volume information to output file
echo $SUBJID $ICV > $WSTAT/${SUBJID}_icv.txt
