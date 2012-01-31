#!/bin/bash
#$ -S /bin/bash
set -x -e

source $ASHS_CONFIG
source ashs_lib.sh

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_template_qsub.sh
	Root: ${ROOT?}
	Working directory: ${WORK?}
	PATH: ${PATH?}
  Skip Rigid: ${SKIP_RIGID}
  Skip Ants: ${SKIP_ANTS}
BLOCK1

# Ensure directory structure
WFSL=$WORK/flirt_t2_to_t1
WANT=$WORK/ants_t1_to_temp
mkdir -p $WANT $WFSL

# Set some variables
TEMP_T1_FULL=$ROOT/data/template/template.nii.gz

# Histogram match the images to a reference image (used later down the road, but better to do it now)
$BIN/c3d $ROOT/data/ref_hm/ref_tse.nii.gz $WORK/tse.nii.gz -histmatch 5 -o $WORK/tse_histmatch.nii.gz
$BIN/c3d $ROOT/data/ref_hm/ref_mprage.nii.gz $WORK/mprage.nii.gz -histmatch 5 -o $WORK/mprage_histmatch.nii.gz

# --- RIGID ALIGNMENT T1/T2 ---
ashs_align_t1t2 $WORK $WFSL

# Use FLIRT to register T1 to template
if [[ -f $WFSL/flirt_t1_to_template_ITK.txt && $SKIP_RIGID ]]; then

  echo "Skipping Affine Registration"

else

  # Use FLIRT to match T2 to T1
  export FSLOUTPUTTYPE=NIFTI_GZ

  # Run flirt with template as reference
  $BIN_FSL/flirt -v -anglerep quaternion -ref $TEMP_T1_FULL -in $WORK/mprage.nii.gz \
    -o $WFSL/test_flirt_affine.nii.gz \
    -omat $WFSL/flirt_intermediate_affine.mat -cost corratio -searchcost corratio

  # Convert the transform to ITK
  $BIN/c3d_affine_tool $WFSL/flirt_intermediate_affine.mat -ref $TEMP_T1_FULL -src $WORK/mprage.nii.gz \
    -fsl2ras -oitk $WFSL/flirt_t1_to_template_ITK.txt

fi

# Use ANTS to warp the MPRAGE image to the template
if [[ $SKIP_ANTS \
  && -f $WANT/ants_t1_to_tempAffine.txt \
  && -f $WANT/ants_t1_to_tempWarpxvec.nii \
  && -f $WANT/ants_t1_to_tempWarpyvec.nii \
  && -f $WANT/ants_t1_to_tempWarpzvec.nii \
  && -f $WANT/ants_t1_to_tempInverseWarpxvec.nii \
  && -f $WANT/ants_t1_to_tempInverseWarpyvec.nii \
  && -f $WANT/ants_t1_to_tempInverseWarpzvec.nii ]]; then

    echo "SKIPPING ANTS"
  
else

    $BIN_ANTS/ANTS 3 -m PR[$TEMP_T1_FULL,$WORK/mprage.nii.gz,1,4] \
      -x $ROOT/data/template/template_bet_mask.nii.gz \
      -o $WANT/ants_t1_to_temp.nii \
      -a $WFSL/flirt_t1_to_template_ITK.txt \
      -i 200x120x40 -v -t SyN[0.5] | tee $WANT/ants_output.txt
  
fi

# Apply the transformation to the T1 image and to the T1 segmentation
$BIN_ANTS/WarpImageMultiTransform 3 $WORK/mprage.nii.gz \
  $WANT/reslice_mprage_to_template.nii.gz -R $TEMP_T1_FULL \
  $WANT/ants_t1_to_tempWarp.nii $WANT/ants_t1_to_tempAffine.txt

# Apply the transformation to the T2 image 
$BIN_ANTS/WarpImageMultiTransform 3 $WORK/tse.nii.gz $WANT/reslice_tse_to_template.nii.gz -R $TEMP_T1_FULL \
  $WANT/ants_t1_to_tempWarp.nii $WANT/ants_t1_to_tempAffine.txt $WFSL/flirt_t2_to_t1_ITK.txt

# Transform all of the images into the chunk template space
ashs_reslice_to_template $WORK $ROOT/data
