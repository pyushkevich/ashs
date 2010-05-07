#!/bin/bash
#$ -S /bin/bash
set -x -e

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
if [[ -f $WFSL/flirt_t2_to_t1_ITK.txt && $SKIP_RIGID ]]; then
  
  echo "Skipping Rigid Registration"

else

  # Use FLIRT to match T2 to T1
  export FSLOUTPUTTYPE=NIFTI_GZ

  # Make the TSE image isotropic and extract a chunk
  $BIN/c3d $WORK/tse.nii.gz -resample 100x100x500% -region 20x20x0% 60x60x100% -o $WFSL/tse_iso.nii.gz

  # Reslice T1 into space of T2 chunk
  $BIN/c3d $WFSL/tse_iso.nii.gz $WORK/mprage.nii.gz -reslice-identity -o $WFSL/mprage_to_tse_iso.nii.gz

  # Run flirt with T2 as reference (does it matter at this point?)
  $BIN_FSL/flirt -v -ref $WFSL/tse_iso.nii.gz -in $WFSL/mprage_to_tse_iso.nii.gz -o $WFSL/test_flirt.nii.gz \
    -omat $WFSL/flirt_intermediate.mat -cost normmi -dof 6 \
    -searchrx -5 5 -searchry -5 5 -searchrz -5 5 -coarsesearch 3 -finesearch 1 -searchcost normmi

  # Convert the T1-T2 transform to ITK
  $BIN/c3d_affine_tool $WFSL/flirt_intermediate.mat -ref $WFSL/tse_iso.nii.gz -src $WFSL/mprage_to_tse_iso.nii.gz \
    -fsl2ras -inv -oitk $WFSL/flirt_t2_to_t1_ITK.txt

fi

# Use ANTS to warp the MPRAGE image to the template
if [[ -f $WANT/ants_t1_to_tempAffine.txt && $SKIP_ANTS ]]; then

    echo "SKIPPING ANTS"
  
else

    $BIN_ANTS/ANTS 3 -m PR[$TEMP_T1_FULL,$WORK/mprage.nii.gz,1,4] \
      -x $ROOT/data/template/template_bet_mask.nii.gz \
      -o $WANT/ants_t1_to_temp.nii \
      -i 200x120x40 -v -t SyN[0.5] | tee $WANT/ants_output.txt
  
fi

# Apply the transformation to the T1 image and to the T1 segmentation
$BIN_ANTS/WarpImageMultiTransform 3 $WORK/mprage.nii.gz \
  $WANT/reslice_mprage_to_template.nii.gz -R $TEMP_T1_FULL \
  $WANT/ants_t1_to_tempWarp.nii $WANT/ants_t1_to_tempAffine.txt

# Apply the transformation to the T2 image 
$BIN_ANTS/WarpImageMultiTransform 3 $WORK/tse.nii.gz $WANT/reslice_tse_to_template.nii.gz -R $TEMP_T1_FULL \
  $WANT/ants_t1_to_tempWarp.nii $WANT/ants_t1_to_tempAffine.txt $WFSL/flirt_t2_to_t1_ITK.txt

# Apply the transformation to the masks
for side in left right; do

  # Define the reference space
  REFSPACE=$ROOT/data/template/refspace_${side}.nii.gz

  # Map the image to the target space
  $BIN_ANTS/WarpImageMultiTransform 3 $WORK/tse_histmatch.nii.gz \
    $WORK/tse_to_chunktemp_${side}.nii.gz -R $REFSPACE \
    $WANT/ants_t1_to_tempWarp.nii $WANT/ants_t1_to_tempAffine.txt $WFSL/flirt_t2_to_t1_ITK.txt

  # Map the image to the target space
  $BIN_ANTS/WarpImageMultiTransform 3 $WORK/mprage_histmatch.nii.gz \
    $WORK/mprage_to_chunktemp_${side}.nii.gz -R $REFSPACE \
    $WANT/ants_t1_to_tempWarp.nii $WANT/ants_t1_to_tempAffine.txt 

  # Create a custom mask for the TSE image
  $BIN/c3d $WORK/tse_to_chunktemp_${side}.nii.gz -verbose -pim r -thresh 0.001% inf 1 0 \
    -erode 0 4x4x4 $REFSPACE -times -type uchar -o $WORK/tse_to_chunktemp_${side}_regmask.nii.gz

  # Create a combined warp from chunk template to T2 native space
  $BIN_ANTS/ComposeMultiTransform 3 $WANT/ants_t2_to_temp_fullWarp.nii -R $REFSPACE \
    $WANT/ants_t1_to_tempWarp.nii $WANT/ants_t1_to_tempAffine.txt $WFSL/flirt_t2_to_t1_ITK.txt

done

