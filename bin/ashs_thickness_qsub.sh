#!/bin/bash
#$ -S /bin/bash
set -x -e

# Determine side based on the TASK ID
sid=$(((SGE_TASK_ID - 1) % 2))
SIDES=(left right)
side=${SIDES[$sid]}

# Determine the anatomy based on the TASK ID
aid=$(((SGE_TASK_ID - 1) / 2))
ANATS=(caphg dg)
anat=${ANATS[$aid]}

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_thickness_qsub.sh
	Root: ${ROOT?}
	Working directory: ${WORK?}
	PATH: ${PATH?}
  SUBJECT: ${SUBJID?}
  Side: $side
  Anatomy: $anat
BLOCK1

# Directory where to compute thickness
WTHK=$WORK/thickness/$anat/$side
mkdir -p $WTHK 

# Image on which to compute the thickness map
SBC=$WORK/subfields/bcfh_heuristic_wgtavg_${side}_native.nii.gz
if [[ ! -f $SBC ]]; then
  echo "Missing segmentation result"
  exit;
fi

# First, create appropriate binary image to fit model to
PARAM=("-replace 3 0 5 0 6 0 7 0 -thresh 1 inf 1 0" "-thresh 3 3 1 0")
$BIN/c3d $SBC ${PARAM[$aid]} \
  -int Gaussian 0.2x0.2x1.0mm -trim 10x10x2vox -resample 100x100x400% -thresh 0.5 inf 1 0 \
  -o $WTHK/target.nii.gz

# Generate a mesh for this model
$BIN/vtklevelset $WTHK/target.nii.gz $WTHK/target_mesh.vtk 0.5

# Compute the skeleton and thickness maps 
$BIN/cmrep_vskel -Q $BIN/qvoronoi -p 1.6 -c 1 \
  -I $WTHK/target.nii.gz $WTHK/thickness.nii.gz $WTHK/depthmap.nii.gz \
  $WTHK/target_mesh.vtk $WTHK/skeleton_vtk

# Register target to the binary template so that we can have a thickness map
# in template space. 

# First get a T2-to-template affine transform (from earlier regs)
$BIN/c3d_affine_tool -itk $WORK/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
  -itk $WORK/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt \
  -mult -oitk $WTHK/test_t2_to_tempAffine.txt

# Set the target template
TEMPLATE=$ROOT/data/template/roi/consensus_${side}_${anat}.nii.gz

# Run ANTS (unless previously run)
if [[ -f $WTHK/roi_antsAffine.txt && $SKIP_ANTS ]]; then
    echo "SKIPPING ANTS"
else
  $BIN_ANTS/ANTS 3 -m MSQ[$TEMPLATE,$WTHK/target.nii.gz,1] -o $WTHK/roi_ants.nii.gz \
    -i 100x100x40 -t SyN[0.5] -a $WTHK/test_t2_to_tempAffine.txt --continue-affine FALSE
fi

# Warp the mask to the template space
$BIN_ANTS/WarpImageMultiTransform 3 $WTHK/target.nii.gz \
  $WTHK/target_to_template.nii.gz -R $TEMPLATE \
  $WTHK/roi_antsWarp.nii.gz $WTHK/roi_antsAffine.txt

# Smooth thickness by 1.5 voxels (why this amount?)
$BIN/c3d $WTHK/thickness.nii.gz -smooth 0.6mm $WTHK/target.nii.gz -as M -smooth 0.6mm \
  -reciprocal -replace inf 0 -times -o $WTHK/thickness_sm.nii.gz

# Warp the thickness map into template space
$BIN_ANTS/WarpImageMultiTransform 3 $WTHK/thickness_sm.nii.gz \
  $WTHK/thickness_to_template_intermediate.nii.gz -R $TEMPLATE \
  $WTHK/roi_antsWarp.nii.gz $WTHK/roi_antsAffine.txt

# Divide thickness by the mask, constrain to thickness mask
$BIN/c3d $WTHK/thickness_to_template_intermediate.nii.gz \
  $TEMPLATE -as T -times \
  -o $WORK/final/${SUBJID}_thickness_${side}_${anat}.nii.gz

