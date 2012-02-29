#!/bin/bash
#$ -S /bin/bash
set -x -e

# Read the library
source $ROOT/bin/ashs_lib.sh

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
TEMP_T1_FULL=$ASHS_ATLAS/template/template.nii.gz

# Copy the images into the working directory
if [[ $MPRAGE -nt $WORK/mprage.nii.gz ]]; then
  c3d $MPRAGE -type ushort -o $WORK/mprage.nii.gz
fi

if [[ $TSE -nt $WORK/tse.nii.gz ]]; then
  c3d $TSE -type ushort -o $WORK/tse.nii.gz 
fi


# Histogram match the images to a reference image (used later down the road, but better to do it now)
for kind in tse mprage; do
  c3d $ASHS_ATLAS/ref_hm/ref_${kind}.nii.gz $WORK/${kind}.nii.gz \
    -histmatch 5 -o $WORK/${kind}_histmatch.nii.gz
done

# --- RIGID ALIGNMENT T1/T2 ---
ashs_align_t1t2 $WORK $WFSL

# Use FLIRT to register T1 to template
WAFF=$WORK/affine_t1_to_template
mkdir -p $WAFF

if [[ -f $WAFF/t1_to_template_ITK.txt && $SKIP_RIGID ]]; then

  echo "Skipping Affine Registration"

else

	# Try using FLIRT
  export FSLOUTPUTTYPE=NIFTI_GZ

  # Run flirt with template as reference
  flirt -v -anglerep quaternion -ref $TEMP_T1_FULL -in $WORK/mprage.nii.gz \
    -datatype short -o $WAFF/test_flirt_affine.nii.gz \
    -omat $WAFF/flirt_intermediate_affine.mat -cost corratio -searchcost corratio

  # Convert the transform to ITK
  c3d_affine_tool $WAFF/flirt_intermediate_affine.mat -ref $TEMP_T1_FULL -src $WORK/mprage.nii.gz \
    -fsl2ras -oitk $WAFF/flirt_t1_to_template_ITK.txt

  # Try using ANTS
	ANTS 3 -m MI[$TEMP_T1_FULL,$WORK/mprage.nii.gz,1,32] -o $WAFF/antsaffineonly.nii.gz -i 0 
	WarpImageMultiTransform 3 $WORK/mprage.nii.gz $WAFF/test_ants_affine.nii.gz \
		-R $TEMP_T1_FULL $WAFF/antsaffineonlyAffine.txt

	# Compare the two images
	SIM_FLIRT=$(c3d $TEMP_T1_FULL $WAFF/test_flirt_affine.nii.gz -nmi | awk '{print int(1000*$3)}');
	SIM_ANTS=$(c3d $TEMP_T1_FULL $WAFF/test_ants_affine.nii.gz -nmi | awk '{print int(1000*$3)}');
	if [[ $SIM_FLIRT -gt $SIM_ANTS ]]; then
		cp -a $WAFF/flirt_t1_to_template_ITK.txt $WAFF/t1_to_template_ITK.txt
	else
		cp -a $WAFF/antsaffineonlyAffine.txt $WAFF/t1_to_template_ITK.txt
	fi

fi

# Use ANTS to warp the MPRAGE image to the template
if [[ $SKIP_ANTS \
  && -f $WANT/ants_t1_to_tempAffine.txt \
  && -f $WANT/ants_t1_to_tempWarp.nii \
  && -f $WANT/ants_t1_to_tempInverseWarp.nii ]]; then

    echo "SKIPPING ANTS"
  
else

    ANTS 3 -m PR[$TEMP_T1_FULL,$WORK/mprage.nii.gz,1,4] \
      -x $ROOT/data/template/template_bet_mask.nii.gz \
      -o $WANT/ants_t1_to_temp.nii \
      -a $WAFF/t1_to_template_ITK.txt \
      -i ${ASHS_TEMPLATE_ANTS_ITER} -v | tee $WANT/ants_output.txt
  
fi

# Apply the transformation to the T1 image and to the T1 segmentation
WarpImageMultiTransform 3 $WORK/mprage.nii.gz \
  $WANT/reslice_mprage_to_template.nii.gz -R $TEMP_T1_FULL \
  $WANT/ants_t1_to_tempWarp.nii $WANT/ants_t1_to_tempAffine.txt

# Apply the transformation to the T2 image 
WarpImageMultiTransform 3 $WORK/tse.nii.gz $WANT/reslice_tse_to_template.nii.gz -R $TEMP_T1_FULL \
  $WANT/ants_t1_to_tempWarp.nii $WANT/ants_t1_to_tempAffine.txt $WFSL/flirt_t2_to_t1_ITK.txt

# Transform all of the images into the chunk template space
ashs_reslice_to_template $WORK $ASHS_ATLAS

# We don't need to keep the histmatch images around
if [[ $ASHS_TIDY ]]; then
	rm -rf $WORK/tse_histmatch.nii.gz $WORK/mprage_histmatch.nii.gz
fi
