#!/bin/bash
#$ -S /bin/bash
set -x -e

# Determine training case and side based on the TASK ID
SIDES=(left right)
side=${SIDES[$(((SGE_TASK_ID - 1) % 2))]}
tid=$(printf %02d $(((SGE_TASK_ID - 1)/ 2)))

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_multiatlas_qsub.sh
	Root: ${ROOT?}
	Working directory: ${WORK?}
	PATH: ${PATH?}
	Subjob ID: ${SGE_TASK_ID}
	Side: ${side?}
	Training subject: ${tid?}
BLOCK1

# Existing directory structure
WFSL=$WORK/flirt_t2_to_t1
WANT=$WORK/ants_t1_to_temp

# Training directory and training data
TDIR=$ROOT/data/train/train${tid}
TSEG=$TDIR/segjp_native_${side}.nii.gz

# Create directory for this registration
WREG=$WORK/tseg_${side}_train${tid}
mkdir -p $WREG

# Run ANTS with current image as fixed, training image as moving
if [[ $SKIP_ANTS \
  && -f $WREG/antsregAffine.txt \
  && -f $WREG/antsregWarpxvec.nii.gz \
  && -f $WREG/antsregWarpyvec.nii.gz \
  && -f $WREG/antsregWarpzvec.nii.gz \
  && -f $WREG/antsregInverseWarpxvec.nii.gz \
  && -f $WREG/antsregInverseWarpyvec.nii.gz \
  && -f $WREG/antsregInverseWarpzvec.nii.gz ]]; then

	# If registration exists, skip this step
	echo "Skipping ANTS registration $side/$tid"

else

	$BIN_ANTS/ANTS 3 \
		-x $WORK/tse_to_chunktemp_${side}_regmask.nii.gz \
		-m PR[$WORK/mprage_to_chunktemp_${side}.nii.gz,$TDIR/mprage_to_chunktemp_${side}.nii.gz,0.5,4] \
		-m PR[$WORK/tse_to_chunktemp_${side}.nii.gz,$TDIR/tse_to_chunktemp_${side}.nii.gz,0.5,4] \
		-o $WREG/antsreg.nii.gz \
		-i 120x120x40 -v -t SyN[0.5] --continue-affine false | tee $WREG/ants_output.txt	

fi

# Apply the warp to the moving image(s). This is still working in template space
$BIN_ANTS/WarpImageMultiTransform 3 $TDIR/mprage_to_chunktemp_${side}.nii.gz \
	$WREG/mprage_reslice_to_chunk.nii.gz -R $WORK/mprage_to_chunktemp_${side}.nii.gz \
	$WREG/antsregWarp.nii.gz

$BIN_ANTS/WarpImageMultiTransform 3 $TDIR/tse_to_chunktemp_${side}.nii.gz \
	$WREG/tse_reslice_to_chunk.nii.gz -R $WORK/tse_to_chunktemp_${side}.nii.gz \
	$WREG/antsregWarp.nii.gz

# Apply the warp, but with just a single level of interpolation
$BIN_ANTS/WarpImageMultiTransform 3 $TDIR/tse_native.nii.gz \
	$WREG/tse_reslice_to_chunk_oneinterp.nii.gz	\
	-R $WORK/tse_to_chunktemp_${side}.nii.gz \
	$WREG/antsregWarp.nii.gz \
	$TDIR/ants_t1_to_tempWarp.nii.gz \
	$TDIR/ants_t1_to_tempAffine.txt \
	$TDIR/flirt_t2_to_t1_ITK.txt

# Apply the warps to each of the labels, going directly from native space to native space
mkdir -p $WREG/sf
for((LAB=0;LAB<=10;LAB++)); do

	SFID=$(printf "sf%02i" $LAB)

	$BIN/c3d $TSEG -thresh $LAB $LAB 1 0 -smooth 0.24mm -o $WREG/sf/${SFID}_src.nii.gz

	# TODO: replace training full-warps with already combined warps residing in template space
	$BIN_ANTS/WarpImageMultiTransform 3 $WREG/sf/${SFID}_src.nii.gz \
    $WREG/sf/${SFID}_to_native.nii.gz -R $WORK/tse.nii.gz \
		-i $WORK/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt \
		-i $WORK/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
		$WORK/ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii \
		$WREG/antsregWarp.nii.gz \
		$TDIR/ants_t1_to_tempWarp.nii.gz \
		$TDIR/ants_t1_to_tempAffine.txt \
		$TDIR/flirt_t2_to_t1_ITK.txt

	$BIN_ANTS/WarpImageMultiTransform 3 $WREG/sf/${SFID}_src.nii.gz \
		$WREG/sf/${SFID}_to_chunk.nii.gz -R $WORK/tse_to_chunktemp_${side}.nii.gz \
		$WREG/antsregWarp.nii.gz \
		$TDIR/ants_t1_to_tempWarp.nii.gz \
		$TDIR/ants_t1_to_tempAffine.txt \
		$TDIR/flirt_t2_to_t1_ITK.txt

done

$BIN/c3d $WORK/tse_to_chunktemp_${side}.nii.gz $WREG/tse_reslice_to_chunk_oneinterp.nii.gz \
	-ncc 4x4x4 $WORK/tse_to_chunktemp_${side}_regmask.nii.gz -times \
	-o $WREG/tse_nccmap.nii.gz

