#!/bin/bash
#$ -S /bin/bash
set -x -e

# Determine side based on the TASK ID
SIDES=(left right)
side=${SIDES[$(((SGE_TASK_ID - 1) % 2))]}

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_voting_qsub.sh
	Root: ${ROOT?}
	Working directory: ${WORK?}
	Subjob ID: ${SGE_TASK_ID}
	Side: ${side?}
	PATH: ${PATH?}
BLOCK1

# Create directory for the subfields (separate for different parameter values)
WSUB=$WORK/subfields
mkdir -p $WSUB

# Rank the metric images and create weight maps
ALPHA=1
RANKIN=($(ls $WORK/tseg_${side}_*/tse_nccmap.nii.gz))
RANKOUT=($(echo ${RANKIN[*]} | sed -e "s/tse_nccmap/tse_nccrank/g"))
$BIN/c3d -verbose ${RANKIN[*]} -foreach -replace nan -1 -endfor -rank \
	-foreach -scale -1 -scale $ALPHA -exp -smooth 1.2mm -endfor -oo ${RANKOUT[*]}

# Compute the consensus labeling for each subfield
for((LAB=0;LAB<=10;LAB++)); do

	SFID=$(printf "sf%02i" $LAB)

	# Average the label masks
	$BIN/c3d $WORK/tseg_${side}_*/sf/${SFID}_to_chunk.nii.gz -mean -o $WSUB/${SFID}_to_chunk_${side}_avg.nii.gz

	# Now compute a weighted average of the label masks
	$BIN/c3d -verbose \
		$WORK/tseg_${side}_*/tse_nccrank.nii.gz \
		$WORK/tseg_${side}_*/sf/${SFID}_to_chunk.nii.gz \
		-reorder 0.5 -weighted-sum-voxelwise -o $WSUB/${SFID}_to_chunk_${side}_wgtavg.nii.gz
		
	# Map the label masks into subject space
	for type in avg wgtavg; do
		$BIN_ANTS/WarpImageMultiTransform 3 $WSUB/${SFID}_to_chunk_${side}_${type}.nii.gz \
			$WSUB/${SFID}_to_native_${side}_${type}.nii.gz \
			-R $WORK/tse.nii.gz \
			-i $WORK/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt \
			-i $WORK/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
			$WORK/ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii
	done

done

for kind in avg wgtavg; do

	# Perform voting in chunk template space
	$BIN/c3d $WSUB/sf??_to_chunk_${side}_${kind}.nii.gz -vote \
		-o $WSUB/consensus_${kind}_${side}_chunk.nii.gz
			
	# Peform voting in native space
	$BIN/c3d $WSUB/sf??_to_native_${side}_${kind}.nii.gz -vote \
		-o $WSUB/consensus_${kind}_${side}_native.nii.gz

done 


