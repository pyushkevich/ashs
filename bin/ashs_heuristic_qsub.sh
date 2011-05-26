#!/bin/bash
#$ -S /bin/bash
set -x -e

# Determine side based on the TASK ID
sid=$(((SGE_TASK_ID - 1) % 2))
if [[ $sid == 0 ]]; then
	side=left
	SL=$SLL
else
	side=right
	SL=$SLR
fi

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_voting_qsub.sh
	Root: ${ROOT?}
	Working directory: ${WORK?}
	Subjob ID: ${SGE_TASK_ID}
	Side: ${side?}
	Slice Labeling: ${SL}
	PATH: ${PATH?}
BLOCK1

# directory for the subfields (separate for different parameter values)
WSUB=$WORK/subfields

# Quit if slice labeling does not exist
if [[ ! -f $SL ]]; then

  # Apply heuristics to the extents of slices
  $BIN/subfield_mask_slices $WSUB/consensus_wgtavg_${side}_native.nii.gz \
    $WSUB/heuristic_mask_${side}.nii.gz

else

  # Apply heuristics to the extents of slices
  $BIN/subfield_mask_slices $SL $WSUB/heuristic_mask_${side}.nii.gz

fi

# Subfield relabeling
for kind in avg wgtavg; do

	# Perform heuristic voting
	PROB=($(ls $WSUB/sf*_to_native_${side}_${kind}.nii.gz))

	$BIN/c3d -verbose $WSUB/heuristic_mask_${side}.nii.gz -popas M \
	-push M -shift 0 -replace 5 1 6 2 -popas Min \
	-push Min -thresh 2 2 1 0 -popas M_body \
	-push Min -thresh 1 1 1 0 -popas M_head \
	-push Min -thresh 3 3 1 0 -popas M_tail \
	${PROB[0]} \
	${PROB[1]} -push M_body -times \
	${PROB[2]} -push M_body -times \
	${PROB[3]} -push M_body -times \
	${PROB[4]} -push M_body -times \
	${PROB[5]} -push M_head -times \
	${PROB[6]} -push M_tail -times \
	${PROB[7]} -scale 0 \
	${PROB[8]} -push Min -thresh 2 3 1 0 -times \
	${PROB[9]}  -push M -thresh 5 6 1 0 -times \
	${PROB[10]} -push M -thresh 5 6 0 1 -times \
	-vote -o $WSUB/consensus_heuristic_${kind}_${side}_native.nii.gz

done

