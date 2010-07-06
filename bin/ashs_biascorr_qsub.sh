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
	Script: ashs_biascorr_qsub.sh
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
if [[ ! -f $WSUB/consensus_heuristic_wgtavg_${side}_native.nii.gz ]]; then
	echo "Segmentation from MASV not available. Quitting"
	exit 0;
fi

# Call Hongzhi's BC code
$BIN/bc $WORK/tse.nii.gz $WSUB/consensus_heuristic_wgtavg_${side}_native.nii.gz \
  $sid $ROOT/data/adaboost $WSUB/bcfh_wgtavg_${side}_native.nii.gz

# Now apply subfield remapping (maintain heuristics)
$BIN/subfield_leveler $SL $WSUB/bcfh_wgtavg_${side}_native.nii.gz $WSUB/bcfh_heuristic_wgtavg_${side}_native.nii.gz

# Copy files into the 'final' directory
c3d $WSUB/bcfh_heuristic_wgtavg_${side}_native.nii.gz -type ushort \
  -o $WORK/final/${SUBJID}_${side}_subfields_final.nii.gz

