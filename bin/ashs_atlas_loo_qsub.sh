#!/bin/bash
#$ -S /bin/bash
set -x -e

source ashs_lib.sh

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_atlas_pairwise.sh
	Root: ${ROOT?}
	Working directory: ${WORK?}
	PATH: ${PATH?}
	Subject: ${id?}
  Train Set: ${TRAIN?}
  Crossval Id: ${XID?}
	Side: ${side?}
BLOCK1

# Create output directory
WXVAL=$WORK/xval/${XID}_${side}
mkdir -p $WXVAL

# Outputs 
SEGUNC=$WXVAL/${id}_${side}_xval_seg_uncorr.nii.gz
SEGCOR=$WXVAL/${id}_${side}_xval_seg_corr.nii.gz
SEGTRUTH=$WORK/atlas/$id/tse_native_chunk_${side}_seg.nii.gz

# Call the label fusion library command
ashs_label_fusion $id $SEGUNC

# Apply the correction
pushd $WORK/train/${XID}_${side}

# If there are heuristics, make sure they are supplied to the LF program
if [[ $ASHS_HEURISTICS ]]; then
  pushd $WORK/atlas/$id
  EXCLCMD=$(for fn in $(ls $WORK/atlas/$id/heurex/heurex_${side}_*.nii.gz); do \
    echo "-x $(echo $fn | sed -e "s/.*_//g" | awk '{print 1*$1}') $fn"; \
    done)
  popd
fi

sa $WORK/atlas/$id/tse_native_chunk_${side}.nii.gz $SEGUNC adaboost $SEGCOR $EXCLCMD

popd

# Compute overlaps
LSET=($(c3d $SEGTRUTH -dup -lstat | awk 'NR > 1 {print $1}'))

# Compute overlaps 
for type in corr uncorr; do
  for lab in ${LSET[*]}; do
    OVL=$(c3d $WXVAL/${id}_${side}_xval_seg_${type}.nii.gz $SEGTRUTH -overlap $lab | awk '$1=="OVL:" {print 1.0 * $6}')
    echo $idtest $side $ixval $lab $OVL
  done > $WXVAL/${id}_${side}_xval_ovl_${type}.txt

  # For both corrected and uncorrected, resample the segmentation to native
  # space, full size.
  c3d $WORK/atlas/$id/tse.nii.gz $WXVAL/${id}_${side}_xval_seg_${type}.nii.gz \
    -int 0 -reslice-identity  -o $WXVAL/${id}_${side}_xval_seg_${type}_fullsize.nii.gz

done
