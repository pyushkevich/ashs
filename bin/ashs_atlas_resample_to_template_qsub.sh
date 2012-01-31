#!/bin/bash
#$ -S /bin/bash
set -x -e

# Include the common file
source ashs_lib.sh

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_atlas_resample_to_template.sh
	Root: ${ROOT?}
	Working directory: ${WORK?}
	PATH: ${PATH?}
  Id: ${id?}
BLOCK1

# Go to my directory
pushd $WORK/atlas/$id

<<'SKIP'

# Copy the warps from the template building directory to the atlas directory
mkdir -p ants_t1_to_temp
for what in Warp.nii.gz InverseWarp.nii.gz Affine.txt; do
  cp -av $WORK/template_build/atlas${id}_mprage${what} ants_t1_to_temp/ants_t1_to_temp${what}
done

# Histogram match the images to a reference image (used later down the road, but better to do it now)
c3d $WORK/final/ref_hm/ref_tse.nii.gz tse.nii.gz \
  -histmatch ${ASHS_HISTMATCH_CONTROLS} -o tse_histmatch.nii.gz

c3d $WORK/final/ref_hm/ref_mprage.nii.gz mprage.nii.gz \
  -histmatch ${ASHS_HISTMATCH_CONTROLS} -o mprage_histmatch.nii.gz

SKIP

# Reslice the atlas to the template
ashs_reslice_to_template . $WORK/final

# If heuristics specified, create exclusion maps for the labels
if [[ $ASHS_HEURISTICS ]]; then

  for side in left right; do

    mkdir -p heurex
    subfield_slice_rules tse_native_chunk_${side}_seg.nii.gz $ASHS_HEURISTICS heurex/heurex_${side}_%04d.nii.gz

  done
fi


popd
