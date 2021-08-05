#!/bin/bash

#######################################################################
#
#  Program:   ASHS (Automatic Segmentation of Hippocampal Subfields)
#  Module:    $Id$
#  Language:  BASH Shell Script
#  Copyright (c) 2012 Paul A. Yushkevich, University of Pennsylvania
#  
#  This file is part of ASHS
#
#  ASHS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details. 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################


# -------------------------------------
# ASHS configuration and parameter file
# -------------------------------------
#
# This file contains default settings for ASHS. Users can make a copy of this
# file and change parameters to better suit their data. To use this file, make
# a copy of the $ASHS_ROOT/bin/ashs_config.sh, edit the copy, and pass the
# copy to the ashs_main or ashs_train scripts using the -C flag. 

# -------------------------------------
# ASHS configuration and parameter file
# -------------------------------------
#
# The hemispheres in which ASHS segmentation is performed. By default this is 
# "left right" but can be set to "left" or "right" if you only need one side.
# The user must still supply the left and right segmentations to the training 
# in the manifest file
ASHS_SIDES="left right"

# -------------------------------------
# ASHS preprocessing parameters
# -------------------------------------
#
# ASHS provide options for non-local mean denoising and super resolution
# upsampling (code provided by Dr. Jose Manjón http://personales.upv.es/jmanjon/).
# Currently, the non-local mean upsampling only support postive interger upsample factors.
# Upsample factors are specified in "lef-right (LR), anterior-posterior (AP), and superior-inferior (SI)" directions in the physical space rather than in x, y, z coordinates
ASHS_MPRAGE_DENOISE=0
ASHS_MPRAGE_SRUPSAMPLE=0
ASHS_MPRAGE_SRUPSAMPLE_FACTOR="1 1 1"
ASHS_TSE_DENOISE=0
ASHS_TSE_SRUPSAMPLE=0
ASHS_TSE_SRUPSAMPLE_FACTOR="1 1 1"

# ---------------------------------
# ASHS_TSE resolution-related parameters
# ---------------------------------
#
# The default parameters are configured for ASHS_TSE images with 0.4 x 0.4 x 2.0
# voxel resolution, with a field of view that extends past the skull in the
# oblique coronal plane. If your data is different, you may want to change 
# these parameters

# The resampling factor to make data isotropic. The data does not have to be
# exactly isotropic. I suggest keeping all numbers multiples of 100. 
ASHS_TSE_ISO_FACTOR="100x100x500%"

# The cropping applied to the ASHS_TSE volume before rigid registration to the T1
# volume. If the field of view in the ASHS_TSE images is limited, you may want to
# change this. Default crops 20% of the image in the oblique coronal plane
ASHS_TSE_ISO_REGION_CROP="20x20x0% 60x60x100%"

# ----------------------------------
# ASHS_TSE-ASHS_MPRAGE registration parameters
# ----------------------------------

# The search parameters for FLIRT. You may want to play with these if the
# rigid alignment of T1 and T2 is failing. You can always override the results
# for any given image manually by performing the registration yourself and
# calling ASHS with the -N flag (to not rerun existing registrations)
ASHS_FLIRT_MULTIMODAL_OPTS="-searchrx -5 5 -searchry -5 5 -searchrz -5 5 -coarsesearch 3 -finesearch 1 -searchcost normmi"

# These are the flags for the additional affine registration between modalities.
# The purpose of additional affine registration is to take care the residule misalignment
# of each side (mainly for MTL segmentation). This is by default disabled.
ASHS_CHUNK_AFFINE=0
ASHS_SKIP_CHUNK_AFFINE=1

# ------------------------------------------------
# ASHS_MPRAGE template creation/registration parameters
# ------------------------------------------------
#
# These parameters affect registration to the template space in ashs_main and
# template construction in ashs_train

# The number of iterations for ANTS when registering ASHS_MPRAGE to the template.
# This is one of the main parameters affecting the runtime of the program,
# since this is the only whole-brain registration performed by ASHS. See ANTS
# documentation for the meaning of this.
ASHS_TEMPLATE_ANTS_ITER="60x20x0"

# Number of stages for template building, total and for rigid registration
ASHS_TEMPLATE_STAGES_TOTAL=4
ASHS_TEMPLATE_STAGES_RIGID=1

# The amount of dilation applied to the average hippocampus mask in order to
# create a registration mask.
ASHS_TEMPLATE_ROI_DILATION="10x10x10vox"

# The size of the margin applied to the above registration mask when creating
# the ROI-specific template. This option only affects ashs_train. You can
# specify this in vox or mm units.
ASHS_TEMPLATE_ROI_MARGIN="4x4x4vox"

# The target resolution of the template (only when building a template with
# ashs_train). This is specified in units of mm. The target resolulion should
# roughly match the in-plane resolution of the ASHS_TSE data, but there is no hard
# rule for this. Setting it too high will slow everything down.
ASHS_TEMPLATE_TARGET_RESOLUTION="0.4688x0.4688x0.4688mm"

# When creating a mask to define the template, we average the input segmentations
# and set a threshold. When there are lots of atlases, 0.5 is reasonable, but for
# smaller atlas sets, this should be reduced
ASHS_TEMPLATE_MASK_THRESHOLD=0.5

# -----------------------------
# Histogram matching parameters
# -----------------------------
#

# The atlas (index into the list of atlases) used as the target for histogram
# matching. This is only used in ashs_train, and the selected atlas is copied
# into the ref_hm directory. By default, the first atlas is used
ASHS_TARGET_ATLAS_FOR_HISTMATCH=0

# The number of control points for histogram matching. See c3d -histmatch
# command. This is used for atlas building and application. In some cases 
# (for example 7T data with large intensity range) histogram matching seems
# to fail. You can set ASHS_HISTMATCH_CONTROLS to 0 turn off matching.
ASHS_HISTMATCH_CONTROLS=5


# -----------------------------------------------
# Pairwise multi-modality registration parameters
# -----------------------------------------------

# The number of Greedy iterations for running pairwise registration
#   - These parameters are for voxel size of 0.4x0.4x2.0 and should be
#   - changed when voxel size is larger and the ROI is similar size. For
#   - whole hippocampus T1 MRI for example, suggest 60x60
ASHS_PAIRWISE_AFFINE_ITER="60x60x0"
ASHS_PAIRWISE_DEFORM_ITER="60x60x20"

# The radius of the cross-correlation patch for pairwise registration
ASHS_PAIRWISE_CROSSCORR_RADIUS="2x2x2"

# The step size for ANTS
ASHS_PAIRWISE_ANTS_STEPSIZE="0.25"

# The relative weight given to the T1 image in the pairwise registration
# Setting this to 0 makes the registration only use the ASHS_TSE images, which
# may be the best option. This has not been tested extensively. Should be 
# a floating point number between 0 and 1
ASHS_PAIRWISE_ANTS_T1_WEIGHT=0

# The amount of smoothing applied to label images before warping them,
# can either be in millimeters (mm) or voxels (vox)
ASHS_LABEL_SMOOTHING="0.24mm"

# -----------------------------------------------
# Label Fusion Parameters
# -----------------------------------------------

# Label fusion strategy. This is the -m parameter to the label_fusion command
# I would stick to the defaults
ASHS_MALF_STRATEGY="Joint[0.1,2]"

# Patch size (specified as radius) for patch search on the MALF method. This is
# the size of the patch used to compute similarity between voxels.
ASHS_MALF_PATCHRAD="3x3x1"

# Neighborhood size (specified as radius) for patch search. This is the range
# where we search for matching patches
ASHS_MALF_SEARCHRAD="3x3x1"

# Padding for the T2 images - this is useful to set to something like 0x0x2 for
# situations where you expect segmentation to go all the way to the very last
# slice of the T2 image. This adds a modest cost to the segmentation
ASHS_MALF_PADDING="0x0x2"


# -----------------------------------------------
# AdaBoost Bias Correction Parameters
# -----------------------------------------------

# Dilation radius. This is the radius applied to the automatic segmentation 
# result to determine the ROI where error correction takes place. We assume
# that all the errors made by the automatic method are in this ROI. 
ASHS_EC_DILATION=1

# Target number of samples for AdaBoost training. If there are more samples,
# only a random fraction of the samples will be used for training. Lower
# values speed up training
ASHS_EC_TARGET_SAMPLES=100000

# The sampling ratio used during the AdaBoost training dry run to determine 
# the total number of samples. This should only be changed if your segmentations
# are large and you are having problems with running out of memory (bad_alloc)
# in the adaboost training stage. In that case, lower this value
ASHS_EC_DRYRUN_SAMPLE_RATE="0.01"

# Number of iterations for AdaBoost training.
ASHS_EC_ITERATIONS=500

# Size of the neighborhood used to derive features for AdaBoost.
ASHS_EC_PATCH_RADIUS=6x6x0

# Minimum number of samples to do bias correction
ASHS_EC_MINIMUM_SAMPLES=100

# Additional QSUB options for error correction training stage. Error correction
# is memory intensive, and some systems have limits on the memory allowed to a
# single qsub job. These qsub options can be used to run jobs in a parallel environment
# An example value of this option might be "-pe mpi 2"
ASHS_EC_QSUB_EXTRA_OPTIONS=
