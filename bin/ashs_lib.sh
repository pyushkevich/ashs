#!/bin/bash - 

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


# ----------------------------------------------------
# A library of routines shared by all the ASHS scripts
# ----------------------------------------------------

# Source the master script in order to get the PATH variable right
source ${ASHS_ROOT?}/bin/ashs_common_master.sh

# If a custom config file specified, first read the default config file
if [[ ! ${ASHS_CONFIG?} == ${ASHS_ROOT?}/bin/ashs_config.sh ]]; then
  source ${ASHS_ROOT?}/bin/ashs_config.sh
fi

# Read the config file
source ${ASHS_CONFIG?}

# Limit the number of threads to one if using QSUB
if [[ $ASHS_USE_QSUB ]]; then
  export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
fi

# Determine the TMDDIR parameter for the child scripts
function get_tmpdir()
{
  # If TMPDIR is already set (i.e., ashs_main is run from qsub)
  # then this will create a subdirectory in TMPDIR. Otherwise
  # this will create a subdirectory in /tmp
  echo $(mktemp -d -t ashs.XXXXXXXX)
}

# Simulate qsub environment in bash (i.e., create and set tempdir)
# This function expects the name of the job as first parameter
function fake_qsub()
{
  local MYNAME=$1
  shift 1

  local PARENTTMPDIR=$TMPDIR
  TMPDIR=$(get_tmpdir)
  export TMPDIR

  bash $* 2>&1 | tee $ASHS_WORK/dump/${MYNAME}.o$(date +%Y%m%d_%H%M%S)

  rm -rf $TMPDIR
  TMPDIR=$PARENTTMPDIR
}

# Function to run a job in parallels
function gnu_parallel_qsub()
{
  local MYNAME=$1
  shift 1

  local PARENTTMPDIR=$TMPDIR
  TMPDIR=$(get_tmpdir)
  export TMPDIR

  export NSLOTS=1

  echo "Starting parallel job $MYNAME"
  bash $* &> $ASHS_WORK/dump/${MYNAME}.o$(date +%Y%m%d_%H%M%S)
  echo "Parallel job $MYNAME done"

  rm -rf $TMPDIR
  TMPDIR=$PARENTTMPDIR
}


# Submit a job to the queue (or just run job) and wait until it finishes
function qsubmit_sync()
{
  local MYNAME=$1
  shift 1

  if [[ $ASHS_USE_QSUB ]]; then
    qsub $QOPTS -sync y -j y -o $ASHS_WORK/dump -cwd -V -N $MYNAME $*
  else
    fake_qsub $MYNAME $*
  fi
}

# Submit an array of jobs parameterized by a single parameter. The parameter is passed in 
# to the job after all other positional parameters
function qsubmit_single_array()
{
  local NAME=$1
  local PARAM=$2
  shift 2;

  # Generate unique name to prevent clashing with qe
  local UNIQ_NAME=${NAME}_${$}

  # Special handling for GNU parallel
  if [[ $ASHS_USE_PARALLEL ]]; then

    export -f gnu_parallel_qsub
    parallel -u gnu_parallel_qsub ${NAME}_{} $* {} ::: $PARAM

  else

    for p1 in $PARAM; do

      if [[ $ASHS_USE_QSUB ]]; then

        qsub $QOPTS -j y -o $ASHS_WORK/dump -cwd -V -N ${UNIQ_NAME}_${p1} $* $p1

      else

        fake_qsub ${NAME}_${p1} $* $p1

      fi
    done

  fi

  # Wait for the jobs to be done
  qwait "${UNIQ_NAME}_*"
}


# Submit an array of jobs parameterized by a single parameter. The parameter is passed in 
# to the job after all other positional parameters
function qsubmit_double_array()
{
  local NAME=$1
  local PARAM1=$2
  local PARAM2=$3
  shift 3;

  # Generate unique name to prevent clashing with qe
  local UNIQ_NAME=${NAME}_${$}

  if [[ $ASHS_USE_PARALLEL ]]; then

    export -f gnu_parallel_qsub
    export -f get_tmpdir
    parallel -u "gnu_parallel_qsub ${NAME}_{1}_{2} $* {1} {2}" ::: $PARAM1 ::: $PARAM2

  else

    for p1 in $PARAM1; do
      for p2 in $PARAM2; do

        if [[ $ASHS_USE_QSUB ]]; then

          qsub $QOPTS -j y -o $ASHS_WORK/dump -cwd -V -N ${UNIQ_NAME}_${p1}_${p2} $* $p1 $p2

        else

          fake_qsub ${NAME}_${p1}_${p2} $* $p1 $p2

        fi
      done
    done

  fi

  # Wait for the jobs to be done
  qwait "${UNIQ_NAME}_*"
}


      



# Submit an array of jobs to the queue
function qsubmit_array()
{
  local NAME=$1
  local SIZE=$2
  local CMD="$3 $4 $5 $6 $7 $8 $9"

  if [[ $ASHS_USE_QSUB ]]; then
    qsub $QOPTS -t 1-${SIZE} -sync y -j y -o $ASHS_WORK/dump -cwd -V -N $NAME $CMD
  else
    for ((i=1; i<=$SIZE;i++)); do

      local PARENTTMPDIR=$TMPDIR
      TMPDIR=$(get_tmpdir)
      export TMPDIR

      local SGE_TASK_ID=$i
      export SGE_TASK_ID

      bash $CMD 2>&1 | tee $ASHS_WORK/dump/${NAME}.o$(date +%Y%m%d_%H%M%S)

      rm -rf $TMPDIR
      TMPDIR=$PARENTTMPDIR
    done
  fi
}

# Wait for qsub to finish
function qwait() 
{
  if [[ $ASHS_USE_QSUB ]]; then
    qsub -b y -sync y -j y -o /dev/null -cwd -hold_jid "$1" /bin/sleep 1
  fi
}

# Report the version number
function vers() 
{
  ASHS_VERSION_SVN=$(cat $ASHS_ROOT/bin/ashs_version.txt)
  echo $ASHS_VERSION_SVN
}


# This function defines common file variables for ASHS subjects. It is an 
# alternative to using filenames
function ashs_subj_vars()
{
  # Get the directory in which these variables are defined
  local WORK=${1?}
  local WANT=$WORK/ants_t1_to_temp
  local WAFF=$WORK/affine_t1_to_template
  local WFSL=$WORK/flirt_t2_to_t1

  # Main images
  SUBJ_MPRAGE=$WORK/mprage.nii.gz
  SUBJ_TSE=$WORK/tse.nii.gz

  # Transformations between T1 and T2
  SUBJ_AFF_T2T1_MAT=$WFSL/flirt_t2_to_t1.mat
  SUBJ_AFF_T2T1_INVMAT=$WFSL/flirt_t2_to_t1_inv.mat

  # Unit transformations between T1 and template
  SUBJ_AFF_T1TEMP_MAT=$WAFF/t1_to_template_affine.mat
  SUBJ_AFF_T1TEMP_INVMAT=$WAFF/t1_to_template_affine_inv.mat

  SUBJ_T1TEMP_WARP=$WANT/greedy_t1_to_template_warp.nii.gz
  SUBJ_T1TEMP_INVWARP=$WANT/greedy_t1_to_template_invwarp.nii.gz

  # Composite transformations from T1 and template
  SUBJ_T1TEMP_TRAN="$SUBJ_T1TEMP_WARP $SUBJ_AFF_T1TEMP_MAT"
  SUBJ_T1TEMP_INVTRAN="$SUBJ_AFF_T1TEMP_INVMAT $SUBJ_T1TEMP_INVWARP"

  # Composite transformations from T2 and template
  SUBJ_T2TEMP_TRAN="$SUBJ_T1TEMP_WARP $SUBJ_AFF_T1TEMP_MAT $SUBJ_AFF_T2T1_MAT"
  SUBJ_T2TEMP_INVTRAN="$SUBJ_AFF_T2T1_INVMAT $SUBJ_AFF_T1TEMP_INVMAT $SUBJ_T1TEMP_INVWARP"
}

# ASHS side-specific variables
function ashs_subj_side_vars()
{
  local WORK=${1?}
  side=${2?}

  ashs_subj_vars $WORK 

  # Side-specific stuff
  SUBJ_SIDE_TSE_TO_CHUNKTEMP=$WORK/tse_to_chunktemp_${side}.nii.gz
  SUBJ_SIDE_MPRAGE_TO_CHUNKTEMP=$WORK/mprage_to_chunktemp_${side}.nii.gz
  SUBJ_SIDE_TSE_TO_CHUNKTEMP_REGMASK=$WORK/tse_to_chunktemp_${side}_regmask.nii.gz
  SUBJ_SIDE_TSE_NATCHUNK=$WORK/tse_native_chunk_${side}.nii.gz
}

# ASHS atlas-specific variables
function ashs_atlas_vars()
{
  set -e
  local tid=${1?}
  local ATLAS_MODE=${2?}

  # If in atlas mode complain
  if [[ $ATLAS_MODE -eq 1 ]]; then
    TDIR=$ASHS_WORK/atlas/${tid}
    ### ATLAS_T1TEMP_WARP=$TDIR/ants_t1_to_chunktemp_${side}Warp.nii.gz
    ### ATLAS_T1TEMP_INVWARP=$TDIR/ants_t1_to_chunktemp_${side}InverseWarp.nii.gz
  else
    TDIR=$ASHS_ATLAS/train/train${tid} 
  fi

  # Basic images
  ATLAS_MPRAGE=$TDIR/mprage.nii.gz
  ATLAS_TSE=$TDIR/tse.nii.gz

  # Transformations between T1 and T2
  ATLAS_AFF_T2T1_MAT=$TDIR/flirt_t2_to_t1_ITK.txt
  ATLAS_AFF_T2T1_INVMAT="$ATLAS_AFF_T2T1_MAT,-1"

  # Unit transformations between T1 and template
  ATLAS_AFF_T1TEMP_MAT=$TDIR/ants_t1_to_tempAffine.txt
  ATLAS_AFF_T1TEMP_INVMAT="$ATLAS_AFF_T1TEMP_MAT,-1"
}  

# ASHS atlas and side-specific variables
function ashs_atlas_side_vars()
{
  local tid=${1?}
  local side=${2?}
  local ATLAS_MODE=${3?}
  local TDIR=$ASHS_ATLAS/train/train${tid}

  ashs_atlas_vars $tid $ATLAS_MODE

  # To save space, the atlas only stores chunk-wise warps to template
  if [[ $ATLAS_MODE -eq 0 ]]; then
    ATLAS_T1TEMP_WARP=$TDIR/ants_t1_to_chunktemp_${side}Warp.nii.gz
    ATLAS_T1TEMP_INVWARP=$TDIR/ants_t1_to_chunktemp_${side}InverseWarp.nii.gz
  fi

  # Composite transformations from T1 and template
  ATLAS_T1TEMP_TRAN="$ATLAS_T1TEMP_WARP $ATLAS_AFF_T1TEMP_MAT"
  ATLAS_T1TEMP_INVTRAN="$ATLAS_AFF_T1TEMP_INVMAT $ATLAS_T1TEMP_INVWARP"

  # Composite transformations from T2 and template
  ATLAS_T2TEMP_TRAN="$ATLAS_T1TEMP_WARP $ATLAS_AFF_T1TEMP_MAT $ATLAS_AFF_T2T1_MAT"
  ATLAS_T2TEMP_INVTRAN="$ATLAS_AFF_T2T1_INVMAT $ATLAS_AFF_T1TEMP_INVMAT $ATLAS_T1TEMP_INVWARP"

  # Template chunk pieces
  ATLAS_SIDE_TSE_TO_CHUNKTEMP=$TDIR/tse_to_chunktemp_${side}.nii.gz
  ATLAS_SIDE_MPRAGE_TO_CHUNKTEMP=$TDIR/mprage_to_chunktemp_${side}.nii.gz
  ATLAS_SIDE_TSE_TO_CHUNKTEMP_REGMASK=$TDIR/tse_to_chunktemp_${side}_regmask.nii.gz
  ATLAS_SIDE_TSE_NATCHUNK=$TDIR/tse_native_chunk_${side}.nii.gz

  # Atlas segmentations in native space
  ATLAS_SEG=$TDIR/tse_native_chunk_${side}_seg.nii.gz
}


# This function aligns the T1 and T2 images together. It takes two parameters: the 
# path to a directory containing images mprage.niigz and tse.nii.gz, and a path to 
# the directory where the output of the registration is stored.
function ashs_align_t1t2()
{
  local ASHS_WORK=${1?}
  local WFSL=${2?}

  # Load the subject variables
  ashs_subj_vars $ASHS_WORK

  if [[ -f $SUBJ_AFF_T2T1_MAT && $ASHS_SKIP_RIGID ]]; then
    
    echo "Skipping Rigid Registration"

  else

    # Use FLIRT to match T2 to T1
    export FSLOUTPUTTYPE=NIFTI_GZ

    # Make the ASHS_TSE image isotropic and extract a chunk
    c3d $SUBJ_TSE -resample ${ASHS_TSE_ISO_FACTOR?} \
      -region ${ASHS_TSE_ISO_REGION_CROP?} \
			-type short -o $TMPDIR/tse_iso.nii.gz

    # Reslice T1 into space of T2 chunk
    c3d $TMPDIR/tse_iso.nii.gz $SUBJ_MPRAGE \
      -reslice-identity -type short -o $TMPDIR/mprage_to_tse_iso.nii.gz

    # Use greedy affine mode to perform registration 
    greedy -d 3 -a -dof 6 -m MI -n 100x100x10 \
      -i $TMPDIR/tse_iso.nii.gz $TMPDIR/mprage_to_tse_iso.nii.gz \
      -o $SUBJ_AFF_T2T1_MAT 

    # Invert the matrix
    c3d_affine_tool $SUBJ_AFF_T2T1_MAT -inv -o $SUBJ_AFF_T2T1_INVMAT

<<'FLIRT_OLD'
    # Run flirt with T2 as reference (does it matter at this point?)
    flirt -v -ref $TMPDIR/tse_iso.nii.gz -in $TMPDIR/mprage_to_tse_iso.nii.gz \
      -omat $WFSL/flirt_intermediate.mat -cost normmi -dof 6 ${ASHS_FLIRT_MULTIMODAL_OPTS?}

    # Convert the T1-T2 transform to ITK
    c3d_affine_tool $WFSL/flirt_intermediate.mat -ref $TMPDIR/tse_iso.nii.gz \
			-src $TMPDIR/mprage_to_tse_iso.nii.gz \
      -fsl2ras -inv -oitk $WFSL/flirt_t2_to_t1_ITK.txt \
      -o $SUBJ_AFF_T2T1_MAT -inv -o $SUBJ_AFF_T2T1_INVMAT
FLIRT_OLD

  fi
}

# This function performs multi-modality ANTS registration between an atlas and the target image
# See below for the list of variables that should be defined.
# Lastly, there is a parameter, whether this is being run in altas building mode (1 yes, 0 no)
function ashs_ants_pairwise()
{
  # Check required variables
  local side=${1?}
  local tid=${2?}
  local WREG=${3?}
	local ATLAS_MODE=${4?}

  # Get the subject data
  ashs_subj_vars $ASHS_WORK
  ashs_subj_side_vars $ASHS_WORK $side
  ashs_atlas_side_vars $tid $side $ATLAS_MODE

  # Define the local transform variables
  local ATLAS_SUBJ_AFF_MAT=$WREG/greedy_atlas_to_subj_affine.mat
  local ATLAS_SUBJ_WARP=$WREG/greedy_atlas_to_subj_warp.nii.gz

  # Run ANTS with current image as fixed, training image as moving
  if [[ $ASHS_SKIP_ANTS && -f $ATLAS_SUBJ_AFF_MAT && -f $ATLAS_SUBJ_WARP ]]; then

    # If registration exists, skip this step
    echo "Skipping Greedy registration $side/$tid"

  else

    # Are we running multi-component registration
    if [[ $(echo $ASHS_PAIRWISE_ANTS_T1_WEIGHT | awk '{print ($1 == 0.0)}') -eq 1 ]]; then

      # T1 has a zero weight
      local METRIC_TERM="-i $SUBJ_SIDE_TSE_TO_CHUNKTEMP $ATLAS_SIDE_TSE_TO_CHUNKTEMP"
      
    else

      # T1 has non-zero weight
      local T2WGT=$(echo $ASHS_PAIRWISE_ANTS_T1_WEIGHT | awk '{print (1.0 - $1)}')
      local METRIC_TERM=\
        "-w $T2WGT \
         -i $SUBJ_SIDE_TSE_TO_CHUNKTEMP $ATLAS_SIDE_TSE_TO_CHUNKTEMP \
         -w $ASHS_PAIRWISE_ANTS_T1_WEIGHT \
         -i $SUBJ_SIDE_MPRAGE_TO_CHUNKTEMP $ATLAS_SIDE_MPRAGE_TO_CHUNKTEMP"
    fi 

    # Perform greedy affine registration with mask and NCC metric
    time greedy -d 3 -a \
      -gm $SUBJ_SIDE_TSE_TO_CHUNKTEMP_REGMASK $METRIC_TERM -o $ATLAS_SUBJ_AFF_MAT \
      -m NCC 2x2x2 -n 60x60x0 -float 

    # Perform greedy deformable registration with NCC metric
    time greedy -d 3 -it $ATLAS_SUBJ_AFF_MAT \
      -gm $SUBJ_SIDE_TSE_TO_CHUNKTEMP_REGMASK $METRIC_TERM -o $ATLAS_SUBJ_WARP \
      -m NCC 2x2x2 -n 60x60x20 -e 0.5 -float

    # TODO: restore the config stuff
    ### -n $ASHS_PAIRWISE_ANTS_ITER -e $ASHS_PAIRWISE_ANTS_STEPSIZE

    # Old code
    ## ANTS 3 \
    ##  -x tse_to_chunktemp_${side}_regmask.nii.gz $ANTS_METRIC_TERM -o $WREG/antsreg.nii.gz \
		##	-i $ASHS_PAIRWISE_ANTS_ITER -t SyN[$ASHS_PAIRWISE_ANTS_STEPSIZE] -v

  fi

  # Define the resliced images
  local ATLAS_RESLICE=$WREG/atlas_to_native.nii.gz
  local ATLAS_RESLICE_SEG=$WREG/atlas_to_native_segvote.nii.gz

  # Reslice the atlas into target space
  if [[ $ASHS_SKIP && -f $ATLAS_RESLICE && -f $ATLAS_RESLICE_SEG ]]; then

    echo "Skipping reslicing into native space"
  
  else

    # Apply full composite warp from atlas TSE to subject TSE
    greedy -d 3 -rf $SUBJ_SIDE_TSE_NATCHUNK \
      -rm $ATLAS_TSE $ATLAS_RESLICE \
      -ri LABEL ${ASHS_LABEL_SMOOTHING} -rm $ATLAS_SEG $ATLAS_RESLICE_SEG \
      -r $SUBJ_T2TEMP_INVTRAN $ATLAS_SUBJ_WARP $ATLAS_SUBJ_AFF_MAT $ATLAS_T2TEMP_TRAN

  fi

	# In tidy mode, we can clean up after this step
	if [[ $ASHS_TIDY ]]; then
		rm -rf $ATLAS_SUBJ_WARP $ATLAS_SUBJ_AFF_MAT
	fi
}



# This function maps histogram-corrected whole-brain images into the 
# reference space of the template. Parameters are the atlas directory
# containing the input images, and the atlas directory
function ashs_reslice_to_template()
{
  local ASHS_WORK=${1?}
  local ATLAS=${2?}
  local WANT=$ASHS_WORK/ants_t1_to_temp
  local WFSL=$ASHS_WORK/flirt_t2_to_t1

  # Subject files generated in this script
  ashs_subj_vars $ASHS_WORK

  # Apply the transformation to the masks
  local side
  for side in left right; do

    ashs_subj_side_vars $ASHS_WORK $side

    # Only generate these targets if needed
    if [[ $ASHS_SKIP && \
          -f $SUBJ_SIDE_TSE_TO_CHUNKTEMP && \
          -f $SUBJ_SIDE_TSE_TO_CHUNKTEMP_REGMASK && \
          -f $SUBJ_SIDE_MPRAGE_TO_CHUNKTEMP && \
          -f $SUBJ_SIDE_TSE_NATCHUNK ]]
    then
      echo "Skipping reslicing of data to template and native ROI"
    else

      # Define the reference space
      REFSPACE=$ATLAS/template/refspace_${side}.nii.gz

      # Map the TSE image to the template space
      greedy -d 3 -rf $REFSPACE \
        -rm $SUBJ_TSE $SUBJ_SIDE_TSE_TO_CHUNKTEMP \
        -r $SUBJ_T2TEMP_TRAN

      # Map the mprage image to the template space
      greedy -d 3 -rf $REFSPACE \
        -rm $SUBJ_MPRAGE $SUBJ_SIDE_MPRAGE_TO_CHUNKTEMP \
        -r $SUBJ_T1TEMP_TRAN

      # Create a custom mask for the ASHS_TSE image
      c3d $SUBJ_SIDE_TSE_TO_CHUNKTEMP -verbose -pim r -thresh 0.001% inf 1 0 \
        -erode 0 4x4x4 $REFSPACE -times -type uchar \
        -o $SUBJ_SIDE_TSE_TO_CHUNKTEMP_REGMASK 

      # Create a native-space chunk of the ASHS_TSE image 
      greedy -d 3 -rf $SUBJ_TSE \
        -rm $SUBJ_SIDE_TSE_TO_CHUNKTEMP_REGMASK $TMPDIR/natmask.nii.gz \
        -r $SUBJ_T2TEMP_INVTRAN

      # Notice that we pad a little in the z-direction. This is to make sure that 
      # we get all the slices in the image, otherwise there will be problems with 
      # the voting code.
      c3d $TMPDIR/natmask.nii.gz -thresh 0.5 inf 1 0 -trim 0x0x2vox \
        $SUBJ_TSE -reslice-identity -type short -o $SUBJ_SIDE_TSE_NATCHUNK

    fi

    # We also resample the segmentation (if it exists, i.e., training mode)
    if [[ -f $ASHS_WORK/seg_${side}.nii.gz ]]; then
      if [[ $ASHS_SKIP && -f $ASHS_WORK/tse_native_chunk_${side}_seg.nii.gz ]]
      then
        echo "Skipping reslicing the segmentation to native space"
      else
        c3d $ASHS_WORK/tse_native_chunk_${side}.nii.gz $ASHS_WORK/seg_${side}.nii.gz \
          -int 0 -reslice-identity -type short \
          -o $ASHS_WORK/tse_native_chunk_${side}_seg.nii.gz
      fi
    fi

  done
}

# This function call the label fusion command given the set of training images, the id of the 
# subject to segment, the side, and the output filename
function ashs_label_fusion()
{
  id=${1?}
  FNOUT=${2?}
  POSTERIOR=${3?}

cat <<-BLOCK1
	Script: ashs_atlas_pairwise.sh
	Root: ${ASHS_ROOT?}
	Working directory: ${ASHS_WORK?}
	PATH: ${PATH?}
	Subject: ${id?}
	Train Set: ${TRAIN?}
	Output: ${FNOUT?}
	Side: ${side?}
BLOCK1

  # Go to the atlas directory
  cd $ASHS_WORK/atlas/$id

  # Perform label fusion using the atlases. We check for the existence of the atlases
  # so if one or two of the registrations failed, the whole process does not crash
  local ATLASES=""
  local ATLSEGS=""
  for id in $TRAIN; do
    local ACAND=pairwise/tseg_${side}_train${id}/atlas_to_native.nii.gz
    local SCAND=pairwise/tseg_${side}_train${id}/atlas_to_native_segvote.nii.gz
    if [[ -f $ACAND && -f $SCAND ]]; then
      ATLASES="$ATLASES $ACAND"
      ATLSEGS="$ATLSEGS $SCAND"
    fi
  done

  # If there are heuristics, make sure they are supplied to the LF program
  if [[ $ASHS_HEURISTICS ]]; then
    EXCLCMD=$(for fn in $(ls heurex/heurex_${side}_*.nii.gz); do \
      echo "-x $(echo $fn | sed -e "s/.*_//g" | awk '{print 1*$1}') $fn"; \
      done)
  fi

  # Run the label fusion program
  /usr/bin/time -f "Label Fusion: walltime=%E, memory=%M" \
    label_fusion 3 -g $ATLASES -l $ATLSEGS \
      -m $ASHS_MALF_STRATEGY -rp $ASHS_MALF_PATCHRAD -rs $ASHS_MALF_SEARCHRAD \
      $EXCLCMD \
      -p ${POSTERIOR} \
      tse_native_chunk_${side}.nii.gz $FNOUT
}

# This version of the function is called in ashs_main, where there are no
# ground truth segmentations to go on. It requires two rounds of label fusion
function ashs_label_fusion_apply()
{
	cd $ASHS_WORK

	BOOTSTRAP=${1?}

cat <<-BLOCK1
	Script: ashs_atlas_pairwise.sh
	Root: ${ASHS_ROOT?}
	Working directory: ${ASHS_WORK?}
	PATH: ${PATH?}
	Side: ${side?}
	Bootstrap: ${BOOTSTRAP?}
BLOCK1

	if [[ $BOOTSTRAP -eq 1 ]]; then
		TDIR=$ASHS_WORK/bootstrap
	else
		TDIR=$ASHS_WORK/multiatlas
	fi

	mkdir -p $TDIR/fusion

  # Perform label fusion using the atlases
	local ATLASES=$(ls $TDIR/tseg_${side}_train*/atlas_to_native.nii.gz)
	local ATLSEGS=$(ls $TDIR/tseg_${side}_train*/atlas_to_native_segvote.nii.gz)

  # Run the label fusion program
	local RESULT=$TDIR/fusion/lfseg_raw_${side}.nii.gz
  label_fusion 3 -g $ATLASES -l $ATLSEGS \
    -m $ASHS_MALF_STRATEGY -rp $ASHS_MALF_PATCHRAD -rs $ASHS_MALF_SEARCHRAD \
    -p $TDIR/fusion/posterior_${side}_%03d.nii.gz \
    tse_native_chunk_${side}.nii.gz $RESULT

  # If there are heuristics, make sure they are supplied to the LF program
  if [[ $ASHS_HEURISTICS ]]; then

		# Apply the heuristics to generate exclusion maps
    mkdir -p $TDIR/heurex
    subfield_slice_rules $RESULT $ASHS_HEURISTICS $TDIR/heurex/heurex_${side}_%04d.nii.gz

		# Rerun label fusion
    local EXCLCMD=$(for fn in $(ls $TDIR/heurex/heurex_${side}_*.nii.gz); do \
      echo "-x $(echo $fn | sed -e "s/.*_//g" | awk '{print 1*$1}') $fn"; \
      done)

		label_fusion 3 -g $ATLASES -l $ATLSEGS \
			-m $ASHS_MALF_STRATEGY -rp $ASHS_MALF_PATCHRAD -rs $ASHS_MALF_SEARCHRAD \
			$EXCLCMD \
      -p $TDIR/fusion/posterior_${side}_%03d.nii.gz \
			tse_native_chunk_${side}.nii.gz $TDIR/fusion/lfseg_heur_${side}.nii.gz

	else
		# Just make a copy
		cp -a $TDIR/fusion/lfseg_raw_${side}.nii.gz $TDIR/fusion/lfseg_heur_${side}.nii.gz
	fi

	# Perform AdaBoost correction. In addition to outputing the corrected segmentation,
  # we output posterior probabilities for each label. 
  for kind in usegray nogray; do

    # The part of the command that's different for the usegray and nogray modes
    if [[ $kind = 'usegray' ]]; then GRAYCMD="-g tse_native_chunk_${side}.nii.gz"; else GRAYCMD=""; fi

    sa $TDIR/fusion/lfseg_heur_${side}.nii.gz \
      $ASHS_ATLAS/adaboost/${side}/adaboost_${kind} \
      $TDIR/fusion/lfseg_corr_${kind}_${side}.nii.gz $EXCLCMD \
      $GRAYCMD -p $TDIR/fusion/posterior_${side}_%03d.nii.gz \
      -op $TDIR/fusion/posterior_corr_${kind}_${side}_%03d.nii.gz

  done
  
  # If there are reference segs, we have to repeat this again, but with heuristics from 
  # the reference segmentations. This allows us to make a more fair comparison
  if [[ $ASHS_HEURISTICS && -f $ASHS_WORK/refseg/refseg_${side}.nii.gz ]]; then

		# Rerun label fusion
    local EXCLCMD=$(for fn in $(ls $ASHS_WORK/refseg/heurex/heurex_${side}_*.nii.gz); do \
      echo "-x $(echo $fn | sed -e "s/.*_//g" | awk '{print 1*$1}') $fn"; \
      done)

		label_fusion 3 -g $ATLASES -l $ATLSEGS \
			-m $ASHS_MALF_STRATEGY -rp $ASHS_MALF_PATCHRAD -rs $ASHS_MALF_SEARCHRAD \
			$EXCLCMD \
      -p $TDIR/fusion/posterior_vsref_${side}_%03d.nii.gz \
			tse_native_chunk_${side}.nii.gz $TDIR/fusion/lfseg_vsref_heur_${side}.nii.gz

    # Rerun AdaBoost
    for kind in usegray nogray; do

      # The part of the command that's different for the usegray and nogray modes
      if [[ $kind = 'usegray' ]]; then GRAYCMD="-g tse_native_chunk_${side}.nii.gz"; else GRAYCMD=""; fi

      sa $TDIR/fusion/lfseg_vsref_heur_${side}.nii.gz \
        $ASHS_ATLAS/adaboost/${side}/adaboost_${kind} \
        $TDIR/fusion/lfseg_vsref_corr_${kind}_${side}.nii.gz $EXCLCMD \
        $GRAYCMD -p $TDIR/fusion/posterior_vsref_${side}_%03d.nii.gz \
        -op $TDIR/fusion/posterior_corr_${kind}_vsref_${side}_%03d.nii.gz

    done

  fi
}




# *************************************************************
#  ENTRY-LEVEL ATLAS-BUILDING FUNCTIONS (Called by ashs_train)
# *************************************************************

# Initialize the atlas directory for each atlas
function ashs_atlas_initialize_directory()
{
  # Initialize Directory
  for ((i=0;i<$N;i++)); do
    id=${ATLAS_ID[i]}
    qsub $QOPTS -j y -o $ASHS_WORK/dump -cwd -V -N "ashs_atlas_initdir_${id}" \
      $ASHS_ROOT/bin/ashs_atlas_initdir_qsub.sh \
        $id ${ATLAS_T1[i]} ${ATLAS_T2[i]} ${ATLAS_LS[i]} ${ATLAS_RS[i]}
  done

  # Wait for jobs to complete
  qwait "ashs_atlas_initdir_*"
}

# Build the template using SYN
function ashs_atlas_build_template()
{
  # All work is done in the template directory
  mkdir -p $ASHS_WORK/template_build
  pushd $ASHS_WORK/template_build

  # Populate
  CMDLINE=""
  for id in ${ATLAS_ID[*]}; do
    ln -sf $ASHS_WORK/atlas/${id}/mprage.nii.gz ./${id}_mprage.nii.gz
    CMDLINE="$CMDLINE ${id}_mprage.nii.gz"
  done

  # Run the template code
  if [[ -f atlastemplate.nii.gz && $ASHS_SKIP_ANTS ]]; then
    echo "Skipping template building"
  else
    export ANTSPATH=$ASHS_ANTS/
    export ANTS_QSUB_OPTS=$QOPTS
    buildtemplateparallel.sh -d 3 -o atlas -m ${ASHS_TEMPLATE_ANTS_ITER?} -r 1 -t GR -s CC $CMDLINE
    
    # Compress the warps
    for id in ${ATLAS_ID[*]}; do
      shrink_warp 3 atlas${id}_mprageWarp.nii.gz atlas${id}_mprageWarp.nii.gz
      shrink_warp 3 atlas${id}_mprageInverseWarp.nii.gz atlas${id}_mprageInverseWarp.nii.gz
    done

    # Copy the template into the final folder
    mkdir -p $ASHS_WORK/final/template/
    cp -av atlastemplate.nii.gz  $ASHS_WORK/final/template/template.nii.gz

  fi

  # We should now map everyone's segmentation into the template to build a mask
  for side in left right; do

    if [[ $ASHS_SKIP_ANTS && \
          -f $ASHS_WORK/final/template/refspace_meanseg_${side}.nii.gz && \
          -f $ASHS_WORK/final/template/refspace_mprage_${side}.nii.gz ]]
    then
      echo "Skipping template ROI definition"
    else

      # Create a working directory
      mkdir -p mask_${side}
      pushd mask_${side}

      # Warp each segmentation
      CMDLINE=""
      for id in ${ATLAS_ID[*]}; do

        WarpImageMultiTransform 3 $ASHS_WORK/atlas/${id}/seg_${side}.nii.gz \
          ${id}_seg_${side}.nii.gz -R ../atlastemplate.nii.gz \
          ../atlas${id}_mprageWarp.nii.gz ../atlas${id}_mprageAffine.txt \
          $ASHS_WORK/atlas/${id}/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt --use-NN

        CMDLINE="$CMDLINE ${id}_seg_${side}.nii.gz"
      done

      # Average the segmentations and create a target ROI with desired resolution
      c3d $CMDLINE \
        -foreach -thresh 0.5 inf 1 0 -endfor -mean -as M -thresh $ASHS_TEMPLATE_MASK_THRESHOLD inf 1 0 \
        -o meanseg_${side}.nii.gz -dilate 1 ${ASHS_TEMPLATE_ROI_DILATION} \
        -trim ${ASHS_TEMPLATE_ROI_MARGIN?} \
        -resample-mm ${ASHS_TEMPLATE_TARGET_RESOLUTION?} \
        -o refspace_${side}.nii.gz \
        ../atlastemplate.nii.gz -reslice-identity -o refspace_mprage_${side}.nii.gz \
        -push M -reslice-identity -thresh 0.5 inf 1 0 -o refspace_meanseg_${side}.nii.gz

      cp -a refspace*${side}.nii.gz $ASHS_WORK/final/template/

      popd

    fi

  done

  # Use the first image in the atlas set as the reference for histogram matching
  if [[ $ASHS_SKIP && \
        -f $ASHS_WORK/final/ref_hm/ref_mprage.nii.gz && \
        -f $ASHS_WORK/final/ref_hm/ref_tse.nii.gz ]]
  then
    echo "Skipping defining the histogram matching references"
  else
    mkdir -p $ASHS_WORK/final/ref_hm
    HMID=${ATLAS_ID[${ASHS_TARGET_ATLAS_FOR_HISTMATCH?}]}
    cp -av $ASHS_WORK/atlas/$HMID/mprage.nii.gz $ASHS_WORK/final/ref_hm/ref_mprage.nii.gz
    cp -av $ASHS_WORK/atlas/$HMID/tse.nii.gz $ASHS_WORK/final/ref_hm/ref_tse.nii.gz
  fi

  popd
}

# Resample all atlas images to the template
function ashs_atlas_resample_tse_to_template()
{
  # Now resample each atlas to the template ROI chunk, setting up for n-squared 
  # registration.
  for ((i=0;i<$N;i++)); do
    id=${ATLAS_ID[i]}
    qsub $QOPTS -j y -o $ASHS_WORK/dump -cwd -V -N "ashs_atlas_resample_${id}" \
      $ASHS_ROOT/bin/ashs_atlas_resample_to_template_qsub.sh $id
  done

  # Wait for jobs to complete
  qwait "ashs_atlas_resample_*"
}

# This is a high-level routine called from the atlas code to run the pairwise registration
function ashs_atlas_register_to_rest()
{
  # Launch all the individual registrations
  for id in ${ATLAS_ID[*]}; do
    for side in left right; do
      for tid in ${ATLAS_ID[*]}; do

        if [[ $id != $tid ]]; then

          qsub $QOPTS -j y -o $ASHS_WORK/dump -cwd -V -N "ashs_nsq_${id}_${tid}" \
            $ASHS_ROOT/bin/ashs_atlas_pairwise_qsub.sh $id $tid $side

        fi
      done
    done
  done

  # Wait for the jobs to be done
  qwait "ashs_nsq_*"
}


# Organize the output directory
function ashs_atlas_organize_final()
{
  # The final directory
  FINAL=$ASHS_WORK/final

  # Generate a file that makes this an ASHS atlas
  cat > $FINAL/ashs_atlas_vars.sh <<-CONTENT
		ASHS_ATLAS_VERSION=$(ashs_version)
		ASHS_ATLAS_N=$N
	CONTENT

  # Copy the heuristic rules if they exist
  if [[ $ASHS_HEURISTICS ]]; then
    cp -av $ASHS_HEURISTICS $FINAL/ashs_heuristics.txt
  else
    rm -rf $ASHS_HEURISTICS $FINAL/ashs_heuristics.txt
  fi

  # Generate directories for each of the training images
  for ((i=0;i<$N;i++)); do

    # Use codes for atlas names, to make this publishable online
    CODE=$(printf train%03d $i)
    id=${ATLAS_ID[i]}

		# The output directory for this atlas
    ODIR=$FINAL/train/$CODE
    IDIR=$ASHS_WORK/atlas/$id
    mkdir -p $ODIR

		# Copy the full ASHS_TSE (fix this later!)
		cp -av $IDIR/tse.nii.gz $ODIR

    # Copy the stuff we need into the atlas directory
    for side in left right; do

      # Copy the images 
      cp -av \
        $IDIR/tse_native_chunk_${side}.nii.gz \
        $IDIR/tse_native_chunk_${side}_seg.nii.gz \
        $IDIR/tse_to_chunktemp_${side}.nii.gz \
        $IDIR/tse_to_chunktemp_${side}_regmask.nii.gz \
        $IDIR/mprage_to_chunktemp_${side}.nii.gz \
				$IDIR/seg_${side}.nii.gz \
        $ODIR/

      # Copy the transformation to the template space. We only care about
      # the part of this transformation that involves the template, to we
      # can save a little space here. 
      c3d $IDIR/mprage_to_chunktemp_${side}.nii.gz -popas REF \
        -mcs $IDIR/ants_t1_to_temp/ants_t1_to_tempWarp.nii.gz \
        -foreach -insert REF 1 -reslice-identity -info -endfor \
        -omc $ODIR/ants_t1_to_chunktemp_${side}Warp.nii.gz 

    done

    # Copy the affine transforms
    cp -av \
      $IDIR/ants_t1_to_temp/ants_t1_to_tempAffine.txt \
      $IDIR/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt \
      $ODIR/

  done

  # Copy the SNAP segmentation labels
  mkdir -p $FINAL/snap
  cp -av $ASHS_LABELFILE $FINAL/snap/snaplabels.txt

  # If a custom config file specified, put a copy of it into the atlas directory
  if [[ ! ${ASHS_CONFIG} == ${ASHS_ROOT}/bin/ashs_config.sh ]]; then
    cp -av $ASHS_CONFIG $FINAL/ashs_user_config.sh
  else
    rm -rf $FINAL/ashs_user_config.sh
    touch $FINAL/ashs_user_config.sh
  fi

  # Copy the system config as well
  cp -av $ASHS_ROOT/bin/ashs_config.sh $FINAL/ashs_system_config.sh

  # Copy the adaboost training files
  for side in left right; do
    mkdir -p $FINAL/adaboost/${side}
    cp -av $ASHS_WORK/train/main_${side}/adaboost* $FINAL/adaboost/${side}/
  done

  # Generate a brain mask for the template
  export FSLOUTPUTTYPE=NIFTI_GZ
  cd $FINAL/template
  bet2 template.nii.gz template_bet -m -v 
}


# Organize cross-validation directories for running cross-validation experiments
function ashs_atlas_organize_xval()
{
  # Training needs to be performed separately for the cross-validation experiments and for the actual atlas
  # building. We will do this all in one go to simplify the code. We create BASH arrays for the train/test
  local NXVAL=0
  if [[ -f $ASHS_XVAL ]]; then
    NXVAL=$(cat $ASHS_XVAL | wc -l)
  fi

  # Arrays for the main training
  XV_TRAIN[0]=${ATLAS_ID[*]}
  XV_TEST[0]=""
  XVID[0]="main"

  # Arrays for the cross-validation training
  for ((jx=1;jx<=$NXVAL;jx++)); do

    XV_TEST[$jx]=$(cat $ASHS_XVAL | awk "NR==$jx { print \$0 }")
    XV_TRAIN[$jx]=$(echo $( for((i=0;i<$N;i++)); do echo ${ATLAS_ID[i]}; done | awk "\"{${XV_TEST[$jx]}}\" !~ \$1 {print \$1}"))
    XVID[$jx]=xval$(printf %04d $jx)

  done

  # Organize the directories
  for ((i=1; i<=$NXVAL; i++)); do

    # The x-validation ID
    local XID=${XVID[i]}
    local TRAIN=${XV_TRAIN[i]}

    # Create the atlas for this experiment (based on all the training cases for it)
    local XVATLAS=$ASHS_WORK/xval/${XID}/atlas
    mkdir -p $XVATLAS

    # Create links to stuff from the main atlas
    for fn in $(ls $ASHS_WORK/final | grep -v "\(adaboost\|train\)"); do
      ln -sf $ASHS_WORK/final/$fn $XVATLAS/$fn
    done

    # For the train directory, only choose the atlases that are in the training set
    mkdir -p $XVATLAS/train
    local myidx=0
    for tid in $TRAIN; do

      # Get the index of the training ID among all training subjects
      local srcdir=$(for qid in ${ATLAS_ID[*]}; do echo $qid; done | awk "\$1~/$tid/ {printf \"train%03d\n\",NR-1}")
      local trgdir=$(printf "train%03d" $myidx)

      # Link the two directories
      ln -sf $ASHS_WORK/final/train/$srcdir $XVATLAS/train/$trgdir

      # Increment the index
      myidx=$((myidx+1))
    done

    # For the adaboost directory, link the corresponding cross-validation experiment
    mkdir -p $XVATLAS/adaboost
    for side in left right; do
      ln -sf $ASHS_WORK/train/${XID}_${side} $XVATLAS/adaboost/$side
    done

    # Now, for each test subject, initialize the ASHS directory with links
    for testid in ${XV_TEST[i]}; do

      # Create the directory for this run
      local IDFULL=${XID}_test_${testid}
      local XVSUBJ=$ASHS_WORK/xval/${XID}/test/$IDFULL
      mkdir -p $XVSUBJ

      # The corresponding atlas directory
      local MYATL=$ASHS_WORK/atlas/$testid/

      # Populate the critical results to avoid having to run registrations twice
      mkdir -p $XVSUBJ/affine_t1_to_template
      ln -sf $MYATL/ants_t1_to_temp/ants_t1_to_tempAffine.txt $XVSUBJ/affine_t1_to_template/t1_to_template_ITK.txt

      mkdir -p $XVSUBJ/ants_t1_to_temp
      ln -sf $MYATL/ants_t1_to_temp/ants_t1_to_tempAffine.txt $XVSUBJ/ants_t1_to_temp/ants_t1_to_tempAffine.txt
      ln -sf $MYATL/ants_t1_to_temp/ants_t1_to_tempWarp.nii.gz $XVSUBJ/ants_t1_to_temp/ants_t1_to_tempWarp.nii.gz
      ln -sf $MYATL/ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii.gz $XVSUBJ/ants_t1_to_temp/ants_t1_to_tempInverseWarp.nii.gz

      mkdir -p $XVSUBJ/flirt_t2_to_t1
      ln -sf $MYATL/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt $XVSUBJ/flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt

      # We can also reuse the atlas-to-target stuff
      myidx=0
      for tid in $TRAIN; do
        for side in left right; do

          local tdir=$XVSUBJ/multiatlas/tseg_${side}_train$(printf %03d $myidx)  
          mkdir -p $tdir

          local sdir=$MYATL/pairwise/tseg_${side}_train${tid}

          for fname in antsregAffine.txt antsregInverseWarp.nii.gz antsregWarp.nii.gz; do
            ln -sf $sdir/$fname $tdir/$fname
          done
        done
        myidx=$((myidx+1))
      done

      # Now, we can launch the ashs_main segmentation for this subject!
      qsub $QOPTS -j y -o $ASHS_WORK/dump -cwd -V -N "ashs_xval_${IDFULL}" \
        $ASHS_ROOT/bin/ashs_main.sh \
          -a $XVATLAS -g $MYATL/mprage.nii.gz -f $MYATL/tse.nii.gz -w $XVSUBJ -N -d \
          -r "$MYATL/seg_left.nii.gz $MYATL/seg_right.nii.gz"

    done

  done

  # Wait for it all to be done
  qwait "ashs_xval_*"
}

# This function Qsubs the bl command unless the output is already present
function ashs_check_bl_result()
{
  local RESFILE=$1;

  if [[ -f $RESFILE ]]; then

    local LASTIT=$(cat $RESFILE | tail -n 1 | awk '{print $1}');

    if [[ $LASTIT -eq $ASHS_EC_ITERATIONS ]]; then
      echo 1
      return
    fi
  fi
}


# This function runs the training for a single label in a training directory
function ashs_bl_train_qsub()
{
  local POSTLIST label id mode xvid side GRAYLIST
  xvid=${1?}
  side=${2?}
  label=${3?}

  # Generate the posterior list for this label
  POSTLIST=$(printf posteriorlist_%03d.txt $label)
  rm -rf $POSTLIST
  for id in $(cat trainids.txt); do
    echo $(printf loo_posterior_${id}_${side}_%03d.nii.gz $label) >> $POSTLIST
  done

  # We will generate two sets of training data: one with the grayscale information and
  # another that only uses the posterior maps, and no grayscale information
  for mode in usegray nogray; do

    # The extra parameter to bl
    if [[ $mode = "usegray" ]]; then GRAYLIST="-g graylist.txt"; else GRAYLIST=""; fi

    # Check if the training results already exist
    if [[ $ASHS_SKIP && $(ashs_check_bl_result adaboost_${mode}-AdaBoostResults-Tlabel${label}) -eq 1 ]]; then

      echo "Skipping training for ${xvid}_${side} label $label adaboost_${mode}"

    else

      # Sampling factor used for dry runs
      local DRYFRAC=0.01

      # Execute a dry run, where we only check the number of samples at the maximum sampling rate
      NSAM=$(bl truthlist.txt autolist.txt $label \
        $ASHS_EC_DILATION $ASHS_EC_PATCH_RADIUS $DRYFRAC 0 adaboost_dryrun \
        $GRAYLIST -p $POSTLIST \
          | grep 'of training data' | tail -n 1 | awk -F '[: ]' '{print $13}')

      # Compute the fraction needed to obtain the desired number of samples per image
      local FRAC=$(echo 0 | awk "{ k=$NSAM / ($ASHS_EC_TARGET_SAMPLES * $DRYFRAC); p=(k==int(k) ? k : 1 + int(k)); print p < 1 ? 1 : 1 / p }")

      # Now run for real
      /usr/bin/time -f "BiasLearn $mode: walltime=%E, memory=%M" \
        bl truthlist.txt autolist.txt $label \
          $ASHS_EC_DILATION $ASHS_EC_PATCH_RADIUS $FRAC $ASHS_EC_ITERATIONS adaboost_${mode} \
            $GRAYLIST -p $POSTLIST

    fi

  done
}

# Top level code for AdaBoost training and cross-validation
function ashs_atlas_adaboost_train()
{
  # Training needs to be performed separately for the cross-validation experiments and for the actual atlas
  # building. We will do this all in one go to simplify the code. We create BASH arrays for the train/test
  local NXVAL=0
  if [[ -f $ASHS_XVAL ]]; then 
    NXVAL=$(cat $ASHS_XVAL | wc -l)
  fi

  # Arrays for the main training
  XV_TRAIN[0]=${ATLAS_ID[*]}
  XV_TEST[0]=""
  XVID[0]="main"

  # Arrays for the cross-validation training
  for ((jx=1;jx<=$NXVAL;jx++)); do

    XV_TEST[$jx]=$(cat $ASHS_XVAL | awk "NR==$jx { print \$0 }")
    XV_TRAIN[$jx]=$(echo $( for((i=0;i<$N;i++)); do echo ${ATLAS_ID[i]}; done | awk "\"{${XV_TEST[$jx]}}\" !~ \$1 {print \$1}"))
    XVID[$jx]=xval$(printf %04d $jx)

  done

  # Perform initial multi-atlas segmentation
  for ((i=0; i<=$NXVAL; i++)); do
    for side in left right; do

      WTRAIN=$ASHS_WORK/train/${XVID[i]}_${side}
      mkdir -p $WTRAIN

      # Perform the segmentation in a leave-one-out fashion among the training images
      for id in ${XV_TRAIN[i]}; do

        # The atlas set for the cross-validation
        TRAIN=$(echo ${XV_TRAIN[i]} | sed -e "s/\<$id\>//g")

        # The output path for the segmentation result
        FNOUT=$WTRAIN/loo_seg_${id}_${side}.nii.gz

        # The output path for the posteriors
        POSTERIOR=$WTRAIN/loo_posterior_${id}_${side}_%03d.nii.gz

        # Perform the multi-atlas segmentation if the outputs exist
        if [[ $ASHS_SKIP && -f $FNOUT ]]; then
          echo "Skipping label fusion for ${XVID[i]}_${side}_loo_${id}"
        else
          qsub $QOPTS -j y -o $ASHS_WORK/dump -cwd -V -N "ashs_lf_${XVID[i]}_${side}_loo_${id}" \
             $ASHS_ROOT/bin/ashs_atlas_lf_qsub.sh $id "$TRAIN" $FNOUT $side $POSTERIOR
        fi
      done
    done
  done

  # Wait for all the segmentations to finish  
  qwait "ashs_lf_*_loo_*"

  # Now perform the training
  for ((i=0; i<=$NXVAL; i++)); do
    for side in left right; do

      WTRAIN=$ASHS_WORK/train/${XVID[i]}_${side}
      pushd $WTRAIN

      # Create text files for input to bl
      rm -rf graylist.txt autolist.txt truthlist.txt trainids.txt
      for id in ${XV_TRAIN[i]}; do
        echo $ASHS_WORK/atlas/$id/tse_native_chunk_${side}.nii.gz >> graylist.txt
        echo $ASHS_WORK/atlas/$id/tse_native_chunk_${side}_seg.nii.gz >> truthlist.txt
        echo loo_seg_${id}_${side}.nii.gz >> autolist.txt
        echo $id >> trainids.txt
      done

      # Count the unique labels in the dataset. Note that for label 0 we perform the dilation to estimate
      # the actual number of background voxels
      for fn in $(cat truthlist.txt); do 
        c3d $fn -dup -lstat; 
      done | awk '$1 ~ /[0-9]+/ && $1 > 0 { h[$1]+=$6; h[0]+=$6 } END {for (q in h) print q,h[q] }' | sort -n > counts.txt

      # For each label, launch the training job. The number of samples is scaled to have roughly 1000 per label
      for label in $(cat counts.txt | awk '{print $1}'); do

        qsub $QOPTS $ASHS_EC_QSUB_EXTRA_OPTIONS -j y -o $ASHS_WORK/dump -cwd -V -N "ashs_bl_${XVID[i]}_${side}_${label}" \
           $ASHS_ROOT/bin/ashs_atlas_bl_qsub.sh ${XVID[i]} $side $label

      done

      popd

    done
  done

  # Wait for the training to finish  
  qwait "ashs_bl_*"

}

# This function checks whether all ASHS outputs have been created for a given stage
function ashs_check_train()
{
  local STAGE=$1
  local NERR=0
  local NWARN=0

  pushd $ASHS_WORK

  for id in ${ATLAS_ID[*]}; do

    if [[ $STAGE -ge 0 ]]; then
      for image in tse.nii.gz mprage.nii.gz flirt_t2_to_t1/flirt_t2_to_t1_ITK.txt \
                   seg_left.nii.gz seg_right.nii.gz
      do
        if [[ ! -f atlas/$id/$image ]]; then
          echo "STAGE $STAGE: missing file atlas/$id/$image"
          let NERR=NERR+1
        fi
      done
    fi

    if [[ $STAGE -ge 1 ]]; then
      for kind in Warp InverseWarp; do
        if [[ ! -f template_build/atlas${id}_mprage${kind}.nii.gz ]]; then
          echo "STAGE $STAGE: missing file template_build/atlas${id}_mprage${kind}.nii.gz"
          let NERR=NERR+1
        fi
      done
    fi

    if [[ $STAGE -ge 2 ]]; then
      for side in left right; do
        for kind in \
          tse_native_chunk_${side} mprage_to_chunktemp_${side} \
          tse_to_chunktemp_${side} tse_to_chunktemp_${side}_regmask \
          seg_${side} 
        do
          if [[ ! -f atlas/$id/${kind}.nii.gz ]]; then
            echo "STAGE $STAGE: missing file atlas/$id/${kind}.nii.gz"
            let NERR=NERR+1
          fi
        done
      done
    fi

    if [[ $STAGE -ge 3 ]]; then
      for side in left right; do
        missing=0
        for tid in ${ATLAS_ID[*]}; do
          if [[ $id != $tid ]]; then
            base=atlas/$id/pairwise/tseg_${side}_train${tid}
            if [[ ! -f $base/atlas_to_native_segvote.nii.gz ]]; then
              echo "STAGE $STAGE: missing file $base/atlas_to_native_segvote.nii.gz"
              let NWARN=NWARN+1
              let missing=missing+1
            fi
          fi
        done
        if [[ $missing -ge $((N-1)) ]]; then
          echo "STAGE $STAGE: missing all pairwise results for atlas $id side $side"
          let NERR=NERR+1
        fi
      done
    fi

  done

  # Cross-validation
  local NXVAL=0
  if [[ -f $ASHS_XVAL ]]; then 
    NXVAL=$(cat $ASHS_XVAL | wc -l)
  fi

  # Arrays for the main training
  XV_TRAIN[0]=${ATLAS_ID[*]}
  XV_TEST[0]=""
  XVID[0]="main"

  # Arrays for the cross-validation training
  for ((jx=1;jx<=$NXVAL;jx++)); do

    XV_TEST[$jx]=$(cat $ASHS_XVAL | awk "NR==$jx { print \$0 }")
    XV_TRAIN[$jx]=$(echo $( for((i=0;i<$N;i++)); do echo ${ATLAS_ID[i]}; done | awk "\"{${XV_TEST[$jx]}}\" !~ \$1 {print \$1}"))
    XVID[$jx]=xval$(printf %04d $jx)

  done

  # Loop over the cross-validation experiments
  if [[ $STAGE -ge 4 ]]; then

    for ((i=0; i<=$NXVAL; i++)); do
      for side in left right; do

        WTRAIN=$ASHS_WORK/train/${XVID[i]}_${side}

        # Perform the segmentation in a leave-one-out fashion among the training images
        for id in ${XV_TRAIN[i]}; do

          if [[ ! -f $WTRAIN/loo_seg_${id}_${side}.nii.gz ]]; then
            echo "STAGE $STAGE: missing file $WTRAIN/loo_seg_${id}_${side}.nii.gz"
            let NERR=NERR+1
          fi

        done

        # Look for the adaboost results
        if [[ ! -f $WTRAIN/counts.txt ]]; then
          echo "STAGE $STAGE: missing file $WTRAIN/counts.txt"
          let NERR=NERR+1
        fi

        for label in $(cat $WTRAIN/counts.txt | awk '{print $1}'); do
          for kind in usegray nogray; do
            if [[ ! -f $WTRAIN/adaboost_${kind}-AdaBoostResults-Tlabel${label} ]]; then
              echo "STAGE $STAGE: missing file $WTRAIN/adaboost_${kind}-AdaBoostResults-Tlabel${label}"
              let NERR=NERR+1
            fi
          done
        done

      done
    done 

  fi

  popd

  echo "*****************************"
  echo "VALIDITY CHECK AT STAGE $STAGE FOUND $NERR ERRORS AND $NWARN WARNINGS"
  echo "*****************************"
  if [[ $NERR -gt 0 ]]; then exit -1; fi
}
