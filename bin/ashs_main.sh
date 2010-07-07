#!/bin/bash

# TODO:
#   - Check that the data are in the right orientation, or handle 
#     various orientations

function usage()
{
  cat <<-USAGETEXT
		ashs_main: automatic segmentation of hippocampal subfields
		usage:
		  ashs_main [options]
		required options:
		  -g image          Filename of 3D (g)radient echo MRI (MPRAGE, T1w)
		  -f image          Filename of 2D focal (f)ast spin echo MRI (TSE, T2w)
		  -w path           Working/output directory
		  -r image          Manual slice labeling of the right hippocampus (see below)
		  -l image          Manual slice labeling of the left hippocampus (see below)
		optional:
		  -d                Enable debugging
		  -h                Print help
		  -s integer        Run only one stage (see below); also accepts range (e.g. -s 1-3)
		  -N                No overriding of ANTS/FLIRT results. If a result from an earlier run
	                            exists, don't run ANTS/FLIRT again
      -I string         Subject ID (for stats output). Defaults to last word of working dir.
      -q string         List of additional options to pass to qsub (Sun Grid Engine)
      -V                Display version information and exit
		stages:
		  0:                initialize work directory
		  1:                fit to population template
		  2:                multi-atlas registration
		  3:                consensus segmentation using voting
		  4:                apply slice boundaries from -l/-r inputs
		  5:                learning-based bias correction
		  6:                q/a
		manual slice labeling (-r. -l options):
		  It is recommended that you supply images where slices belonging to the
		  body/head/tail of the hippocampus are labeled 1/5/6 respectively. Only
		  one voxel per slice needs to be labeled. 
    data requirements:
      The TSE image slice direction should be z. In other words, the dimension
      of TSE image should be 400x400x30 or something like that, not 400x30x400
	USAGETEXT
}

function vers() 
{
  ASHS_VERSION_SVN=$(cat $ASHS_ROOT/bin/ashs_version.txt)
  echo $ASHS_VERSION_SVN
}

# The ROOT variable must be set
if [[ ! $ASHS_ROOT ]]; then
  echo "ASHS_ROOT is not set. Please set this variable to point to the root ASHS directory"
  exit -1;
fi

# Get the architecture and check ability to run binaries
ARCH=$(uname);
BIN=$ASHS_ROOT/ext/$ARCH/bin
BIN_ANTS=$BIN/ants
BIN_FSL=$BIN/fsl
if [[ ! $($BIN/c3d -version | grep 'Version') ]]; then
  echo "Can not execute command \'$BIN/c3d -version\'. Wrong architecture?"
  exit -1
fi

# Check that the data directory exists
if [[ ! -f $ASHS_ROOT/data/train/train21/tse_native.nii.gz ]]; then
  echo "Data files appear to be missing. Can't locate $ASHS_ROOT/data/train/train21/tse_native.nii.gz"
  exit -1
fi

# Print usage by default
if [[ $# -lt 1 ]]; then
  echo "Try $0 -h for more information."
  exit 2
fi

# Read the options
while getopts "g:f:w:s:r:l:q:I:NdhV" opt; do
  case $opt in

    g) MPRAGE=$OPTARG;;
    f) TSE=$OPTARG;;
    w) WORK=$OPTARG;;
		s) STAGE_SPEC=$OPTARG;;
		N) SKIP_ANTS=1; SKIP_RIGID=1; ;;
		l) SLL=$OPTARG;;
		r) SLR=$OPTARG;;
    I) SUBJID=$OPTARG;;
    q) QOPTS=$OPTARG;;
    d) set -x -e;;
    h) usage; exit 0;;
    V) vers; exit 0;;
    \?) echo "Unknown option $OPTARG"; exit 2;;
    :) echo "Option $OPTARG requires an argument"; exit 2;;

  esac
done

# Make sure all files exist
if [[ ! $MPRAGE || ! -f $MPRAGE ]]; then
	echo "T1-weighted 3D gradient echo MRI (-g) must be specified"
	exit 2;
elif [[ ! $TSE || ! -f $TSE ]]; then
	echo "T2-weighted 2D fast spin echo MRI (-f) must be specified"
	exit 2;
elif [[ ! $WORK ]]; then
	echo "Working/output directory must be specified"
	exit 2;
fi

# Check that the dimensions of the T2 image are right
DIMS=$($BIN/c3d $TSE -info | cut -d ';' -f 1 | sed -e "s/.*\[//" -e "s/\].*//" -e "s/,//g")
if [[ ${DIMS[2]} > ${DIMS[0]} || ${DIMS[2]} > ${DIMS[1]} ]]; then
  echo "The T2-weighted image has wrong dimensions (fails dim[2] < min(dim[0], dim[1])"
  exit -1
fi

# Subject ID set to work dir last work
if [[ ! $SUBJID ]]; then
  SUBJID=$(basename $WORK)
fi

# Create the working directory and the dump directory
mkdir -p $WORK $WORK/dump $WORK/final

# Run the stages of the script
ROOT=$ASHS_ROOT; 
export ROOT BIN WORK SKIP_ANTS SKIP_RIGID SLL SLR SUBJID BIN_ANTS BIN_FSL

# Set the start and end stages
if [[ $STAGE_SPEC && $STAGE_SPEC =~ "^[0-9]*$" ]]; then
  STAGE_START=$STAGE_SPEC
  STAGE_END=$STAGE_SPEC
elif [[ $STAGE_SPEC && $STAGE_SPEC =~ "^([0-9]*)\-([0-9]*)$" ]]; then
  STAGE_START=${BASH_REMATCH[1]}
  STAGE_END=${BASH_REMATCH[2]}
elif [[ ! $STAGE_SPEC ]]; then
  STAGE_START=0
  STAGE_END=15
else
  echo "Wrong stage specification -s $STAGE_SPEC"
  exit -1;
fi

for ((STAGE=$STAGE_START; STAGE<=$STAGE_END; STAGE++)); do

  case $STAGE in 

    0)
    # Initialize Directory
    echo "Running stage 0: initialize work directory"

    # Copy the input images into the working directory
    $BIN/c3d $MPRAGE -type ushort -o $WORK/mprage.nii.gz
    $BIN/c3d $TSE -type ushort -o $WORK/tse.nii.gz ;;

    1) 
    # Template matching
    echo "Running stage 1: normalize to T1 population template"
    qsub $QOPTS -sync y -j y -o $WORK/dump -cwd -V -N "ashs_stg1" $ROOT/bin/ashs_template_qsub.sh ;;

    2) 
    # Multi-atlas matching 
    echo "Running stage 2: normalize to multiple T1/T2 atlases"
    qsub $QOPTS -t 1-64 -sync y -j y -o $WORK/dump -cwd -V -N "ashs_stg2" $ROOT/bin/ashs_multiatlas_qsub.sh ;;

    3) 
    # Voting
    echo "Running stage 3: voting"
    qsub $QOPTS -t 1-2 -sync y -j y -o $WORK/dump -cwd -V -N "ashs_stg3" $ROOT/bin/ashs_voting_qsub.sh ;;

    4) 
    # Apply heuristic
    echo "Running stage 4: apply slice boundaries to segmentation"
    qsub $QOPTS -t 1-2 -sync y -j y -o $WORK/dump -cwd -V -N "ashs_stg4" $ROOT/bin/ashs_heuristic_qsub.sh ;;

    5) 
    # Apply AdaBoost voting (this is a hack!)
    echo "Running stage 5: AdaBoost bias correction"
    qsub $QOPTS -t 1-2 -sync y -j y -o $WORK/dump -cwd -V -N "ashs_stg5" $ROOT/bin/ashs_biascorr_qsub.sh ;;

    6)
    # Final QA
    echo "Running stage 6: Final QA"
    qsub $QOPTS -sync y -j y -o $WORK/dump -cwd -V -N "ashs_stg6" $ROOT/bin/ashs_finalqa_qsub.sh ;;
  
    7) 
    # Statistics & Volumes
    echo "Running stage 7: Statistics and Volumes"
    qsub $QOPTS -sync y -j y -o $WORK/dump -cwd -V -N "ashs_stg7" $ROOT/bin/ashs_extractstats_qsub.sh ;;

    8) 
    # Fit cm-rep models
    echo "Running stage 8: Thickness analysis"
    qsub $QOPTS -t 1-4 -sync y -j y -o $WORK/dump -cwd -V -N "ashs_stg8" $ROOT/bin/ashs_thickness_qsub.sh ;;

  esac  

done
