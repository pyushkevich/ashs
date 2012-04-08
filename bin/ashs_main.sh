#!/bin/bash
#$ -S /bin/bash

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
		  -a dir            Location of the atlas directory. Can be a full pathname or a
		                    relative directory name under ASHS_ROOT/data directory. 
		  -g image          Filename of 3D (g)radient echo MRI (MPRAGE, T1w)
		  -f image          Filename of 2D focal (f)ast spin echo MRI (TSE, T2w)
		  -w path           Working/output directory

		optional:
		  -d                Enable debugging
		  -h                Print help
		  -s integer        Run only one stage (see below); also accepts range (e.g. -s 1-3)
		  -N                No overriding of ANTS/FLIRT results. If a result from an earlier run
		                    exists, don't run ANTS/FLIRT again
		  -T                Tidy mode. Cleans up files once they are unneeded. The -N option will
		                    have no effect in tidy mode, because ANTS/FLIRT results will be erased.
		  -I string         Subject ID (for stats output). Defaults to last word of working dir.
		  -V                Display version information and exit
		  -C file           Configuration file. If not passed, uses $ASHS_ROOT/bin/ashs_config.sh
		  -Q                Use Sun Grid Engine (SGE) to schedule sub-tasks in each stage. By default,
		                    the whole ashs_main job runs in a single process. If you are doing a lot
		                    of segmentations and have SGE, it is better to run each segmentation 
		                    (ashs_main) in a separate SGE job, rather than use the -q flag. The -q flag
		                    is best for when you have only a few segmentations and want them to run fast.
      -q OPTS           Pass in additional options to SGE's qsub. Also enables -Q option above.

		stages:
		  1:                fit to population template
		  2:                multi-atlas registration
		  3:                consensus segmentation using voting
		  4:                bootstrap registration
		  5:                bootstrap segmentation using voting
		  6:                segmentation Q/A
		  7:                volumes and statistics

		notes:
		  The TSE image slice direction should be z. In other words, the dimension
		  of TSE image should be 400x400x30 or something like that, not 400x30x400
	USAGETEXT
}

# Run parent-level code
source ${ASHS_ROOT?}/bin/ashs_common_master.sh

# Print usage by default
if [[ $# -lt 1 ]]; then
  echo "Try $0 -h for more information."
  exit 2
fi

# Configuration file
ASHS_CONFIG=$ASHS_ROOT/bin/ashs_config.sh

# Load the library
source $ASHS_ROOT/bin/ashs_lib.sh

# Read the options
while getopts "g:f:w:s:a:q:I:C:NTdhVQ" opt; do
  case $opt in

    a) ATLAS=$OPTARG;;
    g) MPRAGE=$OPTARG;;
    f) TSE=$OPTARG;;
    w) WORK=$OPTARG;;
		s) STAGE_SPEC=$OPTARG;;
		N) SKIP_ANTS=1; SKIP_RIGID=1; ;;
		T) ASHS_TIDY=1;;
    I) SUBJID=$OPTARG;;
    Q) ASHS_USE_QSUB=1;;
    q) ASHS_USE_QSUB=1; QOPTS=$OPTARG;;
    C) ASHS_CONFIG=$OPTARG;;
    d) set -x -e;;
    h) usage; exit 0;;
    V) vers; exit 0;;
    \?) echo "Unknown option $OPTARG"; exit 2;;
    :) echo "Option $OPTARG requires an argument"; exit 2;;

  esac
done

# Check if the required parameters were passed in
echo "Atlas    : ${ATLAS?    "Directory for atlas was not specified. See $0 -h"}
echo "T1 Image : ${MPRAGE?   "T1-weighted MRI was not specified. See $0 -h"}
echo "T2 Image : ${TSE?      "T2-weighted MRI was not specified. See $0 -h"}
echo "WorkDir  : ${WORK?     "Working directory was not specified. See $0 -h"}

# Whether we are using QSUB
if [[ $ASHS_USE_QSUB ]]; then
  if [[ ! $SGE_ROOT ]]; then
    echo "-Q flag used, but SGE is not present."
    exit -1;
  fi
  echo "Using SGE with root $SGE_ROOT and options $QOPTS"
else
  echo "Not using SGE"
fi

# Convert the work directory to absolute path
mkdir -p ${WORK?}
WORK=$(cd $WORK; pwd)
if [[ ! -d $WORK ]]; then 
  echo "Work directory $WORK does not exist";
fi

# Check the atlas location
if [[ -f $ATLAS/ashs_atlas_vars.sh ]]; then
  ASHS_ATLAS=$ATLAS;
elif [[ -f $ASHS_ROOT/data/$ATLAS/ashs_atlas_vars.sh ]]; then
  ASHS_ATLAS=$ASHS_ROOT/data/$ATLAS
else
  echo "Atlas directory must be specified"
  exit 2;
fi

# Check the heuristics in the atlas
if [[ -f $ASHS_ATLAS/ashs_heuristics.txt ]]; then
	ASHS_HEURISTICS=$ASHS_ATLAS/ashs_heuristics.txt
fi

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
DIMS=$(c3d $TSE -info | cut -d ';' -f 1 | sed -e "s/.*\[//" -e "s/\].*//" -e "s/,//g")
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

# Read the configuration file
source $ASHS_CONFIG

# Run the stages of the script
ROOT=$ASHS_ROOT; 
export ROOT PATH WORK SKIP_ANTS SKIP_RIGID SUBJID BIN_ANTS BIN_FSL ASHS_CONFIG ASHS_ATLAS
export ASHS_HEURISTICS ASHS_TIDY MPRAGE TSE

# Set the start and end stages
if [[ $STAGE_SPEC && $STAGE_SPEC =~ "^[0-9]*$" ]]; then
  STAGE_START=$STAGE_SPEC
  STAGE_END=$STAGE_SPEC
elif [[ $STAGE_SPEC && $STAGE_SPEC =~ "^([0-9]*)\-([0-9]*)$" ]]; then
  STAGE_START=${BASH_REMATCH[1]}
  STAGE_END=${BASH_REMATCH[2]}
elif [[ ! $STAGE_SPEC ]]; then
  STAGE_START=1
  STAGE_END=15
else
  echo "Wrong stage specification -s $STAGE_SPEC"
  exit -1;
fi

# Get the number of atlases, other information
source $ASHS_ATLAS/ashs_atlas_vars.sh

for ((STAGE=$STAGE_START; STAGE<=$STAGE_END; STAGE++)); do

  case $STAGE in 

    1) 
    # Template matching
    echo "Running stage 1: normalize to T1 population template"
    qsubmit_sync "ashs_stg1" $ROOT/bin/ashs_template_qsub.sh ;;

    2) 
    # Multi-atlas matching 
    echo "Running stage 2: normalize to multiple T1/T2 atlases"
    qsubmit_array "ashs_stg2" $((ASHS_ATLAS_N*2)) $ROOT/bin/ashs_multiatlas_qsub.sh ;;

    3) 
    # Voting
    echo "Running stage 3: Label Fusion"
    qsubmit_array "ashs_stg3" 2 $ROOT/bin/ashs_voting_qsub.sh 0;;

		4)
		# Bootstrapping
		echo "Running stage 4: Bootstrap segmentation"
    qsubmit_array "ashs_stg3" $((ASHS_ATLAS_N*2)) $ROOT/bin/ashs_bootstrap_qsub.sh ;;

	  5)
		# Bootstrap voting
		echo "Running stage 5: Bootstrap label fusion" 
    qsubmit_array "ashs_stg5" 2 $ROOT/bin/ashs_voting_qsub.sh 1;;

    6)
    # Final QA
    echo "Running stage 6: Final QA"
    qsubmit_sync "ashs_stg6" $ROOT/bin/ashs_finalqa_qsub.sh ;;
  
    7) 
    # Statistics & Volumes
    echo "Running stage 7: Statistics and Volumes"
    qsubmit_sync "ashs_stg7" $ROOT/bin/ashs_extractstats_qsub.sh ;;

    ### 8) 
    ### Fit cm-rep models
    ### echo "Running stage 8: Thickness analysis"
    ### qsubmit_array "ashs_stg8" 4 $ROOT/bin/ashs_thickness_qsub.sh ;;

    ## 4)
    ## # Final QA
    ### echo "Running stage 4: Final QA"
    ### qsub $QOPTS -sync y -j y -o $WORK/dump -cwd -V -N "ashs_stg6" $ROOT/bin/ashs_finalqa_qsub.sh ;;

  esac  

done
