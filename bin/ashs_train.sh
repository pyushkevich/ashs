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


function usage()
{
  cat <<-USAGETEXT
		ashs_train: generate a new training set for ASHS segmentation
		usage:
		  ashs_train [options]
		required options:
		  -D file           Input data file. See below.
		  -L file           Label description file. See below.
		  -w path           Working/output directory
		optional:
		  -d                Enable debugging
		  -h                Print help
		  -s integer        Run only one stage (see below); also accepts range (e.g. -s 1-3)
		  -r file           Apply heuristic rules (see below) to segmentations produced by ASHS.
		  -x file           Outer cross-validation loop specification file (see below)
		  -N                No overriding of ANTS/FLIRT results. If a result from an earlier run
		                    exists, don't run ANTS/FLIRT again
		  -Q                Use Sun Grid Engine (SGE) to schedule sub-tasks in each stage. By default,
		                    the whole ashs_train job runs in a single process. If you are on a cluster 
		                    that has SGE, you should really use this flag
		  -q OPTS           Pass in additional options to SGE's qsub. Also enables -Q option above.
		  -z script         Provide a path to an executable script that will be used to retrieve SGE or
		                    GNU parallel options for different stages of ASHS. Takes precendence over -q
		  -P                Use GNU parallel to run on multiple cores on the local machine. You need to
		                    have GNU parallel installed.
		  -C file           Configuration file. If not passed, uses $ASHS_ROOT/bin/ashs_config.sh
		  -V                Display version information and exit

		stages:
		  1                 Initialize atlas directory
		  2                 Build population-specific template
		  3                 Resample each atlas to template space
		  4                 Perform n^2 registration between all atlases
		  5                 Train AdaBoost method
		  6                 Organize final directory
		  7                 Cross-validation      

		data file:
		  The datafile is a comma-separated file listing, with rows corresponding to input images
		  For each image, the following columns are included:
		  
		    ID T1_MRI TSE_MRI SEG_LEFT SEG_RIGHT

		label description file:
		  This file is used to define the labels in the segmentation protocol. It is possible to 
		  have the same strugure (e.g., CA1) to have the same label on the left and on the right,
		  or to have then have separate labels. ASHS will handle it either way. The label file is
		  in the ITK-SNAP label format, i.e., a format of the file obtained by running ITK-SNAP
		  and selecting from the menu "Segmentation->Save Label Descriptions". Each line in the file
		  consists of seven numbers separated by whitespace, followed by a string in quotation marks.
		
		      0     0    0    0        0  0  0    "Clear Label"
		      1   255    0    0        1  1  1    "CA1"
		      2     0  255    0        1  1  1    "CA2"
		
		  The entries are: (1) the index of the label; (2-4) the R/G/B components of the color 
		  corresponding to the label (used to generate some figures in ASHS); (5-7) affect the
		  visibility of the label in SNAP; set them to 0 0 0 for the background and 1 1 1 for the
		  foreground labels. The name of the label is in quotation marks. The file can also include
		  comments (lines that start with a # character). 
		  
		data requirements:
		  The ASHS_TSE image slice direction should be z. In other words, the dimension
		  of ASHS_TSE image should be 400x400x30 or something like that, not 400x30x400
		
		cross-validation file (-x option):
		  This file specifies what kind of a cross-validation experiment to perform to test the
		  segmentation performance across the atlases. Each row corresponds to one experiment.
		  Each column lists the IDS of the test subjects for that experiment (i.e., the subjects
		  left out from atlas building). Cross-validation is optional.
		
		heuristic rules (-r option):
		  Rules can be specified to restrict certain labels in the segmentation to certain slices.
		  This is needed when the underlying segmentation protocol covers a portion of the slices
		  in the hippocampus, and the automatic method needs to be consistent with that protocol.
		  See the Yushkevich et al., 2010 Neuroimage paper for the heuristic rules applied there.
		  The format of the heuristic file can be found in the help for the subfield_slice_rules
		  program.

		SGE Options:
		  You can have detailed control over SGE options by passing a custom shell script to the -z
		  option. ASHS will call this shell script with the working directory as the first parameter
		  and stage as the second parameter. The script should print to stdout the SGE (qsub) options
		  that should be used for this stage. This allows you to allocate resources for each stage.
	USAGETEXT
}

# Dereference a link - different calls on different systems
function dereflink ()
{
  if [[ $(uname) == "Darwin" ]]; then
    local SLTARG=$(readlink $1)
    if [[ $SLTARG ]]; then
      echo $SLTARG
    else
      echo $1
    fi
  else
    readlink -f $1
  fi
}

# Print usage by default
if [[ $# -lt 1 ]]; then
  echo "Try $0 -h for more information."
  exit 2
fi

# Special actions (e.g., print version info, etc.)
unset ASHS_SPECIAL_ACTION

# Read the options
while getopts "C:D:L:w:s:x:q:r:z:NdhVQP" opt; do
  case $opt in

    D) ASHS_TRAIN_MANIFEST=$(dereflink $OPTARG);;
    L) ASHS_LABELFILE=$(dereflink $OPTARG);;
    w) ASHS_WORK=$OPTARG;;
    s) STAGE_SPEC=$OPTARG;;
    N) ASHS_SKIP_ANTS=1; ASHS_SKIP_RIGID=1; ASHS_SKIP=1;;
    Q) ASHS_USE_QSUB=1;;
    P) ASHS_USE_PARALLEL=1;;
    q) ASHS_USE_QSUB=1; ASHS_QSUB_OPTS=$OPTARG;;
    z) ASHS_USE_QSUB=1; ASHS_QSUB_HOOK=$OPTARG;;
    C) ASHS_CONFIG=$(dereflink $OPTARG);;
    r) ASHS_HEURISTICS=$(dereflink $OPTARG);;
    x) ASHS_XVAL=$(dereflink $OPTARG);;
    d) set -x -e;;
    h) usage; exit 0;;
    V) ASHS_SPECIAL_ACTION=vers;;
    \?) echo "Unknown option $OPTARG"; exit 2;;
    :) echo "Option $OPTARG requires an argument"; exit 2;;

  esac
done

# Check the root dir
if [[ ! $ASHS_ROOT ]]; then
  echo "Please set ASHS_ROOT to the ASHS root directory before running $0"
  exit -2
elif [[ $ASHS_ROOT != $(dereflink $ASHS_ROOT) ]]; then
  echo "ASHS_ROOT must point to an absolute path, not a relative path"
  exit -2
fi


# Create the working directory and the dump directory
mkdir -p $ASHS_WORK/dump

# Get rid of symlinks in the work path and make it global
ASHS_WORK=$(dereflink $ASHS_WORK)

# Redirect output/error to a log file in the dump directory
LOCAL_LOG=$(date +ashs_train.o%Y%m%d_%H%M%S)
mkdir -p $ASHS_WORK/dump
exec > >(tee -i $ASHS_WORK/dump/$LOCAL_LOG)
exec 2>&1

# Set the config file
if [[ ! $ASHS_CONFIG ]]; then
  ASHS_CONFIG=$ASHS_ROOT/bin/ashs_config.sh
fi

# Load the library and read the config file in the process
source $ASHS_ROOT/bin/ashs_lib.sh

# Just print version?
if [[ $ASHS_SPECIAL_ACTION == "vers" ]]; then
  vers
  exit 0
fi

# Check the listfile
if [[ ! -f $ASHS_TRAIN_MANIFEST ]]; then
  echo "Missing data list file (-D)"
  exit 1;
fi

if [[ ! -f $ASHS_LABELFILE ]]; then
  echo "Missing label description file (-L)"
  exit -1;
fi

# Check that parallel and qsub are not both on
if [[ $ASHS_USE_PARALLEL && $ASHS_USE_QSUB ]]; then
  echo "Cannot use SGE (-Q) and GNU Parallel (-P) at the same time"
  exit -2
fi


# Get the list of ids
ATLAS_ID=( $(cat $ASHS_TRAIN_MANIFEST | awk '! /^[ \t]*#/ {print $1}') );
ATLAS_T1=( $(cat $ASHS_TRAIN_MANIFEST | awk '! /^[ \t]*#/ {print $2}') );
ATLAS_T2=( $(cat $ASHS_TRAIN_MANIFEST | awk '! /^[ \t]*#/ {print $3}') );
ATLAS_LS=( $(cat $ASHS_TRAIN_MANIFEST | awk '! /^[ \t]*#/ {print $4}') );
ATLAS_RS=( $(cat $ASHS_TRAIN_MANIFEST | awk '! /^[ \t]*#/ {print $5}') );

# Sides - from config file
SIDES="$ASHS_SIDES"

# Get the number of atlases
N=${#ATLAS_ID[*]}

# Check that all the input files exist
for ((i=0; i < $N; i++)); do
  if [[ -f ${ATLAS_T1[i]} && -f ${ATLAS_T2[i]} && -f ${ATLAS_LS[$i]} && -f ${ATLAS_RS[$i]} ]]; then
    echo Verified atlas ${ATLAS_ID[i]}
  else
    echo Bad specification for atlas \"${ATLAS_ID[i]}\"
    exit 1
  fi
done

# Check the heuristic file
if [[ $ASHS_HEURISTICS ]]; then
  subfield_slice_rules --check-rules $ASHS_HEURISTICS
  if [[ $? -ne 0 ]]; then
    echo "Bad heuristics rule file"
    exit 1
  fi
fi
 
# Create the final directory (why?)
mkdir -p $ASHS_WORK/final

# Whether we are using QSUB
if [[ $ASHS_USE_QSUB ]]; then
  if [[ ! $SGE_ROOT ]]; then
    echo "-Q flag used, but SGE is not present."
    exit -1;
  fi
  if [[ $ASHS_QSUB_HOOK ]]; then
    if [[ ! -f $ASHS_QSUB_HOOK ]]; then
      echo "Parameter to -z ($ASHS_QSUB_HOOK) does not point to a file"
      exit 2
    fi
    echo "Using SGE with root $SGE_ROOT and callback script $ASHS_QSUB_HOOK"
  elif [[ $ASHS_QSUB_OPTS ]]; then 
    echo "Using SGE with root $SGE_ROOT and options \"$ASHS_QSUB_OPTS\""
  else
    echo "Using SGE with root $SGE_ROOT and default options"
  fi
elif [[ $ASHS_USE_PARALLEL ]]; then
  echo "Using GNU parallel"
else
  echo "Not using SGE or GNU parallel"
fi

# Run the stages of the script
export ASHS_ROOT ASHS_BIN ASHS_WORK ASHS_SKIP_ANTS ASHS_SKIP_RIGID ASHS_BIN_ANTS ASHS_SKIP
export ASHS_BIN_FSL ASHS_CONFIG ASHS_HEURISTICS ASHS_XVAL ASHS_LABELFILE ASHS_USE_QSUB QOPTS
export ASHS_USE_PARALLEL ASHS_TRAIN_MANIFEST
export SIDES

# Set the start and end stages
if [[ $STAGE_SPEC ]]; then
  STAGE_START=$(echo $STAGE_SPEC | awk -F '-' '$0 ~ /^[0-9]+-*[0-9]*$/ {print $1}')
  STAGE_END=$(echo $STAGE_SPEC | awk -F '-' '$0 ~ /^[0-9]+-*[0-9]*$/ {print $NF}')
else
  STAGE_START=1
  STAGE_END=7
fi

if [[ ! $STAGE_END || ! $STAGE_START || $STAGE_START -le 0 || $STAGE_END -gt 7 ]]; then
  echo "Wrong stage specification -s $STAGE_SPEC"
  exit -1;
fi

declare -F

# Names of the different stages
STAGE_NAMES=(\
  "Initialize work directory" \
  "Build template" \
  "Resample atlases to template" \
  "Pairwise registration between atlases" \
  "AdaBoost training and cross-validation" \
  "Organize final atlas" \
  "Perform cross-validation")

# If starting at stage other than 1, check for the correct output of
# the previous stages
if [[ $STAGE_START -gt 1 ]]; then

  # Run the validity check
  ashs_check_train $((STAGE_START-1)) || exit -1

fi

# Go through the stages
for ((STAGE=$STAGE_START; STAGE<=$STAGE_END; STAGE++)); do

  # The desription of the current stage
  STAGE_TEXT=${STAGE_NAMES[STAGE-1]}
  echo "****************************************"
  echo "Starting stage $STAGE: $STAGE_TEXT"
  echo "****************************************"

  # Put together qsub options for this stage
  if [[ $ASHS_USE_QSUB ]]; then

    # Is there a callback script
    if [[ $ASHS_QSUB_HOOK ]]; then
      QOPTS="$(bash $ASHS_QSUB_HOOK $ASHS_WORK $STAGE)"
      echo "Qsub options for this stage: $QOPTS"
    else
      QOPTS="${ASHS_QSUB_OPTS}"
    fi
  fi

  case $STAGE in 

    1)
    # Initialize Directory
    ashs_atlas_initialize_directory;;

    2)
    # The first step is to build a template from the atlas images using the standard
    # code in ANTS. For this, we got to copy all the atlases to a common directory
    ashs_atlas_build_template;;

    3)
    # Resample atlas to template
    ashs_atlas_resample_tse_to_template;;

    4)
    # Perform pairwise registration between all atlases
    ashs_atlas_register_to_rest;;

    5)
    # Train error correction
    ashs_atlas_adaboost_train;;

    6)
    # Organize everything into an atlas that can be used with the main ASHS script
    ashs_atlas_organize_final;;

    7)
    # Final cross-validation
    ashs_atlas_organize_xval;;

  esac

  # Run the validity check
  ashs_check_train $STAGE || exit -1

done
