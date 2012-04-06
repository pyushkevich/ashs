#!/bin/bash - 

#######################################################################
#
#  Program:   ASHS (Automatic Segmentation of Hippocampal Subfields)
#  Module:    $Id: $
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
      -q string         List of additional options to pass to qsub (Sun Grid Engine)
      -C file           Configuration file. If not passed, uses $ASHS_ROOT/bin/ashs_config.sh
      -V                Display version information and exit
		stages:

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
      The TSE image slice direction should be z. In other words, the dimension
      of TSE image should be 400x400x30 or something like that, not 400x30x400

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
	USAGETEXT
}

# Default config
ASHS_CONFIG=$ASHS_ROOT/bin/ashs_config.sh

# Load the library
source $ASHS_ROOT/bin/ashs_lib.sh

# Common code
source ashs_common_master.sh

# Print usage by default
if [[ $# -lt 1 ]]; then
  echo "Try $0 -h for more information."
  exit 2
fi

# Read the options
while getopts "C:D:L:w:s:x:q:r:NdhV" opt; do
  case $opt in

    D) LISTFILE=$OPTARG;;
    L) LABELFILE=$OPTARG;;
    w) WORK=$OPTARG;;
		s) STAGE_SPEC=$OPTARG;;
		N) SKIP_ANTS=1; SKIP_RIGID=1; ;;
    q) QOPTS=$OPTARG;;
    C) ASHS_CONFIG=$OPTARG;;
    r) ASHS_HEURISTICS=$(readlink -f $OPTARG);;
    x) XVAL=$OPTARG;;
    d) set -x -e;;
    h) usage; exit 0;;
    V) vers; exit 0;;
    \?) echo "Unknown option $OPTARG"; exit 2;;
    :) echo "Option $OPTARG requires an argument"; exit 2;;

  esac
done

# Check the listfile
if [[ ! -f $LISTFILE ]]; then
  echo "Missing data list file (-D)"
  exit 1;
fi

if [[ ! -f $LABELFILE ]]; then
  echo "Missing label description file (-L)"
  exit -1;
fi

# Read the config file
source $ASHS_CONFIG

# Get the list of ids
ATLAS_ID=( $(cat $LISTFILE | awk '! /^[ \t]*#/ {print $1}') );
ATLAS_T1=( $(cat $LISTFILE | awk '! /^[ \t]*#/ {print $2}') );
ATLAS_T2=( $(cat $LISTFILE | awk '! /^[ \t]*#/ {print $3}') );
ATLAS_LS=( $(cat $LISTFILE | awk '! /^[ \t]*#/ {print $4}') );
ATLAS_RS=( $(cat $LISTFILE | awk '! /^[ \t]*#/ {print $5}') );

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
 
# Create the working directory and the dump directory
mkdir -p $WORK $WORK/dump $WORK/final

# Run the stages of the script
ROOT=$ASHS_ROOT;
export ROOT PATH ASHS_BIN WORK SKIP_ANTS SKIP_RIGID ASHS_BIN_ANTS 
export ASHS_BIN_FSL ASHS_CONFIG ASHS_HEURISTICS XVAL LABELFILE

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

declare -F

# Go through the stages
for ((STAGE=$STAGE_START; STAGE<=$STAGE_END; STAGE++)); do

  case $STAGE in 

    0)

    # Initialize Directory
    echo "Running stage 0: initialize work directory"
    ashs_atlas_initialize_directory;;

    1)

    # The first step is to build a template from the atlas images using the standard
    # code in ANTS. For this, we got to copy all the atlases to a common directory
    echo "Running stage 1: build template"
    ashs_atlas_build_template;;

    2)

    echo "Running stage 2: resample all T2 data to template"
    ashs_atlas_resample_tse_to_template;;

    3)

    # Perform pairwise registration between all atlases
    echo "Running stage 3: pairwise registration between atlases"
    ashs_atlas_register_to_rest;;

    4)

    # Perform cross-validation experiments
    echo "Running stage 4: AdaBoost training and cross-validation"
    ashs_atlas_adaboost_train;;

    5)

    # Organize everything into an atlas that can be used with the main ASHS script
    echo "Running stage 5: Organize the output directory"
    ashs_atlas_organize_final;;


  esac

done
