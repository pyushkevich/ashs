#!/bin/bash

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


# Common code executed at startup of all 'master' scripts, i.e., ashs_main,
# ashs_train

# The ROOT variable must be set
if [[ ! $ASHS_ROOT ]]; then
  echo "ASHS_ROOT is not set. Please set this variable to point to the root ASHS directory"
  exit -1;
fi

# Get the architecture and check ability to run binaries
ARCH=$(uname);
ASHS_BIN=$ASHS_ROOT/ext/$ARCH/bin
#ASHS_ANTS=$ASHS_BIN/ants
ASHS_ANTS=$ASHS_BIN/ants_1042
ASHS_FSL=$ASHS_BIN/fsl
if [[ ! $($ASHS_BIN/c3d -version | grep 'Version') ]]; then
  echo "Can not execute command \'$ASHS_BIN/c3d -version\'. Wrong architecture?"
  exit -1
fi

# Set the path for the ASHS programs, to ensure that we don't call some other
# version of ants or c3d. This is nicer than having to prefix every call to
# c3d by the path. 
PATH=$ASHS_BIN:$ASHS_ANTS:$ASHS_FSL:$PATH
