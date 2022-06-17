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
RED='\033[0;31m'
YELLOW='\033[0;33m'
CYAN='\033[0;36m'
NC='\033[0m'

# ----------------------------------------------------
# A dummy hook file that is used by default (without -H option)
# ----------------------------------------------------
if [[ $1 == "progress" ]]; then

  echo "===================================="
  echo "         PROGRESS: $2               "
  echo "===================================="


elif [[ $1 == "error" ]]; then
  echo -e "${RED}"
  echo -e "**************** !! ERROR !! *******************"
  echo $2 | fold -w 48 -s 
  echo -e "************************************************"
  echo -e "${NC}"
elif [[ $1 == "warning" ]]; then
  echo -e "${YELLOW}"
  echo -e "==================  WARNING  ==================="
  echo $2 | fold -w 48 -s 
  echo -e "================================================"
  echo -e "${NC}"
elif [[ $1 == "info" ]]; then
  echo -e "${CYAN}"
  echo -e "-------------------  INFO  ---------------------"
  echo $2 | fold -w 48 -s 
  echo -e "------------------------------------------------"
  echo -e "${NC}"
fi


