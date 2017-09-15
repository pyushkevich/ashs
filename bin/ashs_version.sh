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

# ----------------------------------------------------
# Define ASHS version (semantic and release date)
#   e.g. 1.0.0-beta
# ----------------------------------------------------
ASHS_VERSION_MAJOR=1
ASHS_VERSION_MINOR=0
ASHS_VERSION_PATCH=0
ASHS_VERSION_NOTE=""
ASHS_VERSION_FULL="${ASHS_VERSION_MAJOR}.${ASHS_VERSION_MINOR}.${ASHS_VERSION_PATCH}${ASHS_VERSION_NOTE}"
ASHS_VERSION_DATE=20170915

# ----------------------------------------------------
# Atlas compatibility: atlases created before this
# date are considered outdated
# ----------------------------------------------------
ASHS_OLDEST_COMPAT_DATE=20170810
