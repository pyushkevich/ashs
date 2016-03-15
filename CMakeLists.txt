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
CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(ASHS)

# Build C3D
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/submodules/c3d)
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/submodules/c3d)
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/submodules/greedy)

# Build LabelFusion
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/src/LabelFusion)
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/src/CorrectiveLearning)
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/src/PointSet)
