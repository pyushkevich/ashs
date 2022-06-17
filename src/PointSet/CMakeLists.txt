#######################################################################
#
#  Program:   ASHS (Automatic Segmentation of Hippocampal Subfields)
#  Module:    $Id$
#  Language:  C++ program
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

CMAKE_MINIMUM_REQUIRED(VERSION 2.6.2)
IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)  

PROJECT(POINTSET)

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

INCLUDE_DIRECTORIES(${POINTSET_SOURCE_DIR})

ADD_EXECUTABLE(ml_affine MultiLabelAffine.cxx)

#SET(VTK_LIBRARIES
#  vtkIO vtkHybrid vtkRendering vtkImaging vtkGraphics vtkFiltering vtkCommon)


SET(COMMON_LIBS
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES}
)

TARGET_LINK_LIBRARIES(ml_affine ${COMMON_LIBS} ${VTK_LIBRARIES})
