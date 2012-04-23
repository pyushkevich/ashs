#!/bin/bash

#######################################################################
#
#  Program:   ASHS (Automatic Segmentation of Hippocampal Subfields)
#  Module:    $Id: ashs_atlas_initdir_qsub.sh 81 2012-04-20 14:34:50Z yushkevich $
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
		ashs_util_makepdf: Generate PDF from a grayscale image and segmentation
		usage:
		  ashs_util_makepdf [options] gray.nii seg.nii labelfile.txt outfile.pdf
		options:
		  -h              Print this help information 
		notes:
		  You need a working LaTeX installation!
USAGETEXT
}

# Run parent-level code
ASHS_CONFIG=${ASHS_ROOT?}/bin/ashs_config.sh
source ${ASHS_ROOT?}/bin/ashs_lib.sh

# Read command line
while getopts "h" opt; do
  case $opt in

    h) usage; exit 0;;
    \?) echo "Unknown option $OPTARG"; exit 2;;
    :) echo "Option $OPTARG requires an argument"; exit 2;;

  esac
done

# Read required arguments
shift $((OPTIND-1))
TSE=${1?}
SEG=${2?}
LABELFILE=${3?}
OUTPUT=${4?}

# Create a temp directory
TMP=$(get_tmpdir)
echo $TMP

# Get the number of slices in the image
NSL=$(c3d $SEG -info-full | awk '$1 ~ /Image/ && $2 ~ /Dimensions/ { print $6 * 1.0 }')

# Export all the slices
for ((isl=0;isl<$NSL;isl++)); do

  idx=$(printf %04d $isl)

  # Another lovely C3D program...
  c3d -verbose $SEG -int 0 -trim 10x10x200mm -resample 400x400x100% -as REF \
    $TSE -stretch 1% 99% 0 255 -clip 0 255 -reslice-identity \
    -push REF -foreach -slice z $isl -flip xy -endfor \
    -popas S -as G -type uchar -o $TMP/gray${idx}.png \
    -clear -push G -push S -oli $LABELFILE 0.5 -omc 3 $TMP/segm${idx}.png > /dev/null

done

# Write LaTeX code
pushd $TMP

cat > doc.tex <<-LATEXSCRIPT
		\documentclass{article}
		\usepackage{graphicx}
		\usepackage[landscape,margin=0.5in]{geometry}
		\begin{document}
LATEXSCRIPT

for ((isl=0;isl<$NSL;isl++)); do

  idx=$(printf %04d $isl)

  cat >> doc.tex <<-SCRIPTLET
  		
		  \medskip
		  \begin{center}
		    \section*{Slice $((isl+1)) of ${NSL}}
		    \bigskip
		    \begin{tabular}{cc}
		      \includegraphics[width=0.48\textwidth]{gray${idx}} &
		      \includegraphics[width=0.48\textwidth]{segm${idx}}
		    \end{tabular}
		  \end{center}
		  \cleardoublepage
	SCRIPTLET

done

echo '\end{document}' >> doc.tex 

# run latex
pdflatex doc.tex

popd

if [[ -f $TMP/doc.pdf ]]; then

	cp -a $TMP/doc.pdf $OUTPUT
	# rm -rf $TMP

fi

 

