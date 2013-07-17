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

set -x -e

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_voting_qsub.sh
	Root: ${ASHS_ROOT?}
	Working directory: ${ASHS_WORK?}
	PATH: ${PATH?}
BLOCK1

# directory for the subfields (separate for different parameter values)
WSUB=$ASHS_WORK/bootstrap/fusion

# Directory for QA output
WQA=$ASHS_WORK/qa
mkdir -p $WQA

# SNAP label file
ASHS_LABELFILE=$ASHS_ATLAS/snap/snaplabels.txt

# Generate the 'final' segmentations
for side in left right; do

  for segtype in heur corr_usegray corr_nogray; do

    c3d $ASHS_WORK/tse.nii.gz $WSUB/lfseg_${segtype}_${side}.nii.gz \
      -int 0 -reslice-identity -type ushort \
      -o $ASHS_WORK/final/${ASHS_SUBJID}_${side}_lfseg_${segtype}.nii.gz

  done
done

# Names of segmentations
for side in left right; do

	SMASV=$WSUB/lfseg_heur_${side}.nii.gz
	SBC=$WSUB/lfseg_corr_usegray_${side}.nii.gz

	if [[ -f $SMASV && -f $SBC ]]; then

		# Another lovely C3D program...
		c3d -verbose $SMASV -int 0 -trim 10x10x0mm -resample 400x400x100% -as REF \
			$ASHS_WORK/tse.nii.gz -stretch 1% 99% 0 255 -clip 0 255 -reslice-identity \
			-push REF $SBC -reslice-identity \
			-push REF -foreach -slice z 50% -flip xy -info -endfor \
			-popas S -popas S2 -as G -type uchar -info -o $WQA/gray_${side}.png \
			-clear -push G -push S -oli $ASHS_LABELFILE 0.5 -omc 3 $WQA/masv_${side}.png \
			-clear -push G -push S2 -oli $ASHS_LABELFILE 0.5 -omc 3 $WQA/final_${side}.png \
			-clear -push G -push S -push S2 -scale -1 -add -thresh 0 0 0 1 \
			-oli $ASHS_LABELFILE 0.5 -omc 3 $WQA/bc_${side}.png

	fi

done

# Create some latex
cd $WQA
cat > summary.tex <<LATEX
\\section{Subject: $ASHS_WORK}

\\setlength{\\tabcolsep}{0.5mm}
\\centering
\\begin{tabular}{cccc}
\\multicolumn{4}{c}{Left Hippocampus}\\\\
T2-MRI & InitSeg & BiasDetect & FinalSeg\\\\
\\includegraphics[width=0.23\\textwidth]{$WQA/gray_left}&
\\includegraphics[width=0.23\\textwidth]{$WQA/masv_left}&
\\includegraphics[width=0.23\\textwidth]{$WQA/bc_left}&
\\includegraphics[width=0.23\\textwidth]{$WQA/final_left}\\medskip\\\\
\\multicolumn{4}{c}{Right Hippocampus}\\\\
T2-MRI & InitSeg & BiasDetect & FinalSeg\\\\
\\includegraphics[width=0.23\\textwidth]{$WQA/gray_right}&
\\includegraphics[width=0.23\\textwidth]{$WQA/masv_right}&
\\includegraphics[width=0.23\\textwidth]{$WQA/bc_right}&
\\includegraphics[width=0.23\\textwidth]{$WQA/final_right}\\\\
\\end{tabular}
LATEX

