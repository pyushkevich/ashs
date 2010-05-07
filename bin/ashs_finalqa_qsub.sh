#!/bin/bash
#$ -S /bin/bash
set -x -e

# Verify all the necessary inputs
cat <<-BLOCK1
	Script: ashs_voting_qsub.sh
	Root: ${ROOT?}
	Working directory: ${WORK?}
	PATH: ${PATH?}
BLOCK1

# directory for the subfields (separate for different parameter values)
WSUB=$WORK/subfields

# Directory for QA output
WQA=$WORK/qa
mkdir -p $WQA

# Names of segmentations
for side in left right; do
	SMASV=$WSUB/consensus_heuristic_wgtavg_${side}_native.nii.gz
	SBC=$WSUB/bcfh_heuristic_wgtavg_${side}_native.nii.gz
	if [[ -f $SMASV && -f $SBC ]]; then

		# Another lovely C3D program...
		$BIN/c3d -verbose $SMASV -trim-to-size 32x32x40mm -int 0 -resample 400x400x100% -as REF \
			$WORK/tse.nii.gz -stretch 1% 99% 0 255 -clip 0 255 -reslice-identity \
			-push REF $SBC -reslice-identity \
			-push REF -foreach -slice z 50% -flip xy -info -endfor \
			-popas S -popas S2 -as G -type uchar -info -o $WQA/gray_${side}.png \
			-clear -push G -push S -oli $ROOT/data/snap/snaplabels.txt 0.5 -omc 3 $WQA/masv_${side}.png \
			-clear -push G -push S2 -oli $ROOT/data/snap/snaplabels.txt 0.5 -omc 3 $WQA/final_${side}.png \
			-clear -push G -push S -push S2 -scale -1 -add -thresh 0 0 0 1 \
			-oli $ROOT/data/snap/snaplabels.txt 0.5 -omc 3 $WQA/bc_${side}.png

	fi
done

# Create some latex
cd $WQA
cat > summary.tex <<LATEX
\\section{Subject: $WORK}

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

