# CMake generated Testfile for 
# Source directory: /mnt/build/pauly/ants/advants/Examples
# Build directory: /mnt/build/pauly/ants/gcc64rel
# 
# This file includes the relevent testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(ANTS_CC_1 "ANTS" "2" "-m" "PR[" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii," "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii," "1" ",2" "]" "-r" "Gauss[" "3" "," "0" "]" "-t" "SyN[" "0.5" "]" "-i" "50x50x30" "-o" "/mnt/build/pauly/ants/gcc64rel/TEST.nii.gz")
ADD_TEST(ANTS_CC_1_WARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTWarp.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_CC_1_JPG "ConvertToJpg" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "ANTSCC1.jpg")
ADD_TEST(ANTS_CC_1_WARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "11.9992" "0.05")
ADD_TEST(ANTS_CC_1_WARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.61" "0.05")
ADD_TEST(ANTS_CC_1_WARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.00038593" "0.05")
ADD_TEST(ANTS_CC_1_INVERSEWARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "-i" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "/mnt/build/pauly/ants/gcc64rel/TESTInverseWarp.nii.gz" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_CC_1_INVERSEWARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "12.1606" "0.05")
ADD_TEST(ANTS_CC_1_INVERSEWARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.61" "0.05")
ADD_TEST(ANTS_CC_1_INVERSEWARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.000380545" "0.05")
ADD_TEST(ANTS_CC_2 "ANTS" "2" "-m" "PR[/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii,/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii,1,4]" "-r" "Gauss[3,0]" "-t" "SyN[0.5]" "-i" "50x50x30" "-o" "/mnt/build/pauly/ants/gcc64rel/TEST.nii.gz" "--go-faster" "true")
ADD_TEST(ANTS_CC_2_WARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTWarp.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_CC_2_JPG "ConvertToJpg" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "ANTSCC2.jpg")
ADD_TEST(ANTS_CC_2_WARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "11.7083" "0.05")
ADD_TEST(ANTS_CC_2_WARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.62" "0.05")
ADD_TEST(ANTS_CC_2_WARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.000461792" "0.05")
ADD_TEST(ANTS_CC_2_INVERSEWARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "-i" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "/mnt/build/pauly/ants/gcc64rel/TESTInverseWarp.nii.gz" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_CC_2_INVERSEWARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "11.912" "0.05")
ADD_TEST(ANTS_CC_2_INVERSEWARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.62" "0.05")
ADD_TEST(ANTS_CC_2_INVERSEWARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.000442839" "0.05")
ADD_TEST(ANTS_CC_3 "ANTS" "2" "-m" "CC[/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii,/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii,1,4]" "-r" "Gauss[3,0]" "-t" "SyN[0.5]" "-i" "50x50x30" "-o" "/mnt/build/pauly/ants/gcc64rel/TEST.nii.gz" "--go-faster" "true")
ADD_TEST(ANTS_CC_3_WARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTWarp.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_CC_3_JPG "ConvertToJpg" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "ANTSCC2.jpg")
ADD_TEST(ANTS_CC_3_WARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "11.7083" "0.05")
ADD_TEST(ANTS_CC_3_WARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.62" "0.05")
ADD_TEST(ANTS_CC_3_WARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.000461792" "0.05")
ADD_TEST(ANTS_CC_3_INVERSEWARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "-i" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "/mnt/build/pauly/ants/gcc64rel/TESTInverseWarp.nii.gz" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_CC_3_INVERSEWARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "11.912" "0.05")
ADD_TEST(ANTS_CC_3_INVERSEWARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.62" "0.05")
ADD_TEST(ANTS_CC_3_INVERSEWARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.000442839" "0.05")
ADD_TEST(ANTS_MSQ "ANTS" "2" "-m" "MSQ[/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii,/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii,1,0]" "-r" "Gauss[3,0]" "-t" "SyN[0.25]" "-i" "50x50x30" "-o" "/mnt/build/pauly/ants/gcc64rel/TEST.nii.gz")
ADD_TEST(ANTS_MSQ_WARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTWarp.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_MSQ_JPG "ConvertToJpg" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "ANTSMSQ.jpg")
ADD_TEST(ANTS_MSQ_WARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "11.7416" "0.05")
ADD_TEST(ANTS_MSQ_WARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.6" "0.05")
ADD_TEST(ANTS_MSQ_WARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-1.2" "0.05")
ADD_TEST(ANTS_MSQ_INVERSEWARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "-i" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "/mnt/build/pauly/ants/gcc64rel/TESTInverseWarp.nii.gz" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_MSQ_INVERSEWARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "11.7845" "0.05")
ADD_TEST(ANTS_MSQ_INVERSEWARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.6" "0.05")
ADD_TEST(ANTS_MSQ_INVERSEWARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-1.2" "0.05")
ADD_TEST(ANTS_MI_1 "ANTS" "2" "-m" "MI[/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii,/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii,1,32]" "-r" "Gauss[3,0]" "-t" "SyN[0.25]" "-i" "50x50x30" "-o" "/mnt/build/pauly/ants/gcc64rel/TEST.nii.gz")
ADD_TEST(ANTS_MI_1_WARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTWarp.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_MI_1_JPG "ConvertToJpg" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "ANTSMI1.jpg")
ADD_TEST(ANTS_MI_1_WARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "12.4" "0.05")
ADD_TEST(ANTS_MI_1_WARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.45" "0.05")
ADD_TEST(ANTS_MI_1_WARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.000370686" "0.05")
ADD_TEST(ANTS_MI_1_INVERSEWARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "-i" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "/mnt/build/pauly/ants/gcc64rel/TESTInverseWarp.nii.gz" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_MI_1_INVERSEWARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "12.6" "0.05")
ADD_TEST(ANTS_MI_1_INVERSEWARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.45" "0.05")
ADD_TEST(ANTS_MI_1_INVERSEWARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.000366879" "0.05")
ADD_TEST(ANTS_MI_2 "ANTS" "2" "-m" "MI[/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii,/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii,1,24]" "-r" "Gauss[3,0]" "-t" "SyN[0.25]" "-i" "50x50x20" "-o" "/mnt/build/pauly/ants/gcc64rel/TEST.nii.gz")
ADD_TEST(ANTS_MI_2_WARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTWarp.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_MI_2_JPG "ConvertToJpg" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "ANTSMI2.jpg")
ADD_TEST(ANTS_MI_2_WARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "12.4" "0.05")
ADD_TEST(ANTS_MI_2_WARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.45" "0.05")
ADD_TEST(ANTS_MI_2_WARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.000370686" "0.05")
ADD_TEST(ANTS_MI_2_INVERSEWARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "-i" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "/mnt/build/pauly/ants/gcc64rel/TESTInverseWarp.nii.gz" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_MI_2_INVERSEWARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "12.6" "0.05")
ADD_TEST(ANTS_MI_2_INVERSEWARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.45" "0.05")
ADD_TEST(ANTS_MI_2_INVERSEWARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.000366879" "0.05")
ADD_TEST(ANTS_ELASTIC "ANTS" "2" "-m" "PR[/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii,/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii,1,2]" "-t" "Elast[1]" "-i" "50x50x50" "-r" "Gauss[0,1]" "-o" "/mnt/build/pauly/ants/gcc64rel/TEST.nii.gz")
ADD_TEST(ANTS_ELASTIC_WARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTWarp.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_ELASTIC_JPG "ConvertToJpg" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "ANTSELASTIC.jpg")
ADD_TEST(ANTS_ELASTIC_WARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "12.5" "0.05")
ADD_TEST(ANTS_ELASTIC_WARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.4" "0.05")
ADD_TEST(ANTS_ELASTIC_WARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.000459869" "0.05")
ADD_TEST(ANTS_GSYN "ANTS" "2" "-m" "CC[/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii,/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii,1,3]" "-t" "SyN[0.25]" "-i" "50x50x50" "-r" "Gauss[3,0.]" "-o" "/mnt/build/pauly/ants/gcc64rel/TEST.nii.gz")
ADD_TEST(ANTS_GSYN_WARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTWarp.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_GSYN_JPG "ConvertToJpg" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "ANTSGSYN.jpg")
ADD_TEST(ANTS_GSYN_WARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "11.7734" "0.05")
ADD_TEST(ANTS_GSYN_WARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.6" "0.05")
ADD_TEST(ANTS_GSYN_WARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.000478672" "0.05")
ADD_TEST(ANTS_GSYN_INVERSEWARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "-i" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "/mnt/build/pauly/ants/gcc64rel/TESTInverseWarp.nii.gz" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_GSYN_JPGINV "ConvertToJpg" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "ANTSGSYNINV.jpg")
ADD_TEST(ANTS_GSYN_INVERSEWARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "12.0541" "0.05")
ADD_TEST(ANTS_GSYN_INVERSEWARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.6" "0.05")
ADD_TEST(ANTS_GSYN_INVERSEWARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.000475175" "0.05")
ADD_TEST(ANTS_EXP "ANTS" "2" "-m" "PR[/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii,/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii,1,4]" "-t" "Exp[0.5,2,0.5]" "-i" "50x50x50" "-r" "Gauss[0.5,0.25]" "-o" "/mnt/build/pauly/ants/gcc64rel/TEST.nii.gz")
ADD_TEST(ANTS_EXP_WARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTWarp.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_EXP_JPG "ConvertToJpg" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "ANTSEXP.jpg")
ADD_TEST(ANTS_EXP_WARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "12" "0.05")
ADD_TEST(ANTS_EXP_WARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.5" "0.05")
ADD_TEST(ANTS_EXP_WARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.000423693" "0.05")
ADD_TEST(ANTS_EXP_INVERSEWARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "-i" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "/mnt/build/pauly/ants/gcc64rel/TESTInverseWarp.nii.gz" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_EXP_INVERSEWARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "12.349" "0.05")
ADD_TEST(ANTS_EXP_INVERSEWARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.5" "0.05")
ADD_TEST(ANTS_EXP_INVERSEWARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.000339606" "0.05")
ADD_TEST(ANTS_SYN "ANTS" "2" "-m" "PR[/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii,/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii,1,2]" "-t" "SyN[0.5,2,0.05]" "-i" "50x50x50" "-r" "Gauss[3,0.0,32]" "-o" "/mnt/build/pauly/ants/gcc64rel/TEST.nii.gz")
ADD_TEST(ANTS_SYN_WARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTWarp.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_SYN_JPG "ConvertToJpg" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "ANTSSYN.jpg")
ADD_TEST(ANTS_SYN_WARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "12.0239" "0.05")
ADD_TEST(ANTS_SYN_WARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.6" "0.05")
ADD_TEST(ANTS_SYN_WARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/warped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.000461922" "0.05")
ADD_TEST(ANTS_SYN_INVERSEWARP "WarpImageMultiTransform" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "-i" "/mnt/build/pauly/ants/gcc64rel/TESTAffine.txt" "/mnt/build/pauly/ants/gcc64rel/TESTInverseWarp.nii.gz" "-R" "/mnt/build/pauly/ants/advants/Examples/Data/r16slice.nii")
ADD_TEST(ANTS_SYN_JPGINV "ConvertToJpg" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "ANTSSYNINV.jpg")
ADD_TEST(ANTS_SYN_INVERSEWARP_METRIC_0 "MeasureImageSimilarity" "2" "0" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "12.5104" "0.05")
ADD_TEST(ANTS_SYN_INVERSEWARP_METRIC_1 "MeasureImageSimilarity" "2" "1" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.6" "0.05")
ADD_TEST(ANTS_SYN_INVERSEWARP_METRIC_2 "MeasureImageSimilarity" "2" "2" "/mnt/build/pauly/ants/advants/Examples/Data/r64slice.nii" "/mnt/build/pauly/ants/gcc64rel/inversewarped.nii.gz" "/mnt/build/pauly/ants/gcc64rel/TESTlog.txt" "/mnt/build/pauly/ants/gcc64rel/TESTmetric.nii.gz" "-0.000444279" "0.05")
ADD_TEST(ANTS_PSE_MSQ_TXT "ANTS" "2" "-i" "91x70x55x40x30" "-r" "Gauss[3,0.,32]" "-t" "SyN[0.25]" "-m" "MSQ[/mnt/build/pauly/ants/advants/Examples/Data/Frown.nii,/mnt/build/pauly/ants/advants/Examples/Data/Smile.nii,1,0]" "-m" "PSE[/mnt/build/pauly/ants/advants/Examples/Data/Frown.nii,/mnt/build/pauly/ants/advants/Examples/Data/Smile.nii,/mnt/build/pauly/ants/advants/Examples/Data/Frown.txt,/mnt/build/pauly/ants/advants/Examples/Data/Smile.txt,1,0.33,11,1,25]" "--continue-affine" "0" "--number-of-affine-iterations" "0" "-o" "/mnt/build/pauly/ants/gcc64rel/TEST.nii.gz")
ADD_TEST(ANTS_PSE_MSQ_VTK "ANTS" "2" "-i" "91x70x55x40x30" "-r" "Gauss[3,0.,32]" "-t" "SyN[0.25]" "-m" "MSQ[/mnt/build/pauly/ants/advants/Examples/Data/Frown.nii,/mnt/build/pauly/ants/advants/Examples/Data/Smile.nii,1,0]" "-m" "PSE[/mnt/build/pauly/ants/advants/Examples/Data/Frown.nii,/mnt/build/pauly/ants/advants/Examples/Data/Smile.nii,/mnt/build/pauly/ants/advants/Examples/Data/Frown.vtk,/mnt/build/pauly/ants/advants/Examples/Data/Smile.vtk,1,0.33,11,1,25]" "--continue-affine" "0" "--number-of-affine-iterations" "0" "-o" "/mnt/build/pauly/ants/gcc64rel/TEST.nii.gz")
