# Modeling application to *in vivo* cortical activity
Main analysis code for estimating the contribution V4, LIP, 7A and FEF to the N2pc from laminar field recordings.

## Description 
CSD data used to forward model EEG can be found through Data Dryad (https://doi.org/10.5061/dryad.9ghx3ffm4). CSD session averages are provided in sw_hemi_icsd.mat file.

Lead field matrices used for the EEG forward models of each area (V4, LIP, 7A, and FEF) are provided in the [lead_fields](lead_fields) folder. 

All EEG topographic maps and 2D layouts were created with [Brainstorm](https://neuroimage.usc.edu/brainstorm/Introduction).

### Files
- [cal_spline_iCSD.m](calculate_spline_iCSD.m): code used to calculate the current source density (CSD) maps from the laminar field potential recordings employing the spline iCSD method. The script calls the function [iCSDspline_Pettersen()](iCSDspline_Pettersen.m) which uses the functions from the CSDplotter toolbox.
- [cal_dip_CSD()](cal_dip_CSD.m): function to calculate the current dipole moment from the CSD.
- [cal_dip_EEG_Fig3.m](cal_dip_EEG_Fig3.m): calculates the grand average CSD for each condition (contra and ipsi), its associated current dipole moment, and the resulting EEG distributions. Results are reported in Figure 3.
- [monkey_Y_head_model_surfaces](monkey_Y_head_model_surfaces): surfaces of the cortex, scalp, and skull obtained from the segmentation of the T1-weighted MRI of monkey Y, used to construct the head model in Figure 3C. We utilized this head model on all EEG estimations from the in-vivo laminar recordings in macaques.

## Dependencies
- MATLAB (release 2021b, [The MathWorks](https://www.mathworks.com/?s_tid=gn_logo)). Toolboxes: Signal Processing Toolbox, Statistics and Machine Learning Toolbox.
- [CSDplotter](matlab_ana_scripts/functions/CSDplotter-0.1.1) (https://github.com/espenhgn/CSDplotter)
- [Brainstorm](https://neuroimage.usc.edu/brainstorm/Introduction)
