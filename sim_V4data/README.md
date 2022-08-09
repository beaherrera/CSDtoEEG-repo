# sim_V4data
Main analysis code for estimating the contribution V4, LIP, 7A and FEF to the N2pc from laminar field recordings.

## Description
- [lead_fields](lead_fields): lead fields used to estimate the contribution of V4, LIP, 7A and FEF to the N2pc. 
- [monkey_Y_head_model_surfaces](monkey_Y_head_model_surfaces): surfaces of the cortex, scalp, and skull obtained from the segmentation of the T1-weighted MRI of monkey Y, used to construct the head model in Figure 3C. We utilized this head model on all EEG estimations from the in-vivo laminar recordings in macaques.
- [cal_spline_iCSD.m](calculate_spline_iCSD.m): code used to calculate the current source density (CSD) maps from the laminar field potential recordings in V4 using the spline iCSD method. The script calls the function [iCSDspline_Pettersen()](iCSDspline_Pettersen.m) which uses the methods provided in the CSDplotter toolbox.
- [cal_dip_CSD()](cal_dip_CSD.m): function to calculate the current dipole moment from the CSD.

## Dependencies
- MATLAB (release 2021b, [The MathWorks](https://www.mathworks.com/?s_tid=gn_logo)). Toolboxes: Signal Processing Toolbox, Statistics and Machine Learning Toolbox.
- [CSDplotter](matlab_ana_scripts/functions/CSDplotter-0.1.1) (https://github.com/espenhgn/CSDplotter)
- [Brainstorm](https://neuroimage.usc.edu/brainstorm/Introduction)