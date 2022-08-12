# sim_V4data
Main analysis code for estimating the contribution V4, LIP, 7A and FEF to the N2pc from laminar field recordings.

## Description
- [lead_fields](lead_fields): lead field matrices used for estimating the EEG evoked by each area (V4, LIP, 7A, and FEF). Only the contribution of two dipoles, one in each hemisphere, was considered. 
- [monkey_Y_head_model_surfaces](monkey_Y_head_model_surfaces): surfaces of the cortex, scalp, and skull obtained from the segmentation of the T1-weighted MRI of monkey Y, used to construct the head model in Figure 3C. We utilized this head model on all EEG estimations from the in-vivo laminar recordings in macaques.
- [cal_spline_iCSD.m](calculate_spline_iCSD.m): code used to calculate the current source density (CSD) maps from the laminar field potential recordings employing the spline iCSD method. The script calls the function [iCSDspline_Pettersen()](iCSDspline_Pettersen.m) which calls the functions from the CSDplotter toolbox.
- [cal_dip_CSD()](cal_dip_CSD.m): function to calculate the current dipole moment from the CSD.

## Dependencies
- MATLAB (release 2021b, [The MathWorks](https://www.mathworks.com/?s_tid=gn_logo)). Toolboxes: Signal Processing Toolbox, Statistics and Machine Learning Toolbox.
- [CSDplotter](matlab_ana_scripts/functions/CSDplotter-0.1.1) (https://github.com/espenhgn/CSDplotter)
- [Brainstorm](https://neuroimage.usc.edu/brainstorm/Introduction)
