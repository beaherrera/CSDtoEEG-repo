# matlab_ana_scripts
matlab scripts used to process simulated data and create Figure 2 C-K.

## Description
### Scripts
- [Matlab_ana_scripts.prj](Matlab_ana_scripts.prj): MATLAB project object - run to set up the MATLAB path. Please, notice that you would still need to download and install [Brainstorm](https://neuroimage.usc.edu/brainstorm/Introduction).
- [cal_columnDip.m](cal_columnDip.m): script used to calculate the equivalent current dipole moment at the center of mass of a cylinder of 2,200 L3 and/or 1,000 L5 simulated pyramidal cells from summed transmembrane currents of all neurons (Data used in Figure 2F,H).
- [cal_EEG_V4_L3.m](cal_EEG_V4_L3.m) and [cal_EEG_V4_L3.m](cal_EEG_V4_L5.m) were used to create Figure 2 E-F. 
- [cal_EEG_V4_L3andL5PCsPop.m](cal_EEG_V4_L3andL5PCsPop.m) was used to create Figure 2G-H and calculate the EEG signals for Figure 2I.
- [cal_EEG_randVert.m](cal_EEG_randVert.m): script used to calculate the EEG across 15 random cortical column locations in the macaque brain. Data used in Figure 2 J-K.
- [cal_error_EEG.m](cal_error_EEG.m) calculates the relative diference (RDM) and magnitude ratio (MAG) measures comparing ground-truth compartment-based EEG to that derived from summed transmembrane currents of all neurons and CSD as a function of the depth of the center of mass of cylinder of 2,200 L3 and 1,000 L5 simulated pyramidal cells relative to the scalp. Results shown in Figure 2J-K.

### Folders
- [functions](functions): functions and toolboxes used for the analyses. (DOES NOT INCLUDE BRAINSTORM)
- [L3_simData_10ms](L3_simData_10ms): simulated data of the population of 2,200 L3 simulated pyramidal cells.
- [L5_simData_10ms](L5_simData_10ms): simulated data of the population of 1,000 L3 simulated pyramidal cells.
- [L3andL5_EEG](L3andL5_EEG): simulated data of the combined activity of 2,200 L3 and 1,000 L5 simulated pyramidal cells.
- [leadfields_NMTv2_atlas](leadfields_NMTv2_atlas): V4 lead field matrix (1 dipolar source in each hemisphere, Figure 2B).
- [NMTv2_atlas_surfaces](NMTv2_atlas_surfaces): symmetric surfaces provided in the NIMH Macaque Template version 2.0 [(Jung et al., 2021)](https://doi.org/10.1016/j.neuroimage.2021.117997) used to construct the head model (Figure 2A).
