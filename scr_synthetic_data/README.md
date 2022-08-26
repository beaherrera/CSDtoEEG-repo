# Modeling application to synthetic field potential recordings
Codes for generating and analyzing the synthetic field potential
recordings used to validate our proposed forward modeling approach linking mesoscopic CSD derived from LFP to macroscopic ERPs (Figure 2 - Herrera and Westerberg et al. (2022)).

### Description
- [python_sim_scripts](python_sim_scripts): scripts for running the biophysical simulations to generate the synthetic field potentials shown in Figure 2 of Herrera and Westerberg et al. (2022).
- [matlab_ana_scripts](matlab_ana_scripts): matlab scripts used to process the synthetic field potentials and create Figure 2 C-K.
 
### Dependencies
- MATLAB (release 2021b, [The MathWorks](https://www.mathworks.com/?s_tid=gn_logo)). Toolboxes: Signal Processing Toolbox, Parallel Computing Toolbox, Statistics and Machine Learning Toolbox.
- NEURON simulator (release 8.0, http://www.neuron.yale.edu/neuron)
- Python 3.8
- LFPy 2.2.2 (https://github.com/LFPy/LFPy)
- [CSDplotter](matlab_ana_scripts/functions/CSDplotter-0.1.1) (https://github.com/espenhgn/CSDplotter)
- [Brainstorm](https://neuroimage.usc.edu/brainstorm/Introduction)
- Elephant (https://elephant.readthedocs.io/en/latest/)
- Neo (https://neo.readthedocs.io/en/latest/)
- [eegfilt.m](matlab_ana_scripts/functions/eegfilt.m) function from [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php)
- [RasterPlot](https://www.mathworks.com/matlabcentral/fileexchange/45671-flexible-and-fast-spike-raster-plotting)
- [point2trimesh.m](matlab_ana_scripts/functions/point2trimesh/point2trimesh.m) (ver. 1.0.0, https://www.mathworks.com/matlabcentral/fileexchange/52882-point2trimesh-distance-between-point-and-triangulated-surface)
