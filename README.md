# Biophysical forward modeling of EEG signals from laminar current source density

Main analysis code for figures in Herrera and Westerberg et al. (2022), *"Resolving the mesoscopic missing link: Biophysical modeling of EEG from cortical columns in primates"*.

## Description
- [sim_PCsPops](sim_PCsPops): Biophysical simulations performed for Figure 2. 
- [sim_V4data](sim_V4data): Main analysis code for estimating the contribution V4, LIP, 7A and FEF to the N2pc from laminar field recordings.

## Main Dependencies
- MATLAB (release 2021b, [The MathWorks](https://www.mathworks.com/?s_tid=gn_logo))
- NEURON simulator (release 8.0, http://www.neuron.yale.edu/neuron)
- LFPy 2.2.2 (https://github.com/LFPy/LFPy)
- [CSDplotter](matlab_ana_scripts/functions/CSDplotter-0.1.1) (https://github.com/espenhgn/CSDplotter)
- [Brainstorm](https://neuroimage.usc.edu/brainstorm/Introduction)

## License
Copyright (C) 2021- Beatriz Herrera, Jacob A. Westerberg, Michelle S. Schall, Alexander Maier, Geoffrey F. Woodman, Jeffrey D. Schall, and Jorge J. Riera.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)

## Citation
Herrera, B., Westerberg, J. A., Schall, M. S., Maier, A., Woodman, G. F., Schall, J. D., & Riera, J. J. (2022). Resolving the mesoscopic missing link: Biophysical modeling of EEG from cortical columns in primates. bioRxiv. **doi**: https://doi.org/10.1101/2022.03.16.484595.
