# python_sim_scripts

Code to reproduce biophysical simulations in Figure 2 of Herrera and Westerberg et al. (2022). 

### Run simulation 
0. compile mod mechanisms. 
1. run simulations:
    1. To simulate a population of *N* L3 pyramidal cells described by [Eyal et al. (2018)](https://doi.org/10.3389/fncel.2018.00181) model, run [popL3PCs_CaspikeParadigm.py](popL3PCs_CaspikeParadigm.py).
    2. To simulate a population of *N* L5 pyramidal cells described by [Hay et al. (2011)](https://doi.org/10.1371/journal.pcbi.1002107) model, run [popL5PCs_CaspikeParadigm.py](popL5PCs_CaspikeParadigm.py).

### Generate neurons location
- run [drawL3PCsposition_V4column.py](drawL3PCsposition_V4column.py) to generate the location of L3 pyramidal neurons (Figure2C).
- run [drawL5PCsposition_V4column.py](drawL5PCsposition_V4column.py) to generate the location of L5 pyramidal neurons (Figure2D).