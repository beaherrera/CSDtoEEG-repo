# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 15:54:03 2021

@author: Beatriz Herrera
"""
from __future__ import division
from os.path import join
import numpy as np
import myfunctions
import matplotlib.pyplot as plt
from neuron import h
from matplotlib.collections import PolyCollection
import LFPy

np.random.seed(0)


''' Functions '''


def make_cellEyalEtAl2018(morphology, dt=2**-6, tstop=1000, cell_model=None,
                          active=False, add_synapses=False):
    r"""
    Create LFPy cell object using Eyal et al. 2018 model.

    Copied and addapted from file EEG-master\src\make_data_figure2.py
    in the python package provided in
    Næss, S., Halnes, G., Hagen, E., Hagler Jr, D. J., Dale, A. M.,
    Einevoll, G. T., & Ness, T. V. (2021).
    Biophysically detailed forward modeling of the neural origin of
    EEG and MEG signals.
    NeuroImage, 225, 117467.

    Parameters
    ----------
    morphology : str
        name of the cell morphology file.
    dt : float, optional
        integration time step. The default is 2**-6.
    tstop : float, optional
        length of the simulation. The default is 1000.
    cell_model : str, optional
        The default is None.
    active : boolean, optional
        The default is False.
    add_synapses : boolean, optional
        The default is False.

    Returns
    -------
    cell : LFPy.Cell Obj

    """
    model_folder = join("cell_models", "EyalEtAl2018")
    morph_path = join(model_folder, "Morphs", morphology)

    if cell_model is None:

        # loading mechanisms
        mod_folder = myfunctions.posixpth(join(model_folder, 'mechanisms'))
        if not hasattr(h, 'NMDA'):
            h.nrn_load_dll(mod_folder+"/nrnmech.dll")

        cell_parameters = {
            'v_init': -70,
            'morphology': morph_path,
            # S/cm^2, mV
            'passive_parameters': {'g_pas': 1./30000, 'e_pas': -70},
            'Ra': 150,  # Ω cm
            'cm': 1,  # µF/cm^2
            'nsegs_method': "lambda_f",
            "lambda_f": 100,
            'dt': 2**-4,  # [ms] Should be a power of 2
            'tstart': -10,  # [ms] Simulation start time
            'tstop': 100,  # [ms] Simulation end time
            "pt3d": True,
            'passive': True
        }
    else:

        # loading mechanisms
        mod_folder = myfunctions.posixpth(join(model_folder,
                                               'ActiveMechanisms'))
        if not hasattr(h, 'NaTg'):  # mod_folder in neuron.nrn_dll_loaded:
            h.nrn_load_dll(mod_folder+"/nrnmech.dll")

        # get the template name
        model_path = myfunctions.posixpth(
            join(model_folder, "ActiveModels", cell_model + '_mod.hoc'))
        f = open(model_path, 'r')
        templatename = myfunctions.get_templatename(f)
        f.close()
        if not hasattr(h, templatename):
            # Load main cell template
            h.load_file(1, model_path)

        cell_parameters = {
            'morphology': myfunctions.posixpth(morph_path),
            'templatefile': model_path,
            'templatename': templatename,
            'templateargs': myfunctions.posixpth(morph_path),
            'v_init': -86,
            'passive': False,
            'dt': dt,  # [ms] Should be a power of 2
            'tstart': -10,  # [ms] Simulation start time
            'tstop': tstop,  # [ms] Simulation end time
            "pt3d": True,
            'nsegs_method': "lambda_f",
            "lambda_f": 100,
        }

    # create cell with parameters in dictionary
    if not active:
        cell = LFPy.Cell(**cell_parameters)
    else:
        cell = LFPy.TemplateCell(**cell_parameters)

    print("L2/3-PC inserted.")

    return cell


def drawRandCellPositions(POPULATION_SIZE, populationParameters):
    """
    Draw & distribute cell positions randomly within a cylinder constraints.

    Parameters
    ----------
    POPULATION_SIZE : int
        number of neurons in the population.
    populationParameters : dict
        dictonary containing the min and max coordinates in the z-axis and the
        radius of the cylinder whithin which the cells have to be distributed.

    Returns
    -------
    cellPositions : ndarray
        x,y,z coordiantes of the cells in the population.

    code taken and adapted from LFPy example example_mpi.py.

    """
    cellPositions = []
    for cellindex in range(POPULATION_SIZE):
        r = np.sqrt(np.random.rand()) * \
            populationParameters['radius']
        theta = np.random.rand() * 2 * np.pi
        x = r * np.sin(theta)
        y = r * np.cos(theta)
        z = np.random.rand() * (populationParameters['zbottom'] -
                                populationParameters['ztop']) + \
            populationParameters['ztop']
        cellPositions.append([x, y, z])
    cellPositions = np.array(cellPositions)

    return cellPositions


def columnPlot(num_neurons, positions, electrodeParameters):
    """
    Plot the positions of the cells in the z-x plane.

    Parameters
    ----------
    num_neurons : int
        number of neurons being plotted.
    positions : ndarray
        x-y-z coordinates of the cells.
    electrodeParameters : 1d array.
        z coordinates of the electrodes in the linear probe.

    Returns
    -------
    None.

    """
    # color = 'r'  # same color for all cells

    fig = plt.figure(figsize=[2.5, 3], dpi=1000)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, num_neurons)))

    ax = fig.add_axes([0.01, 0.0, 1, 1.0],
                      aspect='equal', frameon=False,
                      xticks=[], xticklabels=[],
                      yticks=[], yticklabels=[])
    for cellindex in range(num_neurons):
        "L3-PCs"
        cell = make_cellEyalEtAl2018(
            morphology='2013_03_06_cell03_789_H41_03.ASC')

        cellRotations = [-np.pi/2, -np.pi/7, 0]
        cell.set_rotation(x=cellRotations[0])
        cell.set_rotation(y=cellRotations[1])
        cell.set_rotation(z=cellRotations[2])

        cell.set_pos(x=positions[cellindex, 0],
                     y=positions[cellindex, 1],
                     z=positions[cellindex, 2])

        zips = []

        for x, z in cell.get_idx_polygons():
            zips.append(list(zip(x, z)))

        c = next(color)
        polycol2 = PolyCollection(zips,
                                  edgecolors='none',
                                  facecolors=c,
                                  zorder=positions[cellindex,
                                                   1])
        ax.add_collection(polycol2)

    ax.plot(electrodeParameters['x'],
            electrodeParameters['z'],
            marker='o', markersize=0.5, linewidth=0.5, color='k',
            clip_on=False, zorder=450)

    # fig.savefig('L3_neurons_plot.pdf', dpi=1000)
    # fig.savefig('L3_neurons_plot.jpeg', dpi=1000)


if __name__ == "__main__":
    tstop, dt = 100, 2**-2

    # Define electrode geometry corresponding to a laminar probe:
    a = 50  # location of the first electrode relative to the grey matter / CSF
    # boundary
    electrode_spacing = 100  # inter-electrodes space
    Ne = 17  # number of electrodes inside the cortex
    # z coordinates of the electrodes in microns
    z = -np.mgrid[a:(Ne*electrode_spacing+electrode_spacing):electrode_spacing]

    electrodeParameters = {
        'x': np.zeros(z.size),
        'y': np.zeros(z.size),
        'z': z,
        'sigma': 0.33,            # S/m
        'method': 'pointsource'   # method used to compute the LFP
    }

    num_neurons = 2200  # total number of neurons

    # will draw random cell locations within cylinder constraints:
    populationParameters = {
        'radius': 1500,   # [um] radius of the cortical column
        'ztop': -675,    # [um] upper limit of the neurons position
        'zbottom': -750     # [um] lower limit of the neurons position
    }

    cellpositions = drawRandCellPositions(num_neurons, populationParameters)

    # data = np.load('cellsPops_L3PCs.npz')
    # cellpositions = data['cellPositions_L3']
    columnPlot(num_neurons, cellpositions, electrodeParameters)

    np.savez('cellsPops_L3PCs.npz',
             cellPositions_L3=cellpositions)
