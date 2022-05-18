# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 15:54:03 2021

@author: Beatriz Herrera
"""
from __future__ import division
from os.path import join
import numpy as np
from hayetalmodel_active_declarations import active_declarations
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import LFPy

np.random.seed(0)


''' Functions '''


def make_cellHayModelL5PC(tstop, dt):
    """
    Create LFPy cell object using Hay et al. 2011 model.

    Parameters
    ----------
    tstop : TYPE
        DESCRIPTION.
    dt : TYPE
        DESCRIPTION.

    Returns
    -------
    cell : TYPE
        DESCRIPTION.

    """
    # General simulation parameters
    holding_potential = -80  # [mV] resting membrane potential
    tstart = -250  # [ms]

    # define cell parameters used as input to cell-class
    model_pth = join("cell_models", "HayModel")
    cellParameters = {
        "morphology": join(model_pth, "morphologies", "cell1.hoc"),
        "v_init": holding_potential,  # initial crossmembrane potential
        "passive": False,  # switch on passive mechs
        "nsegs_method": "lambda_f",  # method for setting number of segments,
        "lambda_f": 100,  # segments are isopotential at this frequency
        "dt": dt,  # dt of LFP and NEURON simulation.
        "tstart": tstart,  # start time, recorders start at t=0
        "tstop": tstop,  # stop time, end of the simulation
        "custom_code": [join(model_pth, "morphologies", "custom_codes.hoc")],
        "custom_fun": [active_declarations],  # will execute this function
        "custom_fun_args": [{"hold_potential": holding_potential}],
    }

    cell = LFPy.Cell(**cellParameters)
    cell.set_rotation(z=np.pi)

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
        "L5-PCs"
        cell = make_cellHayModelL5PC(tstop, dt)

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

    fig.savefig('neurons_plot.pdf', dpi=1000)
    fig.savefig('neurons_plot.jpeg', dpi=1000)


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

    num_neurons = 1000  # total number of neurons

    # will draw random cell locations within cylinder constraints:
    populationParameters = {
        'radius': 1500,   # [um] radius of the cortical column
        'ztop': -1250,    # [um] upper limit of the neurons position
        'zbottom': -1750     # [um] lower limit of the neurons position
    }

    # generate neurons location
    # cellpositions = drawRandCellPositions(num_neurons, populationParameters)

    # load neurons location if already generated and you're just plotting,
    # otherwise: comment next lines and uncomment call to
    # drawRandCellPositions() above
    data = np.load('cellsPops_L5PCs.npz')
    cellpositions = data['cellPositions_L5']
    # plot neurons location
    columnPlot(num_neurons, cellpositions, electrodeParameters)

    # save neurons location
    # np.savez('cellsPops_L5PCs.npz',
    #          cellPositions_L5=cellpositions)
