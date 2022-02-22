# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 15:44:37 2021

@author: Beatriz Herrera
"""
from __future__ import division
import os
from os.path import join
import sys
import numpy as np
from scipy import io
import neo
import quantities as pq
import matplotlib.pyplot as plt
import neuron
from neuron import h
import LFPy
from hayetalmodel_active_declarations import active_declarations
from lfpykit import CellGeometry
from elephant import spike_train_generation  # , statistics  # , spectral

SEED = 12

""" Functions """


def run_oneStimulus(cell, cellindex, **kwargs):
    """
    Run simulation using Hay et al. Stimulation Battery.

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    cellindex : int
        index of the cell being simulated.
    **kwargs : dict
        dictonary containig the input parameters given the type of stimulation.

    Returns
    -------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    stimulus : obj, LFPy.pointprocess.PointProcess object
        obj for NEURON point processes representing electrode currents.
    synapse : obj, LFPy.Synapse object
        obj for NEURON point processes representing synaptic inputs.
    isyn : np.array
        array containing the values of the synaptic current at each time point.

    """
    # loading package with stimulus
    neuron.h.load_file("custom_codes.hoc")

    # time vector
    # tot_ntsteps = int(round((cell.tstop) / cell.dt + 1))
    # tvec = np.arange(tot_ntsteps)*cell.dt

    # Stimulation protocol
    if kwargs["stimulationProtocol"] == "SomaticCurrentPulse":
        """Somatic stimulation only -> Current pulse"""
        from haystimbattery_l5 import soma_IClamp_stimuli

        cell, stimulus = soma_IClamp_stimuli(cell)

        # plotting input current
        # if cellindex == 0:
        #     fig = plt.figure(dpi=600)
        #     plt.title('Stimulus')
        #     plt.plot(tvec, stimulus.i, 'k-')
        #     plt.legend('Soma', loc='best')
        #     # plt.show()
        # fig.savefig(join(data_folder, 'Fig' +str(runNumb)+ '_'
        #                  +stimulusType['stimulationProtocol']+ '.png'))
        #     # plt.close(fig)

        synapse = []
        isyn = []

        return cell, stimulus, synapse, isyn

    elif kwargs["stimulationProtocol"] == "BAC_firing":
        """BAC firing Paradigm"""
        from haystimbattery_l5 import BAC_firing

        BAC_params = {
            "experiment_type": kwargs["stimulus_subtype"],
            "proximalpoint": 400,  # [um] proximal apical dendritc point
            "distalpoint": 620,  # [um] distal apical dendritic point
        }
        # inserting stimulus
        cell, synapse, isyn, stimulus, idx_proxDen, idx_distDen = BAC_firing(
            cell, **BAC_params
        )

        # plotting results
        # if cellindex == 0:
        #     fig = plt.figure(dpi=600)
        #     plt.plot(tvec, stimulus.i, 'k-')
        #     plt.plot(tvec, np.zeros(np.size(cell.tvec)), 'b-')
        #     plt.plot(tvec, -isyn, 'r-')
        #     plt.ylabel('Input Current (nA)')
        #     plt.xlabel('Time (ms)')
        # plt.legend(('Soma', 'Proximal Dendrites', 'Distal dendrites'),
        #            loc='best')
        # fig.savefig(join(data_folder, 'Fig' +str(runNumb)+ '_'
        #                  +stimulusType['stimulationProtocol']+ '.png'))
        # plt.close(fig)

        return cell, stimulus, synapse, isyn

    elif kwargs["stimulationProtocol"] == "CriticalFrequency":
        """Critical Frequency Paradigm"""
        from haystimbattery_l5 import critical_frequency

        CFparadigm_params = {
            # Stimulus Type: squarePulse_Train or noisy_current
            "stimulus_type": kwargs["stimulus_subtype"],
            "distalpoint": 620,  # [um] distal dendrites recording location
            "freq": 120,  # [Hz] frequency of the pulses 70 and 120
            # 'stimulus_onset': 10              # [ms] stimulus onset
            "iAmp": kwargs["iAmp"],
        }
        # inserting stimulus
        cell, synapse, isyn, idx_distDen = critical_frequency(
            cell, **CFparadigm_params)
        # plotting results
        # if cellindex == 0:
        #     plt.figure(dpi=600)
        #     plt.plot(tvec, -isyn, 'k-')
        #     # plt.ylim((-0.5,5))
        #     # plt.xlim((280,480))
        #     plt.ylabel('Input Current (nA)')
        #     plt.xlabel('Time (ms)')
        #     plt.legend(('Soma'), loc='best')
        #     fig.savefig(join(data_folder, 'Fig' +str(runNumb)+ '_'
        #                      +stimulusType['stimulus_subtype']+ '.png'))
        #     # plt.close(fig)

        stimulus = []

        return cell, stimulus, synapse, isyn


class Population:
    """Population class."""

    """Copied and adapted to our purposes from LFPy examples: example_mpi.py"""

    def __init__(
        self,
        POPULATION_SIZE,
        cellParameters,
        populationParameters,
        electrodeParameters,
        stimulusType,
        runNumb,
        data_folder,
    ):
        """
        Class initialization.

        POPULATION_SIZE:       int, number of cells
        cellParameters:        dict
        populationParameters:  dict
        electrodeParameters:   dict
        stimulusType:     dict
        inputParameters:       dict
        """
        self.POPULATION_SIZE = POPULATION_SIZE
        self.cellParameters = cellParameters
        self.populationParameters = populationParameters
        self.electrodeParameters = electrodeParameters
        self.stimulusType = stimulusType
        self.runNumb = runNumb
        self.data_folder = data_folder

        # get synaptic times and cell positions, rotations, store in
        # self-object
        data = np.load("cellsPops_L5PCs.npz")
        self.cellPositions = data["cellPositions_L5"]
        self.cellRotations = self.drawRandCellRotations()

    def run(self):
        """Execute the proper simulation and collect simulation results."""
        # produce simulation results
        self.results = self.distribute_cellsims()

        # superimpose local LFPs from all cells, then sum
        self.LFP = []
        if POPULATION_SIZE <= 10:
            self.cell_geoSt = []
            self.somav = []
            self.v_mbp = []

        for key, value in list(self.results.items()):
            self.LFP.append(value["LFP"])  # Ncells x Nelectrodes x Ntimes
            if POPULATION_SIZE <= 10:
                self.cell_geoSt.append(
                    value["cell_geo"]
                )  # Ncells -> one object per cell
                self.somav.append(value["somav"])  # Ncells x Ntimes
                self.v_mbp.append(value["v_mbp"])  # Ncells x Ntimes

        self.LFP = np.array(self.LFP).sum(axis=0)  # summing for the cells
        # -> Nelectrodes x Ntimes
        if POPULATION_SIZE <= 10:
            self.somav = np.array(self.somav)
            self.v_mbp = np.array(self.v_mbp)

    def distribute_cellsims(self):
        """Will run cell simulations."""
        # start unique cell simulation, and store the electrode
        # and cell objects in dicts indexed by cellindex
        results = {}

        for cellindex in range(self.POPULATION_SIZE):
            print("Running cell %d" % cellindex)
            results.update({cellindex: self.cellsim(cellindex)})

        return results

    def make_cellHayModelL5PC(self):
        """
        Create LFPy cell object using Hay et al. 2011 model.

        Returns
        -------
        cell : LFPy.Cell Obj
             object built on top of NEURON representing biological neuron.

        """
        # General simulation parameters
        holding_potential = -80  # [mV] resting membrane potential

        # define cell parameters used as input to cell-class
        model_pth = join("cell_models", "HayModel")
        cellParameters = {
            "morphology": join(model_pth, "morphologies", "cell1.hoc"),
            "v_init": holding_potential,  # initial crossmembrane potential
            "passive": False,  # switch on passive mechs
            "nsegs_method": "lambda_f",  # method for setting number of segments,
            "lambda_f": 100,  # segments are isopotential at this frequency
            # dt of LFP and NEURON simulation.
            "dt": self.cellParameters['dt'],
            "tstart": -250,  # start time, recorders start at t=0
            # stop time, end of the simulation
            "tstop": self.cellParameters['tstop'],
            "custom_code": [join(model_pth, "morphologies", "custom_codes.hoc")],
            "custom_fun": [active_declarations],  # will execute this function
            "custom_fun_args": [{}],
        }

        cell = LFPy.Cell(**cellParameters)
        cell.set_rotation(z=np.pi)

        return cell

    def cellsim(self, cellindex):
        """
        Run cell and LFP simulation procedure.

        Parameters
        ----------
        cellindex : int
            Index of the cell in population being simulated.

        Returns
        -------
        None.

        """
        # Initialize cell instance
        cell = self.make_cellHayModelL5PC(**self.cellParameters)

        # set the position of midpoint in soma
        cell.set_pos(
            x=self.cellPositions[cellindex, 0],
            y=self.cellPositions[cellindex, 1],
            z=self.cellPositions[cellindex, 2],
        )

        # cell, synapse, isyn = make_backgroundNoise(cell=cell)
        cell, stimulus, synapse, isyn = run_oneStimulus(
            cell, cellindex, **self.stimulusType
        )

        # create extracellular electrode object
        electrode = LFPy.RecExtElectrode(cell, **self.electrodeParameters)

        # Parameters for the cell.simulate() call,
        # recording membrane- and syn.-currents
        simulationParameters = {
            'probes': [electrode],
            'rec_imem': True,  # Record Membrane currents during simulation
            'rec_vmem': True,  # record membrane voltage
        }

        # perform NEURON simulation, results saved as attributes in cell
        cell.simulate(**simulationParameters)

        syntype = []
        for i in range(len(cell.synapses)):
            syntype.append(str(cell.synapses[i].syntype))

        zips = []
        for x, z in cell.get_idx_polygons():
            zips.append(list(zip(x, z)))

        cell_geo = CellGeometry(x=cell.x, y=cell.y, z=cell.z, d=cell.d)
        cell_geo.zips = zips

        # extract spike times - soma and dendrites
        signal_somav = neo.core.AnalogSignal(
            cell.somav, units="mV", sampling_rate=(1 / (dt * 1e-3)) * pq.Hz
        )
        signal_dendv = neo.core.AnalogSignal(
            cell.vmem[616, :], units="mV",
            sampling_rate=(1 / (dt * 1e-3)) * pq.Hz
        )

        soma_spkTimes = spike_train_generation.peak_detection(
            signal_somav, threshold=0.0 * pq.mV, sign="above", as_array=False
        )

        dend_spkTimes = spike_train_generation.peak_detection(
            signal_dendv, threshold=-30.0 * pq.mV, sign="above", as_array=False
        )

        saveData = {
            'cell_geo': cell_geo,  # geometry L5 PCs strectched
            # [mV] somatic membrane potatential
            'Vs': np.array(cell.somav),
            # [mV] distal dendrites memberane potential
            'v_mbp': np.array(cell.vmem[616, :]),
            # [ms] time of presynaptic spikes Ncells x Nsynp x Ntskp
            'It': np.array(cell.imem),  # transmembrane currents
            # type of synapse associated with pre_synSpkTimes, same order
            'soma_spkTimes': soma_spkTimes,  # [s] times of somatic APs
            'dend_spkTimes': dend_spkTimes,  # [s]
            'LFP_neuron': electrode.data,  # [mV] lfp produced by the neuron
        }

        # mat files
        io.savemat(join(self.data_folder, 'NeuronsData_r' + str(self.runNumb) +
                        '_n#' +
                        str(cellindex) + '.mat'), saveData)

        # return dict with primary results from simulation
        if POPULATION_SIZE <= 10:
            return {'LFP': electrode.data,
                    'cell_geogeo': cell_geo,
                    "somav": cell.somav,
                    "v_mbp": cell.vmem[616, :],
                    }
        else:
            return {'LFP': electrode.data,
                    }

    def save_simData(self):
        """Save simulations data."""
        saveData = {
            'ze': self.electrodeParameters['z'],  # [um] electrodes position
            'LFP': self.LFP,                     # [mV] LFP values
        }

        # mat files
        io.savemat(join(self.data_folder, 'SimData_r' + str(self.runNumb) +
                   '_PS' + str(self.POPULATION_SIZE) + '.mat'),
                   saveData)

    def plotstuff(self):
        """
        Plot LFPs colormap.

        Returns
        -------
        None.

        """
        fig = plt.figure(dpi=600)
        tvec = np.arange(np.size(self.LFP, axis=1)) * self.cellParameters["dt"]
        im = plt.pcolormesh(
            tvec,
            self.electrodeParameters["z"],
            self.LFP,
            cmap="PRGn",
            vmin=-self.LFP.std() * 3,
            vmax=self.LFP.std() * 3,
            shading="auto",
        )
        plt.axis(plt.axis("tight"))
        cbar = plt.colorbar(im)  # , cax=cax)
        cbar.set_label("LFP (mV)")
        plt.xlabel("time (ms)")
        plt.ylabel(r"$z$ ($\mu$m)")
        plt.xlim(500, 800)

        if "win32" in sys.platform:
            fig.savefig(
                join(
                    self.data_folder,
                    "Fig_run#"
                    + str(self.runNumb)
                    + "_PS"
                    + str(self.POPULATION_SIZE)
                    + ".png",
                )
            )
            #
            plt.show()
        else:
            plt.ioff()
            fig.savefig(
                join(
                    self.data_folder,
                    "Fig_run#"
                    + str(self.runNumb)
                    + "_PS"
                    + str(self.POPULATION_SIZE)
                    + ".png",
                )
            )
            plt.close(fig)


if __name__ == "__main__":

    # load some required neuron-interface files
    h.load_file("stdrun.hoc")
    h.load_file("import3d.hoc")

    dt = 2 ** -7

    cellParameters = {
        "tstop": 800,
        "dt": dt,
    }

    # the number of cells in the population
    POPULATION_SIZE = 2

    num_trials = 1  # number of simulation runs

    for runNumb in np.arange(1, num_trials + 1):

        h("forall delete_section()")

        stimulusType = {
            "HayBattery_active": True,
            "stimulationProtocol": "CriticalFrequency",  # SomaticCurrentPulse,
            # BAC_firing, CriticalFrequency
            "stimulus_subtype": "noisy_current",  # CriticalFrequency:
            # squarePulse_Train or noisy_current
            # BAC_firing: BAP, CaBurst, EPSP, or BAC
            "iAmp": 1.85,
            # 1.9 -> amplitude for supra-CF  # [nA] mean amplitude
            # of the pulse
        }

        # Define electrode geometry corresponding to a laminar probe:
        a = 50
        electrode_spacing = 100
        z = -np.mgrid[a:(17*100+100):electrode_spacing]  # microns
        electrodeParameters = {
            "x": np.zeros(z.size),
            "y": np.zeros(z.size),
            "z": z,
            "sigma": 0.33,  # S/m
            "method": "pointsource",  # method used to compute the LFP
        }

        # will draw random cell locations within cylinder constraints:
        populationParameters = {
            'radius': 1500,   # [um] radius of the cortical column
            'ztop': -1250,    # [um] upper limit of the neurons position 600
            # [um] lower limit of the neurons position 1100
            'zbottom': -1750
        }

        data_folder = join("results_L5PCsPop_CaParadigm_randStimWind",
                           stimulusType["stimulationProtocol"])

        if (
            stimulusType["stimulationProtocol"] == "BAC_firing"
            or stimulusType["stimulationProtocol"] == "CriticalFrequency"
        ):
            data_folder = join(
                data_folder, stimulusType["stimulus_subtype"])

        if not os.path.isdir(data_folder):
            try:
                os.makedirs(data_folder)
            except OSError:
                print("Creation of the directory %s failed" % data_folder)
            else:
                print("Successfully created the directory %s " % data_folder)

        # ############ INITIALIZE POPULATION ##################
        population = Population(
            POPULATION_SIZE,
            cellParameters,
            populationParameters,
            electrodeParameters,
            stimulusType,
            runNumb,
            data_folder,
        )
        population.run()
        population.save_simData()
        population.plotstuff()
