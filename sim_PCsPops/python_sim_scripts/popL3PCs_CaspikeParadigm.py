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
import functions
from neuron import h
import LFPy
from lfpykit import CellGeometry
from elephant import spike_train_generation

SEED = 12
np.random.seed(SEED)

""" Functions """


class Population:
    """
    L3 Pyramidal Neuron Population class.

    Create a population (unconnected cells) of 'POPULATION_SIZE' L3 PCs
    described by the Eyal et al. 2018 model, consisting of LFPy.Cell objects.

    Based on the prototype cell population in example_mpi.py from LFPy package
    examples/.
    """

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
        cellParameters:        dict, neuron model and simulation parameters
        populationParameters:  dict, cortical column parameters
        electrodeParameters:   dict, electrodes geometry
        stimulusType:          dict, stimulus specifications
        runNumb:               int, simulated trial number
        data_folder:           str, path to folder where results will be stored
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
        data = np.load("cellsPops_L3PCs.npz")
        self.cellPositions = data["cellPositions_L3"]
        self.cellRotations_z = self.drawRandCellRotations()

    def run(self):
        """Execute the proper simulation and collect simulation results."""
        # produce simulation results
        self.results = self.distribute_cellsims()

        # superimpose local LFPs from all cells, then sum
        self.LFP = []

        for key, value in list(self.results.items()):
            self.LFP.append(value["LFP"])  # Ncells x Nelectrodes x Ntimes

        self.LFP = np.array(self.LFP).sum(axis=0)  # summing for the cells
        # -> Nelectrodes x Ntimes

    def distribute_cellsims(self):
        """Will run cell simulations."""
        # start unique cell simulation, and store the electrode
        # and cell objects in dicts indexed by cellindex
        results = {}

        for cellindex in range(self.POPULATION_SIZE):
            print("Running cell %d" % cellindex)
            results.update({cellindex: self.cellsim(cellindex)})

        return results

    def make_cellEyalEtAl2018(self):
        r"""
        Create LFPy cell object using Eyal et al. 2018 model.

        Copied and addapted from file EEG-master\src\make_data_figure2.py in
        the python package provided in Næss, S., Halnes, G., Hagen, E.,
        Hagler Jr,
        D. J., Dale, A. M., Einevoll, G. T., & Ness, T. V. (2021).
        Biophysically detailed forward modeling of the neural origin of EEG
        and MEG signals. NeuroImage, 225, 117467.

        """
        model_folder = join("cell_models", "EyalEtAl2018")
        morph_path = join(model_folder, "Morphs", self.cellParameters["morphology"])
        add_synapses = False

        if self.cellParameters["cell_model"] is None:
            """Passive Model."""
            # loading mechanisms
            mod_folder = join(model_folder, "mechanisms")
            if not hasattr(h, "NMDA"):
                if "win32" in sys.platform:
                    h.nrn_load_dll(mod_folder + "/nrnmech.dll")
                else:
                    neuron.load_mechanisms(mod_folder)

            cell_parameters = {
                "v_init": -70,
                "morphology": morph_path,
                # S/cm^2, mV
                "passive_parameters": {"g_pas": 1.0 / 30000, "e_pas": -70},
                "Ra": 150,  # Ω cm
                "cm": 1,  # µF/cm^2
                "nsegs_method": "lambda_f",
                "lambda_f": 100,
                "dt": 2 ** -4,  # [ms] Should be a power of 2
                "tstart": -10,  # [ms] Simulation start time
                "tstop": self.cellParameters["tstop"],  # [ms] Simulation
                # end time
                "pt3d": True,
                "passive": True,
            }

            # create cell with parameters in dictionary
            cell = LFPy.Cell(**cell_parameters)

        else:
            """Active Model."""
            # loading mechanisms
            mod_folder = join(model_folder, "ActiveMechanisms")
            if not hasattr(h, "NaTg"):
                if "win32" in sys.platform:
                    h.nrn_load_dll(mod_folder + "/nrnmech.dll")
                else:
                    neuron.load_mechanisms(mod_folder)

            # get the template name
            model_path = functions.posixpth(
                join(
                    model_folder,
                    "ActiveModels",
                    self.cellParameters["cell_model"] + "_mod.hoc",
                )
            )
            f = open(model_path, "r")
            templatename = functions.get_templatename(f)
            f.close()
            if not hasattr(h, templatename):
                # Load main cell template
                h.load_file(1, model_path)

            cell_parameters = {
                "morphology": functions.posixpth(morph_path),
                "templatefile": model_path,
                "templatename": templatename,
                "templateargs": functions.posixpth(morph_path),
                "v_init": -86,
                "passive": False,
                "dt": self.cellParameters["dt"],  # [ms] Should be a power of 2
                "tstart": -10,  # [ms] Simulation start time
                "tstop": self.cellParameters["tstop"],  # [ms] Simulation
                # end time
                "pt3d": True,
                "nsegs_method": "lambda_f",
                "lambda_f": 100,
            }

            # create cell with parameters in dictionary
            cell = LFPy.TemplateCell(**cell_parameters)

        # rotate the morphology
        cellRotations = [-np.pi / 2, -np.pi / 7, 0]
        cell.set_rotation(x=cellRotations[0])
        cell.set_rotation(y=cellRotations[1])
        cell.set_rotation(z=cellRotations[2])

        return cell

    def drawRandCellRotations(self):
        """
        draw random cell rotations around z-axis for all cells in the population.

        Returns
        -------
        ndarray
            rotation angles in radians.

        """
        cellRotations_z = np.random.rand(self.POPULATION_SIZE) * np.pi * 2

        return cellRotations_z

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
        cell = self.make_cellEyalEtAl2018()

        # set the position of midpoint in soma
        cell.set_pos(
            x=self.cellPositions[cellindex, 0],
            y=self.cellPositions[cellindex, 1],
            z=self.cellPositions[cellindex, 2],
        )
        # rotate cell around z-axis
        cell.set_rotation(z=self.cellRotations_z[cellindex])

        # loading package with stimulus
        neuron.h.load_file("custom_codes.hoc")
        from haystimbattery_l3 import critical_frequency

        CFparadigm_params = {
            # Stimulus Type: squarePulse_Train or noisy_current
            "stimulus_type": self.stimulusType["stimulus_subtype"],
            "distalpoint": 258.7071228,  # [um] distal dendrites
            # recording location
            "freq": 120,  # [Hz] frequency of the pulses 70 and 120
            "iAmp": self.stimulusType["iAmp"],
        }
        # inserting stimulus
        cell, synapse, isyn, idx_distDen = critical_frequency(cell, **CFparadigm_params)

        # create extracellular electrode object
        electrode = LFPy.RecExtElectrode(cell, **self.electrodeParameters)

        # Parameters for the cell.simulate() call,
        # recording membrane- and syn.-currents
        simulationParameters = {
            "probes": [electrode],
            "rec_imem": True,  # Record Membrane currents during simulation
            "rec_vmem": True,  # record membrane voltage
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
            cell.vmem[612, :], units="mV", sampling_rate=(1 / (dt * 1e-3)) * pq.Hz
        )

        soma_spkTimes = spike_train_generation.peak_detection(
            signal_somav, threshold=0.0 * pq.mV, sign="above", as_array=False
        )

        dend_spkTimes = spike_train_generation.peak_detection(
            signal_dendv, threshold=-30.0 * pq.mV, sign="above", as_array=False
        )

        saveData = {
            "cell_geo": cell_geo,  # geometry L3 PCs strectched
            # [mV] somatic membrane potatential
            "Vs": np.array(cell.somav),
            # [mV] distal dendrites memberane potential
            "v_mbp": np.array(cell.vmem[612, :]),
            # [ms] time of presynaptic spikes Ncells x Nsynp x Ntskp
            "It": np.array(cell.imem),  # transmembrane currents
            "soma_spkTimes": soma_spkTimes,  # [s] times of somatic APs
            "dend_spkTimes": dend_spkTimes,  # [s]
            "LFP_neuron": electrode.data,  # [mV] lfp produced by the neuron
        }

        # mat files
        io.savemat(
            join(
                self.data_folder,
                "NeuronsData_r" + str(self.runNumb) + "_n#" + str(cellindex) + ".mat",
            ),
            saveData,
        )

        # return dict with primary results from simulation
        return {
            "LFP": electrode.data,
        }

    def save_simData(self):
        """Save simulations data."""
        saveData = {
            "ze": self.electrodeParameters["z"],  # [um] electrodes position
            "LFP": self.LFP,  # [mV] LFP values
        }

        # mat files
        io.savemat(
            join(
                self.data_folder,
                "SimData_r"
                + str(self.runNumb)
                + "_PS"
                + str(self.POPULATION_SIZE)
                + ".mat",
            ),
            saveData,
        )

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
        "cell_model": "cell0603_08_model_602",
        "morphology": "2013_03_06_cell03_789_H41_03.ASC",
    }

    # the number of cells in the population
    POPULATION_SIZE = 2  # 2200

    stimulusType = {
        "stimulus_subtype": "noisy_current",  # CriticalFrequency:
        # squarePulse_Train or noisy_current
        # BAC_firing: BAP, CaBurst, EPSP, or BAC
        "iAmp": 1.9,  # [nA] mean amplitude
        # of the pulse
        # 1.9 -> amplitude for supra
        # 1.85 -> threshold input
    }

    # Define electrode geometry corresponding to a laminar probe:
    a = 50
    electrode_spacing = 100
    z = -np.mgrid[a : (17 * 100 + 100) : electrode_spacing]  # microns
    electrodeParameters = {
        "x": np.zeros(z.size),
        "y": np.zeros(z.size),
        "z": z,
        "sigma": 0.33,  # S/m
        "method": "pointsource",  # method used to compute the LFP
    }

    # will draw random cell locations within cylinder constraints:
    populationParameters = {
        "radius": 1500,  # [um] radius of the cortical column
        "zmin": -675,  # [um] upper limit of the neurons position
        "zmax": -750,  # [um] lower limit of the neurons position
    }

    data_folder = join(
        "CSDtoEEG_paper_sim", "sim_L3PCs", stimulusType["stimulus_subtype"],
    )

    if not os.path.isdir(data_folder):
        try:
            os.makedirs(data_folder)
        except OSError:
            print("Creation of the directory %s failed" % data_folder)
        else:
            print("Successfully created the directory %s " % data_folder)

    for runNumb in np.arange(1, 2):

        h("forall delete_section()")

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
