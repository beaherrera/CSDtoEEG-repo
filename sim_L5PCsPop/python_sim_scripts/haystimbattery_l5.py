# -*- coding: utf-8 -*-

"""
Created on Sat May 23 14:09:38 2020

@author: Beatriz Herrera
"""
from __future__ import division
import numpy as np
import LFPy
import random
import neuron

nrn = neuron.h


def soma_IClamp_stimuli(cell):
    """
    Create a current pulse input for soma stimuli.

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.

    Returns
    -------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    stimulus : obj, LFPy.pointprocess.PointProcess object
        obj for NEURON point processes representing electrode currents.

    """
    PointProcParams = {
        "idx": 0,
        "record_current": True,
        "pptype": "IClamp",
        "amp": 0.793,
        "dur": 150,
        "delay": 600,
    }
    stimulus = LFPy.StimIntElectrode(cell, **PointProcParams)

    return cell, stimulus


def myexp(x):
    """
    Return the exp(x) if x < -100, and 0 otherwise.

    Parameters
    ----------
    x : float

    Returns
    -------
    myexp : float

    """
    if x < -100:
        myexp = 0
    else:
        myexp = np.exp(x)

    return myexp


def _make_epsp(cell, onset=0, tau0=0.2, tau1=3.0, imax=0):
    """
    Create an excitatory post-synaptic (EPSP) like input current.

    Based on Schaefer et al. 2003.

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    onset : float, optional
        [ms] stimulus onset. The default is 0.
    tau0 : float, optional
        [ms] rise time constant of double-exponential function.
        The default is 0.2.
    tau1 : float, optional
        [ms] decay time constant of double-exponential function.
        The default is 3.0.
    imax : float, optional
        [nA] maximum current amplitude. The default is 0.

    Returns
    -------
    It : ndarray
        [nA] input current (np.ndarray).

    """
    It = []  # initializing
    tot_ntsteps = int(round((cell.tstop) / cell.dt + 1))
    tvec = np.arange(tot_ntsteps) * cell.dt
    tpeak = tau0 * tau1 * np.log(tau0 / tau1) / (tau0 - tau1)
    adjust = 1 / ((1 - myexp(-tpeak / tau0)) - (1 - myexp(-tpeak / tau1)))
    amp = adjust * imax

    for t in tvec:
        if t < onset:
            value = 0
            It = np.append(It, value, axis=None)
        else:
            a0 = 1 - myexp(-(t - onset) / tau0)
            a1 = 1 - myexp(-(t - onset) / tau1)
            It = np.append(It, (-amp * (a0 - a1)), axis=None)

    return It


def insert_epsp_stimuli(cell, input_idx, input_seg, **kwargs):
    """
    Insert ESPS-like input current.

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    input_idx : int
        compartment index where input will be inserted.
    input_seg : float
        segment of the section associated with the compartment index input_idx
        where the input will be inserted.
    **kwargs : dict
        dictonary containing input parameters depending on the input type.

    Returns
    -------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    syn : TYPE
        DESCRIPTION.
    isyn : TYPE
        DESCRIPTION.

    """
    onset = kwargs["onset"]
    tau0 = kwargs["tau0"]
    tau1 = kwargs["tau1"]
    imax = kwargs["imax"]

    input_array = _make_epsp(cell, onset, tau0, tau1, imax)
    isyn = neuron.h.Vector(input_array)

    # print("Input inserted in ", HCell.apic[input_idx].name())
    syn = nrn.ISyn(input_seg, sec=nrn.apic[input_idx])
    # print("Dist: ", nrn.distance(HCell.apic[input_idx](input_seg)))
    syn.dur = 1e9
    syn.delay = 0  # cell.tstart

    isyn.play(syn._ref_amp, cell.dt)

    return cell, syn, isyn


def _get_apic_idx_outputMatrix(cell, int_idx, int_seg):
    """
    Compute the idx of the apical segment, apic[int_idx](int_seg), in cell.

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    int_idx : int
        apic sesction number in the NEURON.hoc object.
    int_seg : int
        seg_x value within apic[int_idx].

    Returns
    -------
    i : int
        compartment number associated with apic[int_idx](int_seg) in LFPy.cell
        object.

    """
    i = 0
    for sec in cell.allseclist:
        for seg in sec:
            if sec.name() == nrn.apic[int_idx].name() and nrn.distance(
                seg.x
            ) == nrn.distance(nrn.apic[int_idx](int_seg)):
                # print('Sec %s = %s' % (sec.name(), nrn.apic[int_idx].name()))
                # print('Sec %d = %d' % (nrn.distance(seg.x),
                #                    nrn.distance(nrn.apic[int_idx](int_seg))))
                break
            i += 1
        if sec.name() == nrn.apic[int_idx].name():
            break

    # print('idx_matrix=%d' % i)

    return i


def BAC_firing(cell, **kwargs):
    """
    BAC-firing protocol.

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    **kwargs : dict
        dictonary containing the input parameters depending on the input type.

    Returns
    -------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    syn : TYPE
        DESCRIPTION.
    isyn : TYPE
        DESCRIPTION.
    stimulus : obj, LFPy.pointprocess.PointProcess object
        obj for NEURON point processes representing electrode currents.
    idx_proxDen : int
        index of the proximal dendrites compartment.
    idx_distDen : int
        index of the distal trunk compartment or main bifurcation point.

    Copied and adapted some needed functions to work with LFPy from file
    simulationcode/BAC_firing.hoc

    """
    #   experiment_type: BAP, CaBurst, EPSP, or BAC
    experiment_type = (
        kwargs["experiment_type"] if "experiment_type" in kwargs else "BAC"
    )

    BACdt = 5  # [ms] pulse width

    """ somatic pulse settings """
    squareAmp = 1.5  # [nA] pulse current amplitude
    pulseDur = 5  # [ms] duration of the current pulse
    # 295
    # [ms] delay of the stimulation -> original value in Hay et al. 2011: 295ms
    delayS = 100

    """ EPSP settings """
    risetau = 0.5  # [ms] rise time constant
    decaytau = 5  # [ms] decay time constant
    Imax = 0.5  # [nA] maximum current amplitude

    # [um] proximal apical dendritc point
    proximalpoint = kwargs["proximalpoint"]
    distalpoint = kwargs["distalpoint"]  # [um] distal apical dendritic point

    if experiment_type == "BAP":
        somastimamp = squareAmp  # somatic stimuli amplitude
        EPSPamp = 0  # EPSP amplitude (dendritic stimuli)
    if experiment_type == "CaBurst":
        somastimamp = 0  # somatic stimuli amplitude
        EPSPamp = Imax * 5  # EPSP amplitude (dendritic stimuli)
    if experiment_type == "BAC":
        somastimamp = squareAmp  # somatic stimuli amplitude
        EPSPamp = Imax  # EPSP amplitude (dendritic stimuli)
    if experiment_type == "EPSP":
        somastimamp = 0  # somatic stimuli amplitude
        EPSPamp = Imax  # EPSP amplitude (dendritic stimuli)

    # -- somatic pulse
    PointProcParams = {
        "idx": 0,
        "record_current": True,
        "pptype": "IClamp",
        "amp": somastimamp,
        "dur": pulseDur,
        "delay": delayS,
    }

    stimulus = LFPy.StimIntElectrode(cell, **PointProcParams)

    # -- Dendritic EPSP-like current
    # HCell = cell.template # getattr(neuron.h, cell.name)[0]
    nrn(
        """
        objref sl //synaptic locations list

        sl = new List()
        double siteVec[2]

        sl = locateSites("apic","""
        + str(distalpoint)
        + """)

        maxdiam = 0
        for(i=0;i<sl.count();i+=1){
          dd1 = sl.o[i].x[1]
          dd = apic[sl.o[i].x[0]].diam(dd1)
          if (dd > maxdiam) {
            j = i
            //print j
            maxdiam = dd
          }
        }

        siteVec[0] = sl.o[j].x[0]
        siteVec[1] = sl.o[j].x[1]
    """
    )
    # stimulation location
    input_idx = int(nrn.siteVec[0])
    input_seg = nrn.siteVec[1]
    # stimulus parameters
    EPSP_params = {
        "onset": delayS + BACdt,  # original: substitute delayS for 295
        "tau0": risetau,
        "tau1": decaytau,
        "imax": EPSPamp,
    }
    cell, syn, isyn = insert_epsp_stimuli(
        cell, input_idx, input_seg, **EPSP_params)

    # plt.figure()
    # tot_ntsteps = int(round((cell.tstop) / cell.dt + 1))
    # tvec = np.arange(tot_ntsteps) * cell.dt
    # plt.plot(tvec, -isyn)
    # plt.show()

    # Proximal dendrites - recording location
    nrn(
        """
    sl = new List()

    sl = locateSites("apic","""
        + str(proximalpoint)
        + """)

    maxdiam = 0
    for(i=0;i<sl.count();i+=1){
      dd1 = sl.o[i].x[1]
      dd = apic[sl.o[i].x[0]].diam(dd1)
      if (dd > maxdiam) {
        j = i
        maxdiam = dd
      }
    }

    siteVec[0] = sl.o[j].x[0]
    siteVec[1] = sl.o[j].x[1]
    """
    )
    record_idx = int(nrn.siteVec[0])
    record_seg = nrn.siteVec[1]

    # getting the idx of the apical dendritic locations in the output matrix
    idx_proxDen = _get_apic_idx_outputMatrix(cell, record_idx, record_seg)
    idx_distDen = _get_apic_idx_outputMatrix(cell, input_idx, input_seg)

    return cell, syn, isyn, stimulus, idx_proxDen, idx_distDen


def _make_squarePulse_Train(cell, pulseNum=5, freq=120,
                            pulseDur=5, pulseAmp=1.99):
    """
    Create a Square Pulses input current.

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    pulseNum : int, optional
        Number of square pulses. The default is 5.
    freq : float, optional
        [Hz] frequency of the pulses. The default is 120.
    pulseDur : float, optional
        [ms] duration of the pulse. The default is 5.
    pulseAmp : float, optional
        [nA] pulse amplitude. The default is 1.99.

    Returns
    -------
    It : np.ndarray
        [nA] input current (np.ndarray).

    """
    tot_ntsteps = int(round(cell.tstop / cell.dt + 1))

    It = np.zeros(tot_ntsteps)

    for i in np.arange(pulseNum):
        pulseStart = cell.tstop / 4 + i * 1000 / freq
        pulseEnd = pulseStart + pulseDur
        tStart_idx = int(round(pulseStart / cell.dt + 1))
        tEnd_idx = int(round(pulseEnd / cell.dt + 1))

        It[tStart_idx:tEnd_idx:1] = -pulseAmp * \
            np.ones((tEnd_idx - tStart_idx))

        # print('Pulse dur: %d' % ((tEnd_idx - tStart_idx)*cell.dt))
        # print('stim_start: %d \n stim_end: %d \n' % (pulseStart, pulseEnd))

    # tvec = np.arange(tot_ntsteps) * cell.dt
    # plt.plot(tvec, I, 'b-')
    # # plt.xlim((200, 300))
    # plt.show()

    return It


def Inoise(dt, mu, sigma, tau, Ipre):
    """
    Generate a noisy current drawn from a Gaussian distribution.

    Parameters
    ----------
    dt : float
        time step in ms.
    mu : float
        mean of the noisy current.
    sigma : float
        standard deviation of the Gaussian distribution.
    tau : float
        correlation length.
    Ipre : float
        current at t-dt, previos time, in nA.

    Returns
    -------
    In : np.array
        Noisy input current.

    """
    In = (
        Ipre
        + (mu - Ipre) * dt / tau
        + sigma * np.random.randn() * np.sqrt(2 * dt / tau)
    )

    return In


def _make_noisyCurrentPulse(cell, iAmp=7, isigma=0.2, ionset=10, idur=20,
                            itau=3):
    """
    Noisy Current Pulse - Optogenetic-like stimulation paradigm.

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    iAmp : float, optional
        [nA] mean amplitude of the pulse. The default is 7.
    isigma : float, optional
        [nA] standard deviation of the noise. The default is 0.2.
    ionset : float, optional
        [ms] stimulation onset. The default is 10.
    idur : float, optional
        [ms] correlation length. The default is 20.
    itau : float, optional
        [ms] duration of the stimulation. The default is 3.

    Returns
    -------
    In : np.array
        Noisy input current.

    """
    tot_ntsteps = int(round(cell.tstop / cell.dt + 1))

    In = np.zeros(tot_ntsteps)

    iend = ionset + idur
    idxStart = int(round(ionset / cell.dt + 1))
    idxEnd = int(round(iend / cell.dt + 1))
    In[idxStart:idxEnd:1] = iAmp + isigma * \
        np.random.randn(len(In[idxStart:idxEnd:1]))

    In = -In

    return In


def insert_CF_stimuli(cell, input_idx, input_seg, stimulusType, **kwargs):
    """
    Insert CF stimuli depending on the chosen stimulus type.

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    input_idx : int
        Idx of the compartment in LFPy.Cell object where the input will be
        applied.
    input_seg : float
        Segment of the section onto which the input will be inserted.
    stimulusType : str
        Str containing the stimulus type: square pulses or square noisy input
        current.
    **kwargs : dict
        Dictonary with stimulus parameters.

    Returns
    -------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    syn : ndarray
        Input current [nA].
    isyn : TYPE
        DESCRIPTION.

    """
    if stimulusType == "squarePulse_Train":
        """Parameters"""
        pulseNum = kwargs["pulseNum"]  # number of square pulses
        freq = kwargs["freq"]  # [Hz] frequency of the pulses
        pulseDur = kwargs["pulseDur"]  # [ms] duration of the stimulation
        pulseAmp = kwargs["pulseAmp"]  # [nA] amplitude of the square pulses

        """Creating the stimulus"""
        input_array = _make_squarePulse_Train(
            cell, pulseNum, freq, pulseDur, pulseAmp)

    if stimulusType == "noisy_current":
        iAmp = kwargs["iAmp"]  # [nA] mean amplitude of the pulse
        isigma = kwargs["isigma"]  # [nA] standard deviation of the noise
        ionset = kwargs["ionset"]  # [ms] stimulation onset
        idur = kwargs["idur"]  # [ms] duration of the stimulation

        """Creating the stimulus"""
        if "itau" in kwargs:
            itau = kwargs["itau"]  # [ms] correlation length
            input_array = _make_noisyCurrentPulse(
                cell, iAmp, isigma, ionset, idur, itau
            )
        else:
            input_array = _make_noisyCurrentPulse(
                cell, iAmp, isigma, ionset, idur)

    isyn = neuron.h.Vector(input_array)

    # print("Input inserted in ", HCell.soma[input_idx].name())
    syn = nrn.ISyn(input_seg, sec=nrn.soma[input_idx])
    # print("Dist: ", nrn.distance(nrn.soma[input_idx](input_seg)))
    syn.dur = 1e9
    syn.delay = 0  # cell.tstart

    isyn.play(syn._ref_amp, cell.dt)

    return cell, syn, isyn


def critical_frequency(cell, **kwargs):
    """
    Critical Frequency protocol.

    Protocol 1: standard square train pulses.
    Protocol 2: optogenetic-like input current.

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    **kwargs : dict
        Dictonary with stimulus parameters.

    Returns
    -------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    syn : TYPE
        DESCRIPTION.
    isyn : TYPE
        DESCRIPTION.
    idx_distalDen : int
        Index of the distal apical compartment where input was applied.

    """
    # type of stimulation that will be used
    stimulusType = kwargs["stimulus_type"]

    if stimulusType == "squarePulse_Train":
        """Parameters"""
        stim_Parameters = {
            "pulseNum": 5,  # number of square pulses
            "freq": kwargs["freq"],  # [Hz] frequency of the pulses 70 and 120
            "pulseDur": 5,  # [ms] duration of the squared current pulse
            # [nA] amplitude of the squared current pulses 1.99
            "pulseAmp": 1.99,
        }
    # Note: this pulse amplitude yields the correct behaviour in 64
    # bit NEURON environment.
    # In 32 bit NEURON environements, due to difference in float precision,
    # this amplitude may need to be
    # modified slightly (amps = 1.94 nA).
    if stimulusType == "noisy_current":
        stim_Parameters = {
            "iAmp": kwargs["iAmp"],  # 1.9,  # [nA] mean amplitude
            # of the pulse 1.89
            "isigma": 0.3,  # [nA] standard deviation of the noise 0.3
            # kwargs['stimulus_onset'],     # [ms] stimulation onset 310
            "ionset": 600 + random.randint(0, 30),
            "idur": 30,  # [ms] duration of the stimulation
        }

    # estimating the recording location
    distalpoint = kwargs["distalpoint"]
    # HCell = cell.template # getattr(neuron.h, cell.name)[0]
    nrn(
        """
        // Copied and adapted some needed functions to work with LFPy from file
        // simulationcode/BAC_firing.hoc
        objref sl //synaptic locations list

        sl = new List()
        double siteVec[2]

        sl = locateSites("apic","""
        + str(distalpoint)
        + """)

        maxdiam = 0
        for(i=0;i<sl.count();i+=1){
          dd1 = sl.o[i].x[1]
          dd = apic[sl.o[i].x[0]].diam(dd1)
          if (dd > maxdiam) {
            j = i
            //print j
            maxdiam = dd
          }
        }

        siteVec[0] = sl.o[j].x[0]
        siteVec[1] = sl.o[j].x[1]
    """
    )
    record_idx = int(nrn.siteVec[0])
    record_seg = nrn.siteVec[1]

    # inserting the stimulus
    input_idx = 0  # soma
    input_seg = 0.5  # middle of the soma
    cell, syn, isyn = insert_CF_stimuli(
        cell, input_idx, input_seg, stimulusType, **stim_Parameters
    )

    # getting the idx of the apical dendritic locations in the output matrix
    idx_distalDen = _get_apic_idx_outputMatrix(cell, record_idx, record_seg)

    return cell, syn, isyn, idx_distalDen
