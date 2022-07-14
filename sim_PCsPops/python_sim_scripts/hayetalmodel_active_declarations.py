# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 10:25:51 2021

@author: Beatriz Herrera
"""
from __future__ import division
from os.path import join
import neuron
from neuron import h
import sys

""" Functions """


def biophys_active():
    r"""
    L5PCbiophys3 - Hay et al 2011 definition.

    Copied and addapted from file hay/hay_active_declarations.py in the
    python package provided in Miceli, Ness, Einevoll, Schubert (2017).
    Impedance Spectrum in Cortical Tissue: Implications for Propagation of LFP
    Signals on the Microscopic Level.  Eneuro 4:1-15.

    Parameters
    ----------
    None.

    Returns
    -------
    None.

    """
    # loading mechanisms
    model_pth = join("cell_models", "HayModel")
    pth = join(model_pth, "mod")
    if not hasattr(h, "Ca_LVAst"):
        if "win32" in sys.platform:
            h.nrn_load_dll(pth + "/nrnmech.dll")
        else:
            neuron.load_mechanisms(pth)

    for sec in neuron.h.allsec():
        sec.insert("pas")
        sec.cm = 1.0
        sec.Ra = 100.0
        sec.e_pas = -90.0

    for sec in neuron.h.soma:
        sec.insert("Ca_LVAst")
        sec.insert("Ca_HVA")
        sec.insert("SKv3_1")
        sec.insert("SK_E2")
        sec.insert("K_Tst")
        sec.insert("K_Pst")
        sec.insert("Nap_Et2")
        sec.insert("NaTa_t")
        sec.insert("CaDynamics_E2")
        sec.insert("Ih")
        sec.ek = -85
        sec.ena = 50
        sec.gIhbar_Ih = 0.0002
        sec.g_pas = 0.0000338
        sec.decay_CaDynamics_E2 = 460.0
        sec.gamma_CaDynamics_E2 = 0.000501
        sec.gCa_LVAstbar_Ca_LVAst = 0.00343
        sec.gCa_HVAbar_Ca_HVA = 0.000992
        sec.gSKv3_1bar_SKv3_1 = 0.693
        sec.gSK_E2bar_SK_E2 = 0.0441
        sec.gK_Tstbar_K_Tst = 0.0812
        sec.gK_Pstbar_K_Pst = 0.00223
        sec.gNap_Et2bar_Nap_Et2 = 0.00172
        sec.gNaTa_tbar_NaTa_t = 2.04

    for sec in neuron.h.apic:
        sec.cm = 2
        sec.insert("Ih")
        sec.insert("SK_E2")
        sec.insert("Ca_LVAst")
        sec.insert("Ca_HVA")
        sec.insert("SKv3_1")
        sec.insert("NaTa_t")
        sec.insert("Im")
        sec.insert("CaDynamics_E2")
        sec.ek = -85
        sec.ena = 50
        sec.decay_CaDynamics_E2 = 122
        sec.gamma_CaDynamics_E2 = 0.000509
        sec.gSK_E2bar_SK_E2 = 0.0012
        sec.gSKv3_1bar_SKv3_1 = 0.000261
        sec.gNaTa_tbar_NaTa_t = 0.0213
        sec.gImbar_Im = 0.0000675
        sec.g_pas = 0.0000589

    h.distribute_channels("apic", "gIhbar_Ih", 2, -0.8696, 3.6161, 0.0, 2.087, 0.0002)
    h.distribute_channels(
        "apic", "gCa_LVAstbar_Ca_LVAst", 3, 1.0, 0.010, 685.0, 885.0, 0.0187
    )
    h.distribute_channels(
        "apic", "gCa_HVAbar_Ca_HVA", 3, 1.0, 0.10, 685.00, 885.0, 0.000555
    )

    for sec in neuron.h.dend:
        sec.cm = 2
        sec.insert("Ih")
        sec.gIhbar_Ih = 0.0002
        sec.g_pas = 0.0000467

    for sec in neuron.h.axon:
        sec.g_pas = 0.0000325

    print("L5-PC inserted.")


def active_declarations(cell):
    """
    Set active conductances for Hay model 2011.

    Parameters
    ----------
    cell : LFPy.Cell Obj 


    Returns
    -------
    None.

    Copied and addapted from file hay/hay_active_declarations.py in the python
    package provided in Miceli, Ness, Einevoll, Schubert (2017) Impedance
    Spectrum in Cortical Tissue: Implications for Propagation of LFP Signals
    on the Microscopic Level. Eneuro 4:1-15.

    """
    h.delete_axon()
    h.geom_nseg()
    h.define_shape()
    exec("biophys_active()")


def get_index(cell, seg_x, HCell_sec):
    """
    Return the cell ind associated with the section=HCell_sec and seg.x=seg_x.

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    seg_x : int
        value of seg_x.
    HCell_sec : hoc.HocObject
        reference to the NEURON object representing the section of a neuron.

    Returns
    -------
    i : int
        index of the compartment associated with the section=HCell_sec
        and seg.x=seg_x.

    """
    i = 0
    names = []
    for sec in cell.allseclist:
        # print(sec.name())
        for seg in sec:
            names.append(sec(seg.x))
            if sec.name() == HCell_sec.name() and (
                h.distance(sec(seg.x)) == h.distance(HCell_sec(seg_x))
            ):
                # print('Sec %s = %s' % (sec.name(), HCell_sec.name()))
                # print('Sec %d = %d' % (h.distance(sec(seg.x)),
                # h.distance(HCell_sec(seg_x))))
                break
            i += 1
        if sec.name() == HCell_sec.name():
            break
    return i
