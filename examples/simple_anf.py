#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np
import matplotlib.pyplot as plt

import spiral_ganglion as sg
from spiral_ganglion.anf import _plot_voltages

import neuron
from neuron import h


def main():

    sg.set_fs(500e3)
    sg.set_celsius(37)

    anf = sg.ANF_Axon(record_voltages=True)

    h.topology()
    h.psection(sec=anf.sections['sec'][0])

    anf.vesicles = [2, 5]

    neuron.init()
    anf.ainit()
    neuron.run(10)

    _plot_voltages( anf.get_voltages()[:,0:6] )

    print( "Spikes:", anf.get_spikes())


    print()
    print( "===============================================")
    print( "Electrical stimulation")

    # set ANF
    anf = sg.ANF_Axon(record_voltages=True)
    # anf.set_geometry('straight', x0=250e-6, y0=500e-6, z0=0)
    anf.set_geometry('bent', a=750e-6, b=500e-6, z=0)

    # set electrode
    el = sg.Electrode()
    el.z = 0

    fs = 200e3                  # [Hz]
    stim = np.zeros(20e-3 * fs)
    stim[np.round(5e-3*fs):np.round(6e-3*fs)] = -0.2e-3 # [A]

    el.fs = fs
    el.stim = stim

    # `connect' ANF and electrode
    anf.electrodes = [el]

    # run
    anf.einit()
    neuron.init()
    neuron.run(len(stim) * h.dt)




    ### Conductances
    capacity = 0.0714e-12

    print()
    print( "g_Na ", anf.sections['sec'][0].gnabar_na_schwarz1987)
    print( "g_Kv ", anf.sections['sec'][0].gkbar_k_schwarz1987)
    print( "g_Klt", anf.sections['sec'][0].gkltbar_klt_rothman2003)
    print( "g_h  ", anf.sections['sec'][0].ghbar_ih_rothman2003)
    print( "g_pas", anf.sections['sec'][0].g_pas)
    print()





    # plot
    print( anf.get_spikes())
    _plot_voltages(anf.get_voltages()[:,0:6])


    sg.plot_geometry( [anf,el] )
    plt.show()


if __name__ == "__main__":
    main()
