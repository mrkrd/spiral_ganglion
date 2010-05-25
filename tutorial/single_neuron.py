#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

import thorns.nrn as thn
import spiral_ganglion as sg

import neuron
from neuron import h

def main():
    h.celsius = 37
    h.dt = 0.002


    ### Set the auditory nerve fiber
    anf = sg.ANF_Axon(record_voltages=True)
    # anf.set_geometry('straight', x0=250, y0=500, z0=0)
    anf.set_geometry('bent', a=750, b=500, z=0)


    ### Set the electrodes
    el = sg.Electrode()
    el.x = 300
    el.y = 0
    el.z = 0


    ### Set the stimulus
    fs = 10e3
    stim = np.zeros(300)         # 30 ms
    stim[100:105] = -0.1         # 50 us
    stim[105:110] = +0.1         # 50 us
    el.fs = fs
    el.stim = stim



    ### Combine ANF and electrodes
    anf.electrodes = [el]


    ### Run simulation
    tmax = 1000 * len(stim) / fs  # ms
    print "tmax:", tmax, "ms"
    anf.einit()
    neuron.init()
    neuron.run(tmax)


    ### Results
    print "Spikes:", anf.get_spikes()
    plot = thn.plot_voltages(1000/h.dt, anf.get_voltages().T)
    plot.show()


if __name__ == "__main__":
    main()
