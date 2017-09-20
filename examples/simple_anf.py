#!/usr/bin/env python

from __future__ import division, print_function, absolute_import
from __future__ import unicode_literals

import numpy as np
import matplotlib.pyplot as plt

import spiral_ganglion as sg

from neuron import h


def main():

    tmax = 20e-3
    fs = 500e3

    sg.set_fs(fs)
    sg.set_celsius(37)

    # Acoustic stimmulation
    anf_ac = sg.ANF_Axon(record_voltages=True)

    h.topology()
    h.psection(
        sec=anf_ac.get_sections()[0]
    )

    anf_ac.vesicles = [2e-3, 5e-3]

    # Electrical stimulation

    # set-up ANF
    anf_el = sg.ANF_Axon(record_voltages=True)
    # anf_el.set_geometry('straight', x0=250e-6, y0=500e-6, z0=0)
    anf_el.set_geometry('bent', a=750e-6, b=500e-6, z=0)

    # set-up electrode
    el = sg.Electrode()
    el.z = 0

    stim = np.zeros(int(tmax * fs))
    stim[int(tmax/3*fs):int((tmax/3+1e-3)*fs)] = -0.2e-3  # (A)

    el.fs = fs
    el.stim = stim

    # `connect' ANF and electrode
    anf_el.electrodes = [el]

    # Run
    sg.run(
        tmax,
        [anf_ac, anf_el]
    )

    # Plot
    print(anf_ac.get_spikes(), anf_el.get_spikes())

    sg.plot_vectors(
        anf_ac.get_voltages()[:, 0:6], fs
    )

    sg.plot_vectors(
        anf_el.get_voltages()[:, 0:6], fs
    )

    sg.plot_geometry(
        [anf_el, el]
    )

    plt.show()


if __name__ == "__main__":
    main()
