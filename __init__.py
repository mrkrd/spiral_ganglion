#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

__author__ = "Marek Rudnicki"

import numpy as np

import neuron
from neuron import h

import os
lib_dir = os.path.dirname(__file__)
neuron.load_mechanisms(lib_dir)

from anf import (
    ANF_Axon,
    plot_geometry,
    plot_vectors
)

from electrodes import Electrode

from ci import (
    run_ci_simulation,
    make_anf_electrode,
    find_threshold
)


import signals
import ci


def run(
        tmax,
        anfs
):
    for anf in anfs:
        anf.einit()

    neuron.init()

    for anf in anfs:
        anf.ainit()

    neuron.run(tmax * 1e3)      # s -> ms



def set_celsius(celsius):
    h.celsius = celsius


def set_fs(fs):
    h.dt = 1e3/fs


def show():
    import matplotlib.pyplot as plt
    plt.show()



def find_action_potentials(voltages, fs):

    aps = []
    for v in voltages.T:
        positive = np.array((v > 0), dtype=int)
        diff = np.diff(positive)

        ## Range of the last AP
        imin = np.where(diff==+1)[0][-1]
        imax = np.where(diff==-1)[0][-1]


        icenter = np.argmax(v[imin:imax])

        aps.append(imin + icenter)

    aps = np.array(aps) / fs

    return aps

    # import matplotlib.pyplot as plt
    # plt.imshow(positive, aspect='auto', interpolation='nearest')
    # plt.show()
    # print(voltages)



def electrical_pulse_parameters(
        shape,
        **kwargs
):

    if shape == 'monophasic':

        polarity = kwargs['polarity']
        duration = kwargs['duration']

        amplitude = polarity / durations






    else:
        raise RuntimeError("Unknown pulse shape")
