#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import
from __future__ import unicode_literals

__version__ = "1"

import numpy as np

import neuron
from neuron import h

import os
lib_dir = os.path.dirname(__file__)
neuron.load_mechanisms(lib_dir)

from . anf import (
    ANF_Axon,
    plot_geometry,
    plot_vectors,
    generate_anf_channels,
)

from . electrodes import (
    Electrode,
    calculate_medel_electrode_z_position,
)

from . ci import (
    run_ci_simulation,
    make_anf_electrode,
    find_threshold,
)

from . ear import run_holmberg2007_sg

from cochlea.holmberg2007 import get_nearest_cf as get_nearest_cf_holmberg2007



def run(duration, anfs):
    """Run a simulation of spiral ganglion neurons.

    This function takes care of proper acoustic and electric
    initialization.

    Parameters
    ----------
    duration : float
        Duration of the simulation in seconds.
    anfs : list of sg.ANF objects
        List of sg.ANF objects for initialization.

    """
    for anf in anfs:
        if anf.electrodes:
            anf.einit()

    neuron.init()

    for anf in anfs:
        if len(anf.vesicles) > 0:
            anf.ainit()

    neuron.run(duration * 1e3)      # s -> ms


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

        # Range of the last AP
        imin = np.where(diff == +1)[0][-1]
        imax = np.where(diff == -1)[0][-1]

        icenter = np.argmax(v[imin:imax])

        aps.append(imin + icenter)

    aps = np.array(aps) / fs

    return aps

    # import matplotlib.pyplot as plt
    # plt.imshow(positive, aspect='auto', interpolation='nearest')
    # plt.show()
    # print(voltages)
