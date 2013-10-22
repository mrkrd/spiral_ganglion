#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

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
