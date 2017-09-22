# -*- coding: utf-8 -*-

"""Test geometry of the SG neuron.

"""
from __future__ import division, print_function, absolute_import
from __future__ import unicode_literals

import numpy as np
from numpy.testing import assert_equal

import spiral_ganglion as sg


def test_get_segment_path_positions():

    anf = sg.ANF_Axon(nodes=2)

    anf.set_geometry('straight', x0=0, y0=0, z0=0)

    # default section length (L): [10, 250, 1, 250, ...]

    expected = np.array([
        5 * 1e-6,
        10e-6 + 250e-6/2,
        10e-6 + 250e-6 + 1e-6/2,
        10e-6 + 250e-6 + 1e-6 + 250e-6/2
    ])

    actual = anf.get_segment_path_positions()

    assert_equal(actual, expected)
