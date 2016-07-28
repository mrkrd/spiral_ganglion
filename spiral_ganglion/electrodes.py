#!/usr/bin/env python

from __future__ import division, print_function, absolute_import
from __future__ import unicode_literals

import numpy as np


class Electrode(object):
    """Point electrode located at (x, y, z) coordinates.

    Attributes
    ----------
    x, y, z : float
        Electrode position (m).
    stim : array_like
        Electric current (A).
    fs : scalar
        Sampling frequency of the electric current vector (Hz).

    """
    def __init__(self):
        self.x = 300e-6         # [m]
        self.y = 0              # [m]
        self.z = 0              # [m]

        self.stim = None        # [A]
        self.fs = None          # [Hz]

    def calculate_potentials(self, anf):
        """
        Calculate electrical potential for each segment of the ANF.

        """
        assert self.stim is not None, "Stimulus must be set."
        assert self.fs is not None, "Sampling frequency must be set."

        anf_pos = anf.get_positions()

        # Calculate decay along cochlea
        gain_dB = -1e3          # [dB/m]: -1dB/mm
        z_anf = anf_pos['z'][0]
        assert np.all(anf_pos['z'] == z_anf)

        exponent = gain_dB * np.abs(self.z - z_anf) / 20
        z_decay_factor = 10**exponent
        stim = z_decay_factor * self.stim

        # Calculate homogenious medium (1/r)
        resistivity = 3         # resistivity = 300 ohm*cm = 3 ohm*m
        r = np.sqrt((self.x - anf_pos['x'])**2 + (self.y - anf_pos['y'])**2)
        node_factors = resistivity / (4 * np.pi * r)

        potentials = np.outer(node_factors, stim)

        return potentials


def calculate_medel_electrode_z_position(electrode_number):
    z = np.linspace(4.6e-3, 31e-3, 12)[electrode_number-1]

    return z
