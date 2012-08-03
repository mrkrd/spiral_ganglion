#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

class Electrode(object):
    """
    Represents electrode located at (self.x, self.y, self.z) coordinates.

    """
    def __init__(self, el_num=None):
        """
        el_num: electrode number (1-12)
        """
        self.x = 300e-6         # [m]
        self.y = 0              # [m]

        self.stim = None        # [A]
        self.fs = None          # [Hz]

        if el_num is not None:
            self.z = np.linspace(4.6e-3, 31e-3, 12)[el_num-1]
        else:
            self.z = 0          # [m]


    def calculate_potentials(self, anf):
        """
        Calculate electrical potential for each segment of the ANF.

        """
        assert self.stim is not None, "Stimulus must be set."
        assert self.fs is not None, "Sampling frequency must be set."


        ### Calculate decay along cochlea
        gain_dB = -1e3          # [dB/m]: -1dB/mm
        z = anf._z[0]
        assert np.all(anf._z == z)

        exponent = gain_dB * np.abs(self.z - z) / 20
        z_decay_factor = 10**exponent
        stim = z_decay_factor * self.stim



        ### Calculate homogenious medium (1/r)
        # resistivity = 300 ohm*cm = 3 ohm*m
        r = np.sqrt((self.x - anf._x)**2 + (self.y - anf._y)**2)
        node_factors = 3 / (4 * np.pi * r)

        # import matplotlib.pyplot as plt
        # plt.plot(node_factors)
        # plt.show()

        potentials = np.outer(node_factors, stim)

        return potentials


def main():
    pass

if __name__ == "__main__":
    main()
