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
        self.x = 300            # um
        self.y = 0              # um
        self.z = None           # um

        self.stim = None        # mA
        self.fs = None          # Hz

        if el_num is not None:
            self.z = np.linspace(4600, 31000, 12)[el_num-1]

    def calculate_potentials(self, anf):
        """
        Calculate electrical potential for each segment of the ANF.

        """
        assert self.stim is not None, "Stimulus must be set."
        assert self.fs is not None, "Sampling frequency must be set."

        ### Calculate decay along cochlea
        gain_dB = -3            # -3dB / mm
        assert np.all(anf._z == anf._z[0])
        exponent = gain_dB * np.abs(self.z - anf._z[0]) / 1000 / 20
        z_decay = 10**exponent
        stim = z_decay * self.stim

        ### Calculate homogenious medium (1/r)
        r = np.sqrt((self.x - anf._x)**2 + (self.y - anf._y)**2)
        node_factors = 300 * 1e4 / (4 * np.pi * r)

        # import matplotlib.pyplot as plt
        # plt.plot(node_factors)
        # plt.show()

        return np.outer(node_factors, stim)


def main():
    pass

if __name__ == "__main__":
    main()
