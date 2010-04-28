#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

class Electrode(object):
    def __init__(self, el_num=None):
        self.x = 0              # um
        self.y = 0              # um
        self.z = None           # um

        self.stim = None        # mA
        self.fs = None          # Hz

        if el_num is not None:
            self.z = np.linspace(4600, 31000, 12)[el_num]

    def calculate_potentials(self, anf):
        """
        Calculate electrical potential for each segment of the ANF.

        """
        ### Calculate decay along cochlea
        gain_dB = -3            # -3dB / mm
        assert np.all(anf._z == anf._z[0])
        exponent = gain_dB * np.abs(self.z - anf._z[0]) / 1000 / 20
        z_decay = 10**exponent
        stim = z_decay * self.stim

        ### Calculate homogenious medium (1/r)
        r = np.sqrt(anf._x**2 + anf._y**2)
        node_factors = 1e4 / (4 * np.pi * r * 300)

        return np.outer(node_factors, stim)


def main():
    pass

if __name__ == "__main__":
    main()
