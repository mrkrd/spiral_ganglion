#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

def generate_biphaseic_pulse(
        fs,                     # Hz
        first_pulse_width,      # s
        second_pulse_width,     # s
        gap_width,              # s
        polarity='ac',
        pad_width=0):           # s
    """Generate a biphasic pulse used in CI simulations.

    The output pulses will be normalized for the charge of 1e-9C (1nC)

    """
    polarities = {'ac': [+1, -1],
                  'ca': [-1, +1]}

    assert polarity in polarities


    unit_charge = 1e-9          # = 1 nC

    widths = np.array([first_pulse_width,
                       second_pulse_width])
    amps = unit_charge/2 / widths

    pulses = [
        np.ones( np.round(width * fs) )
        for width in widths
    ]
    pulses *= amps
    pulses *= polarities[polarity]


    gap = np.zeros( np.round(gap_width * fs) )

    pad = np.zeros( np.round(pad_width * fs) )




    signal = np.concatenate((
        pad,
        pulses[0],
        gap,
        pulses[1],
        pad
    ))


    return signal




def main():

    s = generate_biphaseic_pulse(
        fs=200e3,
        first_pulse_width=40e-6,
        second_pulse_width=20e-6,
        gap_width=10e-6,
        polarity='ca',
        pad_width=100e-6)

    plt.plot(s)
    plt.show()



if __name__ == "__main__":
    import matplotlib.pyplot as plt

    main()
