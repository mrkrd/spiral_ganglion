#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

def generate_biphaseic_pulse(
        fs,                     # Hz
        widths,                 # [s, s]
        gap_width,              # s
        polarity='ac',
        pad_width=[0]):         # s
    """Generate a biphasic pulse used in CI simulations.

    The output pulses will be normalized for the charge of 1e-9 C (1nC).

    The output is in [mA].

    """
    polarities = {'ac': [+1, -1],
                  'ca': [-1, +1]}

    widths = np.array(widths)
    pad_width = np.array(pad_width)

    assert widths.size == 2
    assert pad_width.size in (1,2)
    assert polarity in polarities


    unit_charge = 1e-9          # = 1 nC


    amps = unit_charge/2 / widths

    pulses = []
    for i,width in enumerate(widths):
        pulse = np.ones( np.round(width * fs) )
        pulse *= amps[i]
        pulse *= polarities[polarity][i]

        pulses.append( pulse )

    gap = np.zeros( np.round(gap_width * fs) )

    if pad_width.size == 1:
        pad_width = np.repeat(pad_width, 2)
    pad = [
        np.zeros( np.round(w * fs) )
        for w in pad_width
    ]




    signal = np.concatenate((
        pad[0],
        pulses[0],
        gap,
        pulses[1],
        pad[1]
    ))


    signal *= 1e3               # convert: A -> mA

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
