#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

def generate_pulse(
        fs,                     # Hz
        widths,                 # [us, ..., us]
        gap_width,              # us
        polarities,             # e.g. 'ac' where len(polarities) == len(widths)
        pad_widths=[0]):         # ms
    """Generate a biphasic pulse used in CI simulations.

    Each of the pulses is normalized for the charge of 1e-9 C (1nC).

    The output is in [mA].

    """
    polarity_map = {'a':+1, 'c':-1}

    widths = np.array(widths) * 1e-6       # us -> s
    gap_width = np.array(gap_width) * 1e-6 # us -> s
    pad_widths = np.array(pad_widths) * 1e-3 # ms -> s

    assert len(polarities) == len(widths)
    assert pad_widths.size in (1,2)


    unit_charge = 1e-9          # = 1 nC


    pulses = []
    for width,polarity in zip(widths, polarities):
        pulse = np.ones( np.round(width * fs) )
        pulse *= unit_charge / width # amplitude
        pulse *= polarity_map[polarity]
        pulses.append( pulse )

        gap = np.zeros( np.round(gap_width * fs) )
        pulses.append( gap )

    pulses.pop()                # removes the last gap


    if pad_widths.size == 1:
        pad_widths = np.repeat(pad_widths, 2)
    pulses.insert(0, np.zeros( np.round(pad_widths[0] * fs) ))
    pulses.append(np.zeros( np.round(pad_widths[0] * fs) ))


    signal = np.concatenate(pulses)

    signal *= 1e3               # convert: A -> mA

    return signal




def main():

    s = generate_pulse(
        fs=200e3,
        widths=[40,20],
        gap_width=20,
        polarities='ca',
        pad_widths=1
    )

    plt.plot(s)
    plt.show()



if __name__ == "__main__":
    import matplotlib.pyplot as plt

    main()
