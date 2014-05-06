#!/usr/bin/env python

from __future__ import division, absolute_import, print_function

__author__ = "Marek Rudnicki"

import numpy as np

import matplotlib.pyplot as plt

import mrlib.thorns as th
import spiral_ganglion as sg


def main():
    ### Generate the stimulus
    fs = 10e3                   # [Hz]
    amp = 400e-6                # [A]

    s = np.zeros(30e-3*fs)      # 30 ms
    s[100:105] = -amp           # [A]
    s[105:110] = +amp           # [A]

    stim = {
        3: s,
        8: np.roll(s, 100)
    }


    ### Run CI simulation
    trains = sg.run_ci_simulation(
        stim=stim,
        fs=fs,
        anf_num=10,
        # map_backend='multiprocessing'
    )


    ### Plot results
    fig, ax = plt.subplots(2, 1, sharex=True)

    th.plot_signal(stim[3], fs=fs, ax=ax[0])
    th.plot_signal(stim[8], fs=fs, ax=ax[0])
    ax[0].set_ylabel("Amplitude [A]")

    th.plot_raster(trains, ax=ax[1])

    plt.show()


if __name__ == "__main__":
    main()
