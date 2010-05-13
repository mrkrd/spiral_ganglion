#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

import thorns as th
import spiral_ganglion as sg


def main():
    ### Set the stimulus
    fs = 10e3
    stim = np.zeros(300)         # 30 ms
    stim[100:105] = -0.4         # mA
    stim[105:110] = +0.4         # mA

    stim_dict = {3: stim,
                 8: np.roll(stim, 100)}

    trains = sg.run_ci_simulation(fs, stim_dict, anf_num=50)

    ### Plot results
    plot = th.plot_raster(trains, symboltype='circle' )
    plot.xrange = (0, 1000*len(stim)/fs)
    plot.show()



if __name__ == "__main__":
    main()
