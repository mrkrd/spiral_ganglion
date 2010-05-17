#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import sys
import numpy as np
import scipy.io

import spiral_ganglion as sg

def main(fin, fout):

    mat = scipy.io.loadmat(fin)
    fs = mat['fs'][0][0]
    stim = mat['stim']

    assert isinstance(fs, float)
    assert stim.ndim == 2

    trains = sg.run_ci_simulation(fs, stim, anf_num=100)
    trains = np.array(trains, dtype=object)

    scipy.io.savemat(fout, {'spike_trains': trains})

    # Plot results
    try:
        import thorns as th
        plot = th.plot_raster(trains) #, symboltype='circle')
        plot.xrange = (0, 1000*len(stim)/fs)
        plot.show()
    except ImportError:
        pass


if __name__ == "__main__":
    fin = sys.argv[1]
    fout = sys.argv[2]

    main(fin, fout)
