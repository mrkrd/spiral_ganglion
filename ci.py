#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import biggles

from anf import ANF_Axon
from electrodes import Electrode

import neuron
from neuron import h

import thorns.nrn as thn

def find_threshold(fs, stim, anf=None):
    h.dt = 0.005
    h.celsius = 37

    if anf == None:
        electrode = Electrode()
        electrode.x = 0
        electrode.z = 0

        anf = ANF_Axon()
        # anf.set_geometry('bent', a=750, b=500, z=0)
        anf.set_geometry('straight', x0=0, y0=500, z0=0)
        anf.electrodes = [electrode]
    else:
        assert len(anf.electrodes) == 1
        electrode = anf.electrodes[0]


    # v = thn.record_voltages(anf.sections['sec'])


    def run_sim(fs, stim, amp):
        electrode.fs = fs
        electrode.stim = stim * amp

        anf.einit()
        neuron.init()
        tmax = len(stim) * 1000 / fs
        neuron.run(tmax)

        return len(anf.spikes)


    lo = 0
    hi = 0.5

    # find initial range: lo/hi
    while run_sim(fs, stim, hi) == 0:
        lo = hi
        hi = hi * 2


    while (hi-lo) > 0.01*(hi+lo)/2:
        amp = (hi+lo)/2
        print amp
        cnt = run_sim(fs, stim, amp)
        if cnt == 0:
            lo = amp
        elif cnt > 0:
            hi = amp
        else:
            assert False, "Spike count should never be < 0"


    # thn.plot_voltages(1000/h.dt, v).show()
    return amp



def main():
    fs = 200000
    stim = np.zeros(1000)
    stim[300:320] = -1

    find_threshold(fs, stim)


if __name__ == "__main__":
    main()
