#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import multiprocessing
import os
#import biggles

from anf import ANF_Axon
from electrodes import Electrode

import neuron
from neuron import h

# import thorns.nrn as thn


def _simulate_anf_at( (z, electrodes) ):

    print
    print os.getpid(), z

    h.dt = 0.005
    h.celsius = 37

    anf = ANF_Axon()
    anf.set_geometry('straight', x0=0, y0=500, z0=z)

    anf.electrodes = electrodes

    tmax = max([len(el.stim) for el in electrodes])
    tmax = 1000 * tmax / electrodes[0].fs

    # v = thn.record_voltages(anf.sections['sec'])

    anf.einit()
    neuron.init()
    neuron.run(tmax)

    # thn.plot_voltages(1000/h.dt, v).show()

    return np.asarray(anf.spikes)



def run_ci_simulation(fs, stim, anf_num=10, nproc=None):

    electrodes = []

    if isinstance(stim, dict):
        for i in stim:
            el = Electrode(i+1)
            el.x = 300
            el.fs = fs
            el.stim = stim[i]
            electrodes.append(el)

    elif isinstance(stim, np.ndarray):
        assert stim.ndim == 2
        assert stim.shape[1] == 12

        for i,s in enumerate(stim.T):
            el = Electrode(i+1)
            el.x = 300
            el.fs = fs
            el.stim = s
            electrodes.append(el)


    z_anf = np.linspace(0, 35000, anf_num)

    if nproc is None:
        nproc = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=nproc)

    space = [(z, el)
             for z in z_anf
             for el in [electrodes]]

    trains = pool.map(_simulate_anf_at, space)

    return trains



def _find_threshold(anf, electrode):
    stim = electrode.stim
    fs = electrode.fs
    anf.electrodes = [electrode]

    def run_sim(fs, stim, amp):
        electrode.fs = fs
        electrode.stim = stim * amp

        anf.einit()
        neuron.init()
        tmax = 1000 * len(stim) / fs
        neuron.run(tmax)

        return len(anf.spikes)

    lo = 0
    hi = 0.2

    # find initial range: lo/hi
    while run_sim(fs, stim, hi) == 0:
        lo = hi
        hi = hi * 2

    # binary search for amp
    while (hi-lo) > 0.01*(hi+lo)/2:
        amp = (hi+lo)/2
        # print amp
        cnt = run_sim(fs, stim, amp)
        if cnt == 0:
            lo = amp
        elif cnt > 0:
            hi = amp
        else:
            assert False, "Spike count should never be < 0"

    # restore original stim
    electrode.stim = stim
    return amp




def find_threshold(fs, stim):
    h.dt = 0.005
    h.celsius = 37

    electrode = Electrode()
    electrode.x = 300
    electrode.z = 0
    electrode.stim = stim
    electrode.fs = fs

    anf = ANF_Axon()
    # anf.set_geometry('bent', a=750, b=500, z=0)
    anf.set_geometry('straight', x0=0, y0=500, z0=0)

    return _find_threshold(anf, electrode)






def main():
    import thorns as th

    fs = 100000
    stim = np.zeros(1000)
    stim[300:320] = -0.5
    # stim[5000:5010] = -0.5

    find_threshold(fs, stim)
    exit()

    stim_dict = {6: stim}
    trains = run_ci_simulation(fs, stim_dict, anf_num=50)
    p = th.plot_raster(trains, symboltype='circle' )
    p.xrange = (0, 1000*len(stim)/fs)
    p.show()


if __name__ == "__main__":
    main()
