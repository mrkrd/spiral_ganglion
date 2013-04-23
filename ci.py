#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
import multiprocessing
import os
import logging

from spiral_ganglion.anf import ANF_Axon
from spiral_ganglion.electrodes import Electrode

import neuron
from neuron import h



def _simulate_anf_at( (z, electrodes, return_voltages) ):

    print
    print os.getpid(), z

    if h.dt > 0.002:
        h.dt = 0.002
    h.celsius = 37

    anf = ANF_Axon(record_voltages=return_voltages)
    anf.set_geometry('straight', x0=0, y0=500e-6, z0=z)

    anf.electrodes = electrodes

    tmax = max([len(el.stim) for el in electrodes])
    tmax = 1000 * tmax / electrodes[0].fs

    anf.einit()
    neuron.init()
    neuron.run(tmax)


    if return_voltages:
        return anf.get_voltages()
    else:
        return np.asarray(anf.spikes)



def run_ci_simulation(fs, stim, anf_num=10, nproc=None, return_voltages=False):
    """
    Run CI simulation where 12 electrodes are located along the
    cochlea and `anf_num' of auditory nerve fibers is spread
    uniformly.  Returns a list of spike trains or recorded membrane
    potentials (return_voltages=True).

    Simulation is run in multiprocessing mode, i.e. each ANF is
    assigned to a different processor.  The more processor the faster
    the simulation!

    fs: sampling frequency of the simulus
    stim: electrical stmulation (dict/np.array)
    anf_num: number of ANF to compute (are equally spread over cochlea)
    nproc: number of processes to spawn
    return_voltages: if True voltage traces are returned instead of spike trains

    """
    electrodes = []

    if isinstance(stim, dict):
        for i in stim:
            el = Electrode(i+1)
            el.x = 300e-6
            el.fs = fs
            el.stim = stim[i]
            electrodes.append(el)

    elif isinstance(stim, np.ndarray):
        assert stim.ndim == 2
        assert stim.shape[1] == 12

        for i,s in enumerate(stim.T):
            el = Electrode(i+1)
            el.x = 300e-6
            el.fs = fs
            el.stim = s
            electrodes.append(el)


    z_anf = np.linspace(0, 35e-3, anf_num)

    if nproc is None:
        nproc = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=nproc)

    space = [(z, el, v)
             for z in z_anf
             for el in [electrodes]
             for v in [return_voltages]]

    trains = pool.map(_simulate_anf_at, space)

    return trains



def _run_single_electrode(
        anf,
        electrode,
        amplitude,
        stimulus,
        fs,
        pre_stimulus,
        pad):

    if pad:
        pre_pad = np.zeros(10e-3 * fs)
        post_pad = np.zeros(5e-3 * fs)
    else:
        pre_pad = np.array([])
        post_pad = np.array([])

    electrode.stim = np.concatenate(
        (pre_pad, pre_stimulus, stimulus*amplitude, post_pad)
    )

    anf.einit()
    neuron.init()
    tmax = 1000 * len(electrode.stim) / electrode.fs
    neuron.run(tmax)

    spikes = anf.get_spikes()['spikes'][0]

    return spikes



def find_threshold(
        anf,
        electrode,
        stimulus,
        fs,
        pre_stimulus=None,
        pad=False,
        error=1e-6):


    # h.dt = 0.002                # [ms]

    electrode.fs = fs
    anf.electrodes = [electrode]

    if pre_stimulus is None:
        pre_stimulus = np.array([])


    lo = 0
    hi = 1e-12

    # find initial range: lo/hi
    while True:
        spikes = _run_single_electrode(
            anf=anf,
            electrode=electrode,
            amplitude=hi,
            stimulus=stimulus,
            fs=fs,
            pre_stimulus=pre_stimulus,
            pad=pad
        )
        if spikes.size > 0:
            break

        logging.debug(" {:>20}  {:<20}".format(lo, hi))

        lo = hi
        hi = hi * 2

    logging.debug("Maximum value found: {hi}".format(hi=hi))

    # binary search for amp
    while (hi-lo) > error*(hi+lo)/2:
        amp = (hi+lo)/2

        spikes = _run_single_electrode(
            anf=anf,
            electrode=electrode,
            amplitude=amp,
            stimulus=stimulus,
            fs=fs,
            pre_stimulus=pre_stimulus,
            pad=pad
        )

        # print(spikes)
        # import spiral_ganglion as sg
        # sg.anf._plot_voltages(anf.get_voltages()[:,-5:-1])
        # import matplotlib.pyplot as plt
        # plt.show()

        logging.debug(" {:>20}  {:<20}".format(lo, hi))

        if spikes.size > 0:
            hi = amp
        else:
            lo = amp

    return amp



def make_anf_electrode():
    h.celsius = 37

    electrode = Electrode()
    electrode.x = 300e-6

    anf = ANF_Axon()
    # anf.set_geometry('bent', a=750, b=500, z=0)
    anf.set_geometry('straight', x0=0, y0=500e-6, z0=0)

    return anf, electrode




def main():
    import matplotlib.pyplot as plt
    from anf import plot_geometry
    anf, electrode = make_anf_electrode()

    plot_geometry( [anf,electrode] )
    plt.show()

    from signals import generate_pulse

    fs = 200e3
    stim = generate_pulse(
        fs=fs,
        widths=[40e-6],
        gap_width=0,
        polarities='c',
        pad_widths=[10e-3, 5e-3]
    )


    th = find_threshold(anf, electrode, stimulus=stim, fs=fs)
    print th

    exit()

    stim_dict = {6: stim}
    trains = run_ci_simulation(fs, stim_dict, anf_num=50)
    p = th.plot_raster(trains, symboltype='circle' )
    p.xrange = (0, 1000*len(stim)/fs)
    p.show()


if __name__ == "__main__":
    main()
