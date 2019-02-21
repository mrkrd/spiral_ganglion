#!/usr/bin/env python

from __future__ import division, print_function, absolute_import

import numpy as np
import multiprocessing

import pandas as pd

from neuron import h

import spiral_ganglion as sg

import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def run_ci_simulation(
        stim,
        fs,
        anf_num=10,
        map_backend='serial',
        return_voltages=False
):
    """Run CI simulation where 12 electrodes are located along the
    cochlea and `anf_num' of auditory nerve fibers is spread
    uniformly.  Returns a list of spike trains or recorded membrane
    potentials (return_voltages=True).

    Parameters
    ----------
    stim : dict or array_like
        Electrical stmulation matrix.
    fs : float
        Sampling frequency of the simulus matrix.
    anf_num : int, optional
        Number of ANFs that are equally spread along the cochlea.
    map_backend : {'serial', 'multiprocessing'}
        Use either serial or parallel backend.
    return_voltages : bool
        If True voltage membranes are returned instead of spike trains.

    Returns
    -------
    spike_trains
        ANF spike trains.

    """
    electrodes = []

    if isinstance(stim, dict):
        for i in stim:
            el = sg.Electrode()
            el.x = 300e-6
            el.z = sg.calculate_medel_electrode_z_position(i+1)
            el.fs = fs
            el.stim = stim[i]
            electrodes.append(el)

    elif isinstance(stim, np.ndarray):
        assert stim.ndim == 2
        assert stim.shape[1] == 12

        for i, s in enumerate(stim.T):
            el = sg.Electrode()
            el.x = 300e-6
            el.z = sg.calculate_medel_electrode_z_position(i+1)
            el.fs = fs
            el.stim = s
            electrodes.append(el)

    else:
        raise NotImplementedError("Unimplemented stimulus format.")

    anf_zs = np.linspace(0, 35e-3, anf_num)

    if return_voltages:
        raise NotImplementedError()

    space = [
        (z, electrodes, return_voltages)
        for z in anf_zs
    ]

    if map_backend == 'serial':
        trains = map(_simulate_anf_at, space)
    elif map_backend == 'multiprocessing':
        pool = multiprocessing.Pool()
        trains = pool.map(_simulate_anf_at, space)
    else:
        raise NotImplementedError("Map backend not implementation: {}".format(map_backend))

    trains = pd.concat(trains)

    return trains


def _simulate_anf_at((z, electrodes, return_voltages)):

    logger.info("Calculating ANF at z = {} [m]".format(z))

    sg.set_fs(500e3)
    sg.set_celsius(37)

    anf = sg.ANF_Axon(record_voltages=return_voltages)
    anf.set_geometry('straight', x0=0, y0=500e-6, z0=z)

    anf.electrodes = electrodes

    tmax = max([len(el.stim) for el in electrodes])
    duration = tmax / electrodes[0].fs

    sg.run(
        duration=duration,
        anfs=[anf]
    )

    if return_voltages:
        return anf.get_voltages()
    else:
        return anf.get_trains()


def find_threshold(
        anf,
        electrode,
        stimulus,
        fs,
        pre_stimulus=None,
        pad=False,
        error=1e-6
):
    """Find firing threshold of an spiral ganglion neuron stimulated by an
    electrode.

    Parameters
    ----------
    anf : ANF_Axon
        Spiral ganglion neuron.
    electrode : Electrode
        Stimulating electrode.
    stimulus : array_like
        Electrical stimulus.
    fs : float
        Sampling frequency of the stimulus.
    pre_stimulus : array_like, optional
        Conditioning stimulus.  The amplitude of `pre_stimulus` is
        constant during simulations.
    pad : bool, optional
        If True, then add extra padding to the stimulus.
    error : float, optional
        Error of the threshold search at which break the search loop.

    Returns
    -------
    float
        Firing threshold.

    """

    # h.dt = 0.002                # [ms]

    electrode.fs = fs
    anf.electrodes = [electrode]

    if pre_stimulus is None:
        pre_stimulus = np.array([])

    lo = 0
    hi = 1e-12

    # find initial range: lo/hi
    while True:
        above_threshold = _is_above_threshold(
            anf=anf,
            electrode=electrode,
            amplitude=hi,
            stimulus=stimulus,
            fs=fs,
            pre_stimulus=pre_stimulus,
            pad=pad
        )
        if above_threshold:
            break

        logger.debug(" {:>20}  {:<20}".format(lo, hi))

        lo = hi
        hi = hi * 1.5

    logger.debug("Maximum value found: {hi}".format(hi=hi))

    # binary search for amp
    while (hi-lo) > error*(hi+lo)/2:
        amp = (hi+lo)/2

        above_threshold = _is_above_threshold(
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

        logger.debug(" {:>20}  {:<20}, {}, {}".format(lo, hi, amp, above_threshold))

        if above_threshold:
            hi = amp
        else:
            lo = amp

    return amp


def _is_above_threshold(
        anf,
        electrode,
        amplitude,
        stimulus,
        fs,
        pre_stimulus,
        pad
):
    if pad:
        pre_pad = np.zeros(10e-3 * fs)
        post_pad = np.zeros(5e-3 * fs)
    else:
        pre_pad = np.array([])
        post_pad = np.array([])

    electrode.stim = np.concatenate(
        (pre_pad, pre_stimulus, stimulus*amplitude, post_pad)
    )

    sg.run(
        duration=(len(electrode.stim) / electrode.fs),
        anfs=[anf]
    )

    anf_trains = anf.get_trains()

    if (anf_trains is None) or (anf_trains['spikes'][0].size > 0):
        above_threshold = True
    else:
        above_threshold = False

    return above_threshold


def make_anf_electrode():
    h.celsius = 37

    electrode = sg.Electrode()
    electrode.x = 300e-6

    anf = sg.ANF_Axon()
    # anf.set_geometry('bent', a=750, b=500, z=0)
    anf.set_geometry('straight', x0=0, y0=500e-6, z0=0)

    return anf, electrode


def main():
    import matplotlib.pyplot as plt
    from anf import plot_geometry
    anf, electrode = make_anf_electrode()

    plot_geometry([anf, electrode])
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
    print(th)

    exit()

    stim_dict = {6: stim}
    trains = run_ci_simulation(fs, stim_dict, anf_num=50)
    p = th.plot_raster(trains, symboltype='circle')
    p.xrange = (0, 1000*len(stim)/fs)
    p.show()


if __name__ == "__main__":
    main()
