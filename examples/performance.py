#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import numpy as np

import spiral_ganglion as sg


# cpu: Intel(R) Core(TM) i7-3770 CPU @ 3.40GHz
# command: time python performance.py
# fs=500e3, anf_num=100, tmax=0.1
# time: 391.72s user 0.04s system 99% cpu 6:31.92 total

# expected: fs=500e3, anf_num=10e3, tmax=0.5
# expected time: 100*5*time ~ 196000s ~ 55h

# 1000 neurons => 0.5h (tmax=0.5)

# possible optimizations:
# - lower fs
# - lower compartment count


def main():
    fs = 500e3
    anf_num = 100
    electrode_num = 12
    tmax = 0.1

    sg.set_fs(fs)
    sg.set_celsius(37)

    stim = np.zeros(tmax*fs)


    electrodes = []
    for i in range(electrode_num):
        electrode = sg.Electrode()
        electrode.fs = fs
        electrode.stim = stim



    anfs = []
    for i in range(anf_num):
        anf = sg.ANF_Axon()
        anf.electrodes = electrodes

        anfs.append(anf)


    sg.run(
        tmax,
        anfs
    )



if __name__ == "__main__":
    main()
