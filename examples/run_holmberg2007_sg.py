# -*- coding: utf-8 -*-

"""Demo showing how to run an inner ear model by Holmberg (2007)
integrated with a spiral ganglion models.

"""

from __future__ import division, print_function, absolute_import
from __future__ import unicode_literals

import matplotlib.pyplot as plt

import spiral_ganglion as sg
import cochlea
import thorns as th
import thorns.waves as wv


def main():

    fs = 48e3
    cf = cochlea.get_nearest_cf_holmberg2007(1e3)


    ### Make sound
    sound = wv.ramped_tone(
        fs=fs,
        freq=cf,
        duration=150e-3,
        pad=10e-3,
        dbspl=70,
    )



    ### Run model
    anf_trains = sg.run_holmberg2007_sg(
        sound,
        fs,
        cf=cf,
        anf_num=(10,0,0),
        seed=0,
    )


    print(th.firing_rate(anf_trains))


    ### Plot auditory nerve response
    fig, ax = plt.subplots(2,1)
    th.plot_signal(
        signal=sound,
        fs=fs,
        ax=ax[0]
    )

    th.plot_raster(
        anf_trains,
        ax=ax[1]
    )

    plt.show()





if __name__ == "__main__":
    main()
