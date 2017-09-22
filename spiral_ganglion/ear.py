# -*- coding: utf-8 -*-

# Copyright 2014 Marek Rudnicki
#
# This file is part of spiral_ganglion.
#
# spiral_ganglion is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# spiral_ganglion is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with spiral_ganglion.  If not, see <http://www.gnu.org/licenses/>.

"""Inner ear model.

"""
from __future__ import division, print_function, absolute_import
from __future__ import unicode_literals

import numpy as np
import pandas as pd

import cochlea
import spiral_ganglion as sg
import thorns as th

def run_holmberg2007_sg(
        sound,
        fs,
        anf_num,
        seed,
        cf=None,
        weight=None,
        channels='schwarz1987_pure',
):
    """Run inner ear model by [Holmberg2007]_ in quantal mode combined
    with spiral ganglion neurons.

    Parameters
    ----------
    weight : float
        Synaptic weight of the IHC synapse.
    channels : str or dict
        Channels used in the ANF model.

    """
    vesicle_trains = cochlea.run_holmberg2007_vesicles(
        sound=sound,
        fs=fs,
        anf_num=anf_num,
        seed=seed,
        cf=cf,
    )

    duration = th.get_duration(vesicle_trains)

    sg.set_celsius(37)
    sg.set_fs(100e3)

    anfs = []
    for _, train in vesicle_trains.iterrows():
        anf = sg.ANF_Axon(
            weight=weight,
            channels=channels,
            meta={'type': train['type'], 'cf': train['cf']},
        )

        anf.vesicles = train['vesicles']
        anfs.append(anf)


    sg.run(
        duration=duration,
        anfs=anfs
    )

    trains = []
    for anf in anfs:
        trains.append(anf.get_trains())


    trains = pd.concat(trains)

    return trains



def main():
    pass


if __name__ == "__main__":
    main()
