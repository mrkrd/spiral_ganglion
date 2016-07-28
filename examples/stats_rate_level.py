#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2009-2014 Marek Rudnicki

# This file is part of spiral_ganglion.

# spiral_ganglion is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# spiral_ganglion is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with spiral_ganglion.  If not, see <http://www.gnu.org/licenses/>.


"""This demo shows how to calculate rate-level characteristic of a
model.

"""
from __future__ import division, absolute_import, print_function

import numpy as np
import matplotlib.pyplot as plt

import spiral_ganglion as sg
from cochlea.stats import calc_rate_level
import cochlea


def main():

    rates = calc_rate_level(
        model=sg.run_holmberg2007_sg,
        cf=cochlea.get_nearest_cf_holmberg2007(1000),
        model_pars={'fs': 48e3, 'anf_num': (10,10,10)},
    )

    print(rates)

    rates.plot()

    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
