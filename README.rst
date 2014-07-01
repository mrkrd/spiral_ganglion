spiral_ganglion
===============

*Spiral ganglion neuron model for synaptic and electrical stimulation*

:Name: spiral_ganglion
:Author: Marek Rudnicki
:Email: marek.rudnicki@tum.de
:License: GNU General Public License v3 or later (GPLv3+)


Requirements
------------

- Python (2.7)
- Numpy
- Pandas
- Matplotlib
- NEURON_

  - Python module, i.e. you should be able to type ``import neuron``
    in a Python shell
  - ``nrnivmodl`` must be available as command.  Update your ``PATH``
    environment variable, if necessary.

- thorns_ (only to run example scripts)

.. _NEURON: http://www.neuron.yale.edu/neuron/
.. _thorns: https://github.com/mrkrd/thorns


Installation
------------

::

   # Install in the developer mode
   python setup.py develop --user

   # Compile mod files (channels)
   cd spiral_ganglion
   nrnivmodl
