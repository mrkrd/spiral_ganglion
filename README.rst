spiral_ganglion
===============

*Spiral ganglion neuron model for synaptic and electrical stimulation*

:Name: spiral_ganglion
:Author: Marek Rudnicki
:Email: marek.rudnicki@tum.de
:License: GNU General Public License v3 or later (GPLv3+)


Dependencies
------------

- numpy_
- neuron_
- thorns_ (only to run example scripts)

.. _numpy: http://www.numpy.org/
.. _neuron: http://www.neuron.yale.edu/neuron/
.. _thorns: https://github.com/mrkrd/thorns


Installation
------------

::

   # Install in the developer mode
   python setup.py develop --user

   # Compile mod files (channels)
   cd spiral_ganglion
   nrnivmodl
