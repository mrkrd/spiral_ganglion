spiral_ganglion
===============

Spiral ganglion neuron model for synaptic and electrical stimulation.



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

- cochlea_
- thorns_ (only to run example scripts)

.. _NEURON: http://www.neuron.yale.edu/neuron/
.. _thorns: https://github.com/mrkrd/thorns
.. _cochlea: https://github.com/mrkrd/cochlea


Installation
------------

::

   # Install in the developer mode
   python setup.py develop --user

   # Compile mod files (channels)
   cd spiral_ganglion
   nrnivmodl



Citing
------

Rudnicki, M. (2018) *Computer models of acoustical and electrical stimulation
of neurons in the auditory system*
PhD thesis. Technische Universität München.
https://mediatum.ub.tum.de/1445042



License
-------

The project is licensed under the GNU General Public License v3 or
later (GPLv3+).
