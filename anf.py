#!/usr/bin/env python

from __future__ import division
__author__ = "Marek Rudnicki"

import numpy as np
from scipy.interpolate import interp1d

import neuron
from neuron import h

class ANF_With_Soma(object):
    def __init__(self):
        print "ANF temperature:", h.celsius, "C"

        self.vesicles = []      # vesicle timings for acoustical stimulation
        self._x = None
        self._y = None
        self._z = None

        sections = []

        ### Peripherial Axon Terminal
        term = h.Section()
        term.nseg = 1
        term.L = 10
        term.Ra = 100
        term.diam = 1.5
        term.cm = 0.9
        term.insert('leak_manis')
        term.g_leak_manis = 1e-5
        term.insert('na_manis')
        term.gnabar_na_manis = 1
        term.insert('kht_manis')
        term.gkhtbar_kht_manis = 0.105
        term.insert('klt_manis')
        term.gkltbar_klt_manis = 0.027
        term.insert('ih_manis')
        term.ghbar_ih_manis = 0.016
        sections.append( ('p_term', term) )

        ### Peripherial Axon Nodes and Internodes
        for i in range(5):
            inode = h.Section()
            inode.nseg = 1
            inode.L = 250
            inode.Ra = 100
            inode.diam = 1.5
            inode.cm = 1e-3
            inode.insert('leak_manis')
            inode.g_leak_manis = 1e-5
            sections.append( ('p_inode', inode) )

            node = h.Section()
            node.nseg = 1
            node.L = 1
            node.Ra = 100
            node.diam = 1.5
            node.cm = 0.9
            node.insert('leak_manis')
            node.g_leak_manis = 1e-5
            node.insert('na_manis')
            node.gnabar_na_manis = 1
            node.insert('kht_manis')
            node.gkhtbar_kht_manis = 0.105
            node.insert('klt_manis')
            node.gkltbar_klt_manis = 0.027
            node.insert('ih_manis')
            node.ghbar_ih_manis = 0.016
            sections.append( ('p_node', node) )


        ### Soma
        soma = h.Section()
        soma.nseg = 9
        soma.L = 30
        soma.Ra = 100
        soma.diam = 30
        soma.cm = 1e-2
        self._calc_soma_diams(soma, 30, 1.5, 2.5)
        soma.insert('leak_manis')
        soma.g_leak_manis = 1e-5
        sections.append( ('soma', soma) )



        ### Central Axon Nodes and Internodes
        for i in range(10):
            node = h.Section()
            node.nseg = 1
            node.L = 1
            node.Ra = 100
            node.diam = 2.5
            node.cm = 0.9
            node.insert('leak_manis')
            node.g_leak_manis = 1e-5
            node.insert('na_manis')
            node.gnabar_na_manis = 1
            node.insert('kht_manis')
            node.gkhtbar_kht_manis = 0.105
            node.insert('klt_manis')
            node.gkltbar_klt_manis = 0.027
            node.insert('ih_manis')
            node.ghbar_ih_manis = 0.016
            sections.append( ('p_node', node) )

            inode = h.Section()
            inode.nseg = 1
            inode.L = 300
            inode.Ra = 100
            inode.diam = 2.5
            inode.cm = 1e-3
            inode.insert('leak_manis')
            inode.g_leak_manis = 1e-5
            sections.append( ('p_inode', inode) )



        # Remove the last inode
        if sections[-1][0] == 'p_inode':
            sections.pop()


        self.sections = np.rec.fromrecords(sections, names='typ,sec')
        for sec in self.sections['sec']:
            sec.v = -60

        # Connect sections
        for prev,next in zip(self.sections['sec'][:-1], self.sections['sec'][1:]):
            next.connect(prev)


        ### IHC Synapse
        q10 = np.exp(np.log(129.4 / 81.7) * (10. / (37 - 22)))
        self._syn = h.Exp2Syn(self.sections['sec'][0](0.5))
        self._syn.tau1 = 0.100906 / q10 ** ((h.celsius - 22) / 10)
        self._syn.tau2 = 0.592521 / q10 ** ((h.celsius - 22) / 10)
        self._con = h.NetCon(None, self._syn)
        self._con.weight[0] = 0.001864 * q10 ** ((h.celsius - 22) / 10)
        self._con.delay = 0


        ### Recording spikes from the last section
        last = self.sections['sec'][-1]
        self._probe = h.NetCon(last(0.5)._ref_v, None, 0, 0, 0, sec=last)
        self.spikes = h.Vector()
        self._probe.record(self.spikes)





    def _get_segments(self):
        """ Returns list of all segments along the neuron """
        segments = []
        for sec in self.sections['sec']:
            idx = [seg.x for seg in sec]
            segments.extend([sec(i) for i in idx])
        return segments

    def _get_segment_positions(self):
        """ Returns a list of positions of all segments along the neuron in um """
        pos = []                # vector of segment's positions
        start = 0               # position of the beginning of the
                                # current section
        for sec in self.sections['sec']:
            seg_x = np.array([seg.x for seg in sec])
            pos.extend(start + seg_x*sec.L)
            start = start + sec.L

        return np.array(pos)


    def set_geometry(self, x0, y0, z0):
        """
        Sets position of each segment of the neuron.

        Cartesian coordinates with respect to the base of scala tympani.

        x0, y0, z0: coordinates of the initial segment

        x0, y0: are in the plane of cross section of scala tympani
        z0: coordinate along cochlea

        """
        self._x = x0 + self._get_segment_positions()
        self._y = y0 * np.ones(len(self._get_segment_positions()))
        self._z = z0 * np.ones(len(self._get_segment_positions()))


    def init(self):
        """
        Must be called in acoustical stimulation between:
        neuron.init() and neuron.run(xx)

        """
        for v in self.vesicles:
            self._con.event(float(v))



    def _calc_soma_diams(self, soma, max_diam, p_diam, c_diam):

        pos = np.linspace(0, 1, soma.nseg)

        flin = interp1d([0, 1], [p_diam, c_diam])

        diams = max_diam * np.sin(pos*np.pi)**2 + flin(pos)

        for i,seg in enumerate(soma):
            seg.diam = diams[i]



if __name__ == "__main__":
    import biggles
    import thorns.nrn as thn

    h.dt = 0.005
    h.celsius = 37

    anf = ANF_With_Soma()

    h.topology()

    v = thn.record_voltages(anf.sections['sec'])

    anf.vesicles = [2, 5]

    neuron.init()
    anf.init()
    neuron.run(15)

    print anf.sections['sec'][0].v

    # thn.plot_voltages(1/h.dt, v).show()

    print np.array(anf.spikes)

    print anf._get_segment_positions()

    anf.set_geometry(1, 2, 3)
    print anf._x
    print anf._y
    print anf._z
