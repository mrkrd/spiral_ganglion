#!/usr/bin/env python

from __future__ import division
__author__ = "Marek Rudnicki"

import numpy as np
from scipy.interpolate import interp1d

import neuron
from neuron import h

class ANF(object):
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


    def set_geometry(self, geo_type, **kwargs):
        """
        Sets position of each segment of the neuron.

        'straight'
        x0, y0, z0: coordinates of the initial segment
        x0, y0: are in the plane of cross section of scala tympani
        z0: coordinate along cochlea

        'bent'
        a, b: distances to the neuron along x and y axis respectively
        z: position along cochlea

        """
        if geo_type == 'straight':
            self._x = kwargs['x0'] + self._get_segment_positions()
            self._y = kwargs['y0'] * np.ones(len(self._get_segment_positions()))
            self._z = kwargs['z0'] * np.ones(len(self._get_segment_positions()))

        elif geo_type == 'bent':
            a = kwargs['a']
            b = kwargs['b']
            z = kwargs['z']
            r = b
            shift = np.abs(a-b)

            L = self._get_segment_positions()
            circumference = 2 * np.pi * r
            quarter = circumference / 4

            self._z = kwargs['z'] * np.ones(len(L))

            # Initial straight segment
            sel = L[ (L<shift) ]
            x = sel
            y = b * np.ones(len(sel))

            # Quarter circle segment
            sel = L[ (L>=shift) & (L<(quarter+shift)) ] - shift
            x = np.append(x, shift + r*np.sin(2 * np.pi * sel / circumference))
            y = np.append(y, r*np.cos(2 * np.pi * sel / circumference))

            # Remaining straight segment
            sel = L[ (L>=(quarter+shift)) ] - shift - quarter
            x = np.append(x, a*np.ones(len(sel)))
            y = np.append(y, -sel)

            # import matplotlib.pyplot as plt
            # plt.plot(x,y, 'o')
            # plt.show()

    def ainit(self):
        """
        Initializes acoustical stimulation.

        Note: must be called *after* neuron.init()

        """
        assert h.dt == 0.005
        for v in self.vesicles:
            self._con.event(float(v))



    def einit(self):
        """
        Initializes electrical stimulation.

        Note: must be called *before* neuron.init()

        """
        assert h.dt == 0.005

        self._stim_vectors = []
        if self.electrodes:
            # Make sure that each segment has defined position
            assert len(self._x) == len(self._y) == len(self._z) == \
                len(self._get_segment_positions())

            # Make sure all electrodes have signals of the same fs
            fss = np.array([el.fs for el in self.electrodes])
            assert np.all(fss == fss[0])
            el_dt = float(1000/fss[0])
            assert fss[0] is not None

            potentials = sum([el.calculate_potentials(self) for el in self.electrodes])

            for pot,seg in zip(potentials, self._get_segments()):
                vec = h.Vector(pot)
                vec.play(seg._ref_e_extracellular, el_dt)
                self._stim_vectors.append(vec)




class ANF_With_Soma(ANF):
    def __init__(self):
        print "ANF temperature:", h.celsius, "C"

        self.vesicles = []      # vesicle timings for acoustical stimulation
        self._x = None          # array of segment's x coordinate locations
        self._y = None          # array of segment's y coordinate locations
        self._z = None          # array of segment's z coordinate locations

        self.electrodes = []    # electrodes that stimulate the neuron
                                # (class Electrode)

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
        term.insert('extracellular')
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
            inode.insert('extracellular')
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
            node.insert('extracellular')
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
        soma.insert('extracellular')
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
            node.insert('extracellular')
            sections.append( ('p_node', node) )

            inode = h.Section()
            inode.nseg = 1
            inode.L = 300
            inode.Ra = 100
            inode.diam = 2.5
            inode.cm = 1e-3
            inode.insert('leak_manis')
            inode.g_leak_manis = 1e-5
            inode.insert('extracellular')
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






    def _calc_soma_diams(self, soma, max_diam, p_diam, c_diam):

        pos = np.linspace(0, 1, soma.nseg)

        flin = interp1d([0, 1], [p_diam, c_diam])

        diams = max_diam * np.sin(pos*np.pi)**2 + flin(pos)

        for i,seg in enumerate(soma):
            seg.diam = diams[i]


class ANF_Axon(ANF):
    def __init__(self):
        print "ANF temperature:", h.celsius, "C"

        self.vesicles = []      # vesicle timings for acoustical stimulation
        self._x = None          # array of segment's x coordinate locations
        self._y = None          # array of segment's y coordinate locations
        self._z = None          # array of segment's z coordinate locations

        self.electrodes = []    # electrodes that stimulate the neuron
                                # (class Electrode)

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
        term.insert('extracellular')
        sections.append( ('p_term', term) )

        ### Peripherial Axon Nodes and Internodes
        for i in range(15):
            inode = h.Section()
            inode.nseg = 1
            inode.L = 250
            inode.Ra = 100
            inode.diam = 1.5
            inode.cm = 1e-3
            inode.insert('leak_manis')
            inode.g_leak_manis = 1e-5
            inode.insert('extracellular')
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
            node.insert('extracellular')
            sections.append( ('p_node', node) )



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



if __name__ == "__main__":
    # import biggles
    import thorns.nrn as thn

    h.dt = 0.005
    h.celsius = 37

    anf = ANF_Axon()

    h.topology()
    # h.psection(sec=anf.sections['sec'][0])

    v = thn.record_voltages(anf.sections['sec'])

    anf.vesicles = [2, 5]

    neuron.init()
    anf.ainit()
    neuron.run(10)

    thn.plot_voltages(1/h.dt, v).show()

    print "Spikes:", np.array(anf.spikes)



    print
    print "========================"
    from electrodes import Electrode
    anf.set_geometry('straight', x0=250, y0=500, z0=0)
    anf.set_geometry('bent', a=750, b=500, z=0)
    anf.vesicles = []


    el = Electrode()
    el.z = 0
    print el.x, el.y, el.z

    stim = np.zeros(1000)
    stim[280:300] = -0.3
    stim[300:320] = 0.3

    el.fs = 200000
    el.stim = stim

    anf.electrodes = [el]

    anf.einit()
    neuron.init()
    neuron.run(len(stim) * h.dt)

    thn.plot_voltages(1/h.dt, v).show()


