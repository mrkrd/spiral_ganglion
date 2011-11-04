#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np
from scipy.interpolate import interp1d

import neuron
from neuron import h


class ANF(object):
    def _Tf(self, q10, temp_ref=22):
        return q10 ** ((h.celsius - temp_ref)/10)

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

        Note: units in use micrometers

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

            self._x = x
            self._y = y
            self._z = kwargs['z'] * np.ones(len(L))

            # import matplotlib.pyplot as plt
            # plt.plot(x,y, 'o')
            # plt.show()


    def get_voltages(self):
        """
        Returns time courses of membrane potentials for each section
        of the model.  Potentials are recorded in the from the middle
        of each section: sec(0.5).

        Neuron must be initialized with `record_voltages=True'.

        """
        voltages = [np.asarray(vec) for vec in self._voltages]
        voltages = np.array(voltages).T
        return voltages


    def get_spikes(self):
        """Return an array of spike timings recorded from the last
        section.

        """
        assert h.t != 0, "Time is 0 (did you run the simulation already?)"

        train = np.array( [(np.array(self._spikes), h.t)],
                          dtype=[('spikes', np.ndarray),
                                 ('duration', float)] )

        return train


    def ainit(self):
        """
        Initializes acoustical stimulation.

        Note: must be called *after* neuron.init()

        """
        assert h.dt <= 0.01

        for sec in self.sections['sec']:
            sec.v = -60

        for v in self.vesicles:
            self._con.event(float(v))



    def einit(self, do_assert=True):
        """
        Initializes electrical stimulation.

        Note: must be called *before* neuron.init()

        """
        if do_assert:
            assert h.dt <= 0.002

        for sec in self.sections['sec']:
            sec.v = -60

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





class ANF_Axon(ANF):
    """
    Model of Auditory Nerve Fiber's peripherial axon.  Can be used
    for acoustical and electrical stimulation.

    Initially has be described by Paul Wilhem Bade in his Master's
    Thesis.  This implementation includes modified Na+ channels
    (Rothman 1993) that correct refractory period in acoustical
    stimulation.  Na+ maximum conductance has also been changed (1
    S/cm2 -> 0.3 S/cm2).

    anf = ANF_Axon()

    Acoustical stimulation
    ======================
    anf.vesicles = [12, 24, 55]
    neuron.init()
    anf.ainit()
    neuron.run(60)


    Electrical stimulation
    ======================
    anf.electrodes = [electrode1, electrode2, electrode3]
    anf.einit()
    neuron.init()
    neuron.run(60)

    """
    def __init__(self,
                 nodes=20,
                 na_type='rothman93',
                 gna=0.324,
                 record_voltages=False):
        """nodes: number of nodes in the model.  Total number of
               compartments if 2*nodes.

        na_type: 'rothman93'/'orig' specifies which Na+ channel should
                 be used for the simulation

        record_voltages: when True membrane potentials are recorded
                         internally and can be returned with
                         get_voltages()

        """
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
        if na_type == 'rothman93':
            term.insert('na_rothman93')
            term.gnabar_na_rothman93 = gna
        elif na_type == 'orig':
            term.insert('na_manis')
            term.gnabar_na_manis = gna
        term.insert('kht_manis')
        term.gkhtbar_kht_manis = 0.105
        term.insert('klt_manis')
        term.gkltbar_klt_manis = 0.027
        term.insert('ih_manis')
        term.ghbar_ih_manis = 0.016
        term.insert('extracellular')
        sections.append( ('p_term', term) )

        ### Peripherial Axon Nodes and Internodes
        for i in range(nodes):
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
            if na_type == 'rothman93':
                node.insert('na_rothman93')
                node.gnabar_na_rothman93 = gna
            elif na_type == 'orig':
                node.insert('na_manis')
                node.gnabar_na_manis = gna
            node.insert('kht_manis')
            node.gkhtbar_kht_manis = 0.105
            node.insert('klt_manis')
            node.gkltbar_klt_manis = 0.027
            node.insert('ih_manis')
            node.ghbar_ih_manis = 0.016
            node.insert('extracellular')
            sections.append( ('p_node', node) )



        self.sections = np.rec.fromrecords(sections, names='typ,sec')
        for sec in self.sections['sec']:
            sec.v = -60

        # Connect sections
        for prev,next in zip(self.sections['sec'][:-1], self.sections['sec'][1:]):
            next.connect(prev)


        ### IHC Synapse
        self._syn = h.Exp2Syn(self.sections['sec'][0](0.5))
        self._syn.tau1 = 0.399806796048 / self._Tf(q10=2.4)
        self._syn.tau2 = 0.399889764048 / self._Tf(q10=2.4)
        assert self._syn.tau1 < self._syn.tau2
        self._syn.e = 0
        self._con = h.NetCon(None, self._syn)
        self._con.weight[0] = 0.000716352978448 * self._Tf(q10=1.5262)
        self._con.delay = 0

        ### Recording spikes from the last section
        last = self.sections['sec'][-1]
        self._probe = h.NetCon(last(0.5)._ref_v, None, 0, 0, 0, sec=last)
        self._spikes = h.Vector()
        self._probe.record(self._spikes)


        ### Recording voltages from all sections
        self._voltages = []
        if record_voltages:
            for sec in self.sections['sec']:
                vec = h.Vector()
                vec.record(sec(0.5)._ref_v)
                self._voltages.append(vec)



def _plot_voltages(voltages):
    subplt_num = voltages.shape[1]
    for i,v in enumerate(voltages.T):
        plt.subplot(subplt_num, 1, i+1)
        plt.plot(v)
    plt.show()



if __name__ == "__main__":
    import matplotlib.pyplot as plt

    h.dt = 0.002
    h.celsius = 37

    anf = ANF_Axon(record_voltages=True)

    h.topology()
    h.psection(sec=anf.sections['sec'][0])

    anf.vesicles = [2, 5]

    neuron.init()
    anf.ainit()
    neuron.run(10)

    _plot_voltages( anf.get_voltages()[:,0:6] )

    print "Spikes:", anf.get_spikes()


    print
    print "==============================================="
    print "Electrical stimulation"
    from electrodes import Electrode

    # set ANF
    anf = ANF_Axon(record_voltages=True)
    # anf.set_geometry('straight', x0=250, y0=500, z0=0)
    anf.set_geometry('bent', a=750, b=500, z=0)

    # set electrode
    el = Electrode()
    el.z = 0

    stim = np.zeros(1000)
    stim[280:300] = -0.5
    stim[300:320] = 0.5

    el.fs = 200000
    el.stim = stim

    # `connect' ANF and electrode
    anf.electrodes = [el]

    # run
    anf.einit()
    neuron.init()
    neuron.run(len(stim) * h.dt)

    # plot
    _plot_voltages(anf.get_voltages()[:,0:6])


