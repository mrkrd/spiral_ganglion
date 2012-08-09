#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

import neuron
from neuron import h

from electrodes import Electrode

class ANF(object):
    def _Tf(self, q10, temp_ref=22):
        return q10 ** ((h.celsius - temp_ref)/10)


    def _get_segments(self):
        """ Returns list of all segments along the neuron """
        sections = self.sections['sec']
        segments = []
        for sec in sections:
            idx = [seg.x for seg in sec]
            segments.extend([sec(i) for i in idx])
        return segments


    def _get_segment_positions_along_neuron(self):
        """Returns a list of positions of all segments along the
        neuron in [m]

        """
        positions = []
        start = 0               # beginning of the currenct section
        for sec in self.sections['sec']:
            seg_x = np.array([seg.x for seg in sec]) # [1]
            positions.extend(start + seg_x*1e-6*sec.L) # um -> m
            start = start + 1e-6*sec.L                 # um -> m

        positions = np.array(positions)

        return positions


    def set_geometry(self, geometry, **kwargs):
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
        Ls = self._get_segment_positions_along_neuron()

        if geometry == 'straight':
            self.x = kwargs['x0'] + Ls
            self.y = kwargs['y0'] * np.ones(len(Ls))
            self.z = kwargs['z0'] * np.ones(len(Ls))

        elif geometry == 'bent':
            a = kwargs['a']
            b = kwargs['b']
            z = kwargs['z']
            r = b
            shift = np.abs(a-b)


            circumference = 2 * np.pi * r
            quarter = circumference / 4


            # Initial straight segment
            sel = Ls[ (Ls<shift) ]
            x = sel
            y = b * np.ones(len(sel))

            # Quarter circle segment
            sel = Ls[ (Ls>=shift) & (Ls<(quarter+shift)) ] - shift
            x = np.append(x, shift + r*np.sin(2 * np.pi * sel / circumference))
            y = np.append(y, r*np.cos(2 * np.pi * sel / circumference))

            # Remaining straight segment
            sel = Ls[ (Ls>=(quarter+shift)) ] - shift - quarter
            x = np.append(x, a*np.ones(len(sel)))
            y = np.append(y, -sel)

            self.x = x
            self.y = y
            self.z = kwargs['z'] * np.ones(len(Ls))

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

        voltages *= 1e-3        # mV -> V

        return voltages


    def get_spikes(self):
        """Return an array of spike timings recorded from the last
        section.

        """
        assert h.t != 0, "Time is 0 (did you run the simulation already?)"

        train = np.array(
            [(1e-3*np.array(self._spikes), 1e-3*h.t)], # ms -> s
            dtype=[('spikes', np.ndarray),
                   ('duration', float)]
            )

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



    def einit(self, dt_assert=0.002, active_nodes=15):
        """
        Initializes electrical stimulation.

        Note: must be called *before* neuron.init()

        """
        assert h.dt <= dt_assert

        for sec in self.sections['sec']:
            sec.v = -60

        self._stim_vectors = []
        if self.electrodes:
            # Make sure that each segment has defined position
            assert (
                len(self.x) ==
                len(self.y) ==
                len(self.z) ==
                len(self._get_segment_positions_along_neuron())
            )

            # Make sure all electrodes have signals of the same fs
            fss = np.array([el.fs for el in self.electrodes])
            fs = fss[0]
            assert fs is not None
            assert np.all(fss == fs)
            el_dt = float(1000/fs) # Hz -> ms

            potentials = sum([
                el.calculate_potentials(self) for el in self.electrodes
            ])


            for pot,seg in zip(potentials, self._get_segments()):
                vec = h.Vector(pot*1e3) # V -> mV
                vec.play(seg._ref_e_extracellular, el_dt)
                self._stim_vectors.append(vec)





class ANF_Axon(ANF):
    """
    Model of Auditory Nerve Fiber's peripherial axon.  Can be used
    for acoustical and electrical stimulation.

    anf = ANF_Axon()

    ### Acoustical stimulation
    anf.vesicles = [12, 24, 55]
    neuron.init()
    anf.ainit()
    neuron.run(60)


    ### Electrical stimulation
    anf.electrodes = [electrode1, electrode2, electrode3]
    anf.einit()
    neuron.init()
    neuron.run(60)

    """
    def __init__(
            self,
            nodes=15,
            record_voltages=False,
            debug=True):
        """nodes: number of nodes in the model.  Total number of
               compartments if 2*nodes.

        record_voltages: when True membrane potentials are recorded
                         internally and can be returned with
                         get_voltages()

        """
        if debug:
            print "ANF temperature:", h.celsius, "C"


        self.vesicles = []      # vesicle timings for acoustical stimulation
        self.x = None          # array of segment's x coordinate locations
        self.y = None          # array of segment's y coordinate locations
        self.z = None          # array of segment's z coordinate locations

        self.electrodes = []    # electrodes that stimulate the neuron
                                # (class Electrode)

        ena = 66
        ek = -88
        gna = 0.324
        gpas = 1e-5

        sections = []

        ### Peripherial Axon Terminal
        term = h.Section()
        term.nseg = 1
        term.L = 10
        term.Ra = 100
        term.diam = 1.5
        term.cm = 0.9

        term.insert('pas')
        term.g_pas = gpas

        term.insert('na_schwarz1987')
        term.gnabar_na_schwarz1987 = gna

        term.insert('k_schwarz1987')
        term.gkbar_k_schwarz1987 = 0.105

        term.insert('klt_manis')
        term.gkltbar_klt_manis = 0.027
        term.insert('ih_manis')
        term.ghbar_ih_manis = 0.016
        term.insert('extracellular')

        term.ena = ena
        term.ek = ek

        sections.append( ('p_term', term) )


        ### Peripherial Axon Nodes and Internodes
        for i in range(nodes):
            inode = h.Section()
            inode.nseg = 1
            inode.L = 250
            inode.Ra = 100
            inode.diam = 1.5
            inode.cm = 1e-3

            inode.insert('pas')
            inode.g_pas = gpas

            inode.insert('extracellular')

            sections.append( ('p_inode', inode) )




            node = h.Section()
            node.nseg = 1
            node.L = 1
            node.Ra = 100
            node.diam = 1.5
            node.cm = 0.9

            node.insert('pas')
            node.g_pas = gpas

            node.insert('na_schwarz1987')
            node.gnabar_na_schwarz1987 = gna

            node.insert('k_schwarz1987')
            node.gkbar_k_schwarz1987 = 0.105

            node.insert('klt_manis')
            node.gkltbar_klt_manis = 0.027

            node.insert('ih_manis')
            node.ghbar_ih_manis = 0.016

            node.insert('extracellular')

            node.ena = ena
            node.ek = ek

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
            print "Recording voltages is on"
            for sec in self.sections['sec']:
                vec = h.Vector()
                vec.record(sec(0.5)._ref_v)
                self._voltages.append(vec)



def _plot_voltages(voltages):

    import matplotlib.pyplot as plt

    fig,axes = plt.subplots(
        voltages.shape[1],
        1,
        sharex=True,
        sharey=True
    )

    for v,a in zip(voltages.T, axes):
        a.plot(v)



def plot_geometry(objects):
    import matplotlib.pyplot as plt

    fig = plt.figure()

    for obj in objects:
        if isinstance(obj, Electrode):
            fmt = 'k*'
            x = obj.x
            y = obj.y
            z = obj.z

        elif isinstance (obj, ANF):
            fmt = 'ko'
            x = obj.x
            y = obj.y
            z = obj.z


        _plot_object(x,y,z,fmt, fig)



def _plot_object(x,y,z,fmt, fig):

    xy = fig.add_subplot(221)
    xy.plot(x,y,fmt)
    xy.set_title("XY")
    xy.set_xlabel("X [m]")
    xy.set_ylabel("Y [m]")
    xy.set_aspect('equal', 'datalim')



    zy = fig.add_subplot(222, sharey=xy)
    zy.plot(z,y,fmt)
    zy.set_title("ZY")
    zy.set_xlabel("Z [m]")
    zy.set_ylabel("Y [m]")
    zy.set_aspect('equal', 'datalim')


    xz = fig.add_subplot(223, sharex=xy)
    xz.plot(x,z,fmt)
    xz.set_title("XZ")
    xz.set_xlabel("X [m]")
    xz.set_ylabel("Z [m]")
    xz.set_aspect('equal', 'datalim')


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
    # anf.set_geometry('straight', x0=250e-6, y0=500e-6, z0=0)
    anf.set_geometry('bent', a=750e-6, b=500e-6, z=0)

    # set electrode
    el = Electrode()
    el.z = 0

    fs = 200e3                  # [Hz]
    stim = np.zeros(20e-3 * fs)
    stim[np.round(5e-3*fs):np.round(6e-3*fs)] = -0.1e-3 # [A]

    el.fs = fs
    el.stim = stim

    # `connect' ANF and electrode
    anf.electrodes = [el]

    # run
    anf.einit()
    neuron.init()
    neuron.run(len(stim) * h.dt)

    # plot
    print anf.get_spikes()
    _plot_voltages(anf.get_voltages()[:,0:6])




    plot_geometry( [anf,el] )
    plt.show()
