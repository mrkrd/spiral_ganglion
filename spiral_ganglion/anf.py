#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import logging
import warnings

import numpy as np
import pandas as pd

from neuron import h

from spiral_ganglion.electrodes import Electrode

logger = logging.getLogger(__name__)


def calc_tf(q10, temp_ref=22):
    tf = q10 ** ((h.celsius - temp_ref)/10)
    return tf


def _is_voltage_too_high(voltage):

    v = np.array(voltage) * 1e-3  # mV -> V

    if np.any(-0.16 > v) or np.any(v > 0.1):
        warnings.warn("Voltage was too high!")
        too_high = True
    else:
        too_high = False

    return too_high


class ANF(object):

    def get_sections(self, name=None):
        selected = []
        for typ, sec in zip(self.section_names, self.sections):
            if (name is None) or (typ == name):
                selected.append(sec)
        return selected

    def _get_segments(self):
        """Returns list of all segments along the neuron."""
        sections = self.sections
        segments = []
        for sec in sections:
            idx = [seg.x for seg in sec]
            segments.extend([sec(i) for i in idx])
        return segments

    def get_segment_path_positions(self, section_type=None):
        """Returns a list of positions of all segments along the
        neuron in meter (m).

        Parameters
        ----------
        section_type : {None, 'node', 'internode'}, optional
            Return segment path positions of the given section type
            only.

        Returns
        -------
        ndarray
            Segment positions.


        TODO: implement using h.Section().distance()

        """
        positions = []
        start = 0               # beginning of the currenct section
        for typ, sec in zip(self.section_names, self.sections):
            if (section_type is None) or (typ == section_type):
                seg_x = np.array([seg.x for seg in sec])    # [1]
                positions.extend(start + seg_x*sec.L*1e-6)  # um -> m
                start = start + 1e-6*sec.L                  # um -> m

        positions = np.array(positions)

        return positions

    def set_geometry(self, geometry, **kwargs):
        """Sets position of each segment of the neuron.

        Parameters
        ----------
        geometry : {'straight', 'bent'}
            Type of the geometry to use.

        'straight' (straight neuron along the x-axis)
        x0, y0, z0: coordinates of the initial segment (m)
        x0, y0: are in the plane of cross section of scala tympani
        z0: coordinate along cochlea

         y ^
           |
           |     o==o==o==o==o==o==...
           | (x0,y0,z0)
           |
           +---------------->
        (0,0,0)             x


        'bent'
        a, b: distances to the neuron along x and y axis respectively
        z: position along cochlea

        """
        self._geometry = geometry
        self._geometry_pars = kwargs

    def get_positions(self, section_type=None):
        """Return position of all segments of the neuron.

        Parameters
        ----------
        section_type : {None, 'node', 'internode'}, optional
            Return segment path positions of the given section type
            only.

        Returns
        -------
        ndarray
            (x,y,x) coordinates of every segment in the model.

        """
        ppos = self.get_segment_path_positions(section_type=section_type)

        if self._geometry is None:
            raise RuntimeError(
                "Neuron's geometry not set. Call anf.set_geometry() first."
            )

        elif self._geometry == 'straight':
            x = self._geometry_pars['x0'] + ppos
            y = self._geometry_pars['y0'] * np.ones(len(ppos))
            z = self._geometry_pars['z0'] * np.ones(len(ppos))

        elif self._geometry == 'bent':
            a = self._geometry_pars['a']
            b = self._geometry_pars['b']
            z_par = self._geometry_pars['z']
            r = b
            shift = np.abs(a-b)

            circumference = 2 * np.pi * r
            quarter = circumference / 4

            # Initial straight segment
            sel = ppos[(ppos < shift)]
            x = sel
            y = b * np.ones(len(sel))

            # Quarter circle segment
            sel = ppos[(ppos >= shift) & (ppos < (quarter+shift))] - shift
            x = np.append(x, shift + r*np.sin(2 * np.pi * sel / circumference))
            y = np.append(y, r*np.cos(2 * np.pi * sel / circumference))

            # Remaining straight segment
            sel = ppos[(ppos >= (quarter+shift))] - shift - quarter

            x = np.append(x, a*np.ones(len(sel)))
            y = np.append(y, -sel)
            z = z_par * np.ones(len(ppos))

            # import matplotlib.pyplot as plt
            # plt.plot(x,y, 'o')
            # plt.show()

        positions = np.rec.fromarrays(
            [x, y, z],
            names=['x', 'y', 'z']
        )
        return positions

    def get_voltages(self):
        """Return time courses of membrane potentials for each section of
        the model.  Potentials are recorded in the from the middle of
        each section: sec(0.5).

        The neuron must be initialized with `record_voltages=True'.

        """
        voltages = [np.asarray(vec) for vec in self._voltages]
        voltages = np.array(voltages).T

        voltages *= 1e-3        # mV -> V

        return voltages

    def get_trains(self):
        """Return an array of spike trains recorded from the last
        section.

        """
        if h.t == 0:
            raise RuntimeError("Time is 0 (did you run the simulation?)")

        if _is_voltage_too_high(self._last_voltage):
            train = None
        else:
            row = {
                'spikes': 1e-3*np.array(self._spikes),
                'duration': 1e-3*h.t,
                'type': 'anf'
            }

            row.update(self._meta)

            train = pd.DataFrame([row])

        return train

    def get_spikes(self):
        """Get raw spike times.  See also `get_trains` for spike_trains
        output.

        """
        if h.t == 0:
            raise RuntimeError("Time is 0 (did you run the simulation?)")

        if _is_voltage_too_high(self._last_voltage):
            spikes = None
        else:
            spikes = 1e-3*np.array(self._spikes)

        return spikes

    def ainit(self):
        """Initializes acoustical stimulation.

        Note: must be called *after* neuron.init()

        """
        assert h.dt <= 0.01

        for sec in self.sections:
            sec.v = self._vrest

        for v in self.vesicles:
            self._con.event(float(v*1e3))

    def einit(self, dt_assert=0.002):
        """Initializes electrical stimulation.

        Note: must be called *before* neuron.init()

        """
        if h.dt > dt_assert:
            raise ValueError("h.dt should be smaller than {dt_assert}, current value is {dt}".format(dt_assert=dt_assert, dt=h.dt))

        for sec in self.sections:
            sec.v = self._vrest

        self._stim_vectors = []
        if self.electrodes:
            # Make sure all electrodes have signals of the same fs
            fss = np.array([el.fs for el in self.electrodes])
            fs = fss[0]
            assert fs is not None
            assert np.all(fss == fs)
            el_dt = float(1000/fs)  # Hz -> ms

            potentials = sum([
                el.calculate_potentials(self) for el in self.electrodes
            ])

            segments = self._get_segments()

            assert len(potentials) == len(segments)

            for pot, seg in zip(potentials, segments):
                vec = h.Vector(pot*1e3)  # V -> mV
                vec.play(seg._ref_e_extracellular, el_dt)
                self._stim_vectors.append(vec)


def generate_anf_channels(
        name,
        diam=1.5
):
    capacitance = 0.0714e-12

    # Conductance is written as:
    # (single_channel_conductance * number_of_channels)
    # Ref: Negm and Bruce (2008) Table I
    g_na = calc_conductivity_cm2(
        25.69e-12 * 1000,
        capacitance
    )
    g_kv = calc_conductivity_cm2(
        50.0e-12 * 166,
        capacitance
    )
    g_klt = calc_conductivity_cm2(
        13.0e-12 * 166,
        capacitance
    )
    g_h = calc_conductivity_cm2(
        13.0e-12 * 100,
        capacitance
    )
    g_pas = 1e-5

    ena = 66
    ek = -88

    cfg = {}

    if name == 'schwarz1987_klt':
        vrest = -63.
        e_pas = -78.
        cfg['node_channels'] = [
            'na_schwarz1987',
            'k_schwarz1987',
            'klt_rothman2003',
            'ih_rothman2003',
            'pas',
            'extracellular'
        ]
        cfg['node_vars'] = {
            'nseg': 1,
            'L': 1,
            'Ra': 100,
            'diam': diam,
            'cm': 0.9,
            'gnabar_na_schwarz1987': g_na,
            'gkbar_k_schwarz1987': g_kv,
            'gkltbar_klt_rothman2003': g_klt,
            'ghbar_ih_rothman2003': g_h,
            'g_pas': g_pas,
            'ena': ena,
            'ek': ek,
            'e_pas': e_pas,
        }
        cfg['internode_channels'] = [
            'pas',
            'extracellular'
        ]
        cfg['internode_vars'] = {
            'nseg': 1,
            'L': 250,
            'Ra': 100,
            'diam': diam,
            'cm': 1e-3,
            'g_pas': g_pas,
            'e_pas': e_pas,
        }
        cfg['global_vars'] = {
            'vrest_na_schwarz1987': vrest,
            'vrest_k_schwarz1987': vrest,
            'vrest_klt_rothman2003': vrest,
            'vrest_ih_rothman2003': vrest,
        }
        cfg['vrest'] = vrest

    elif name == 'schwarz1987_pure':
        vrest = -78.
        cfg['node_channels'] = [
            'na_schwarz1987',
            'k_schwarz1987',
            'pas',
            'extracellular'
        ]
        cfg['node_vars'] = {
            'nseg': 1,
            'L': 1,
            'Ra': 100,
            'diam': diam,
            'cm': 0.9,
            'gnabar_na_schwarz1987': g_na,
            'gkbar_k_schwarz1987': g_kv,
            'g_pas': g_pas,
            'ena': ena,
            'ek': ek,
            'e_pas': vrest
        }
        cfg['internode_channels'] = [
            'pas',
            'extracellular'
        ]
        cfg['internode_vars'] = {
            'nseg': 1,
            'L': 250,
            'Ra': 100,
            'diam': diam,
            'cm': 1e-3,
            'g_pas': g_pas,
            'e_pas': vrest
        }
        cfg['global_vars'] = {
            'vrest_na_schwarz1987': vrest,
            'vrest_k_schwarz1987': vrest,
        }
        cfg['vrest'] = vrest

    elif name == 'passive_only':
        vrest = -78
        cfg['node_channels'] = [
            'pas',
            'extracellular'
        ]
        cfg['node_vars'] = {
            'nseg': 1,
            'L': 1,
            'Ra': 100,
            'diam': diam,
            'cm': 0.9,
            'g_pas': g_pas,
            'e_pas': vrest
        }
        cfg['internode_channels'] = [
            'pas',
            'extracellular'
        ]
        cfg['internode_vars'] = {
            'nseg': 1,
            'L': 250,
            'Ra': 100,
            'diam': diam,
            'cm': 1e-3,
            'g_pas': g_pas,
            'e_pas': vrest
        }
        cfg['global_vars'] = {
        }
        cfg['vrest'] = -78

    elif name == 'rothman1993_klt':
        cfg['node_channels'] = [
            'na_rothman1993',
            'kht_rothman2003',
            'klt_rothman2003',
            'ih_rothman2003',
            'pas'
        ]
        cfg['node_vars'] = {
            'nseg': 1,
            'L': 1,
            'Ra': 100,
            'diam': diam,
            'cm': 0.9,
            'gnabar_na_rothman1993': g_na,
            'gkhtbar_kht_rothman2003': g_kv,
            'gkltbar_klt_rothman2003': g_klt,
            'ih_rothman2003': g_h,
            'g_pas': g_pas,
            'ena': ena,
            'ek': ek,
            'e_pas': vrest
        }
        cfg['internode_channels'] = [
            'pas',
            'extracellular'
        ]
        cfg['internode_vars'] = {
            'nseg': 1,
            'L': 250,
            'Ra': 100,
            'diam': diam,
            'cm': 1e-3,
            'g_pas': g_pas,
            'e_pas': vrest
        }
        cfg['global_vars'] = {
            'vrest_na_rothman1993': -64,
            'vrest_kht_rothman2003': -64,
            'vrest_klt_rothman2003': -64,
            'vrest_ih_rothman2003': -64
        }
        cfg['vrest'] = -64

    else:
        raise RuntimeError("Unknown channel config name: {}".format(name))

    return cfg


class ANF_Axon(ANF):
    """Model of spiral ganglion's peripherial axon.  It can be used for
    acoustical and electrical stimulation.

    Parameters
    ----------
    nodes : int
        Number of nodes in the model. Total number of compartments
        (nodes + internodes) is 2*nodes.
    record_voltages : bool
        When True membrane potentials are recorded internally and
        can be returned with `get_voltages`.
    weight : float, optional
        Synaptic weight.
    meta : dict, optional
        Meta data that will be added to the output spike train.


    Example
    -------
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
            nodes=20,
            record_voltages=False,
            channels='schwarz1987_pure',
            terminal_length=10,
            terminal_nseg=1,
            diam=1.5,
            weight=None,
            meta=None,
    ):
        logger.info("ANF temperature: {} C".format(h.celsius))

        # vesicle timings for acoustical stimulation
        self.vesicles = []

        self._geometry = None

        # electrodes that stimulate the neuron (class Electrode)
        self.electrodes = []

        if not isinstance(channels, dict):
            cfg = generate_anf_channels(
                channels,
                diam=diam
            )
        else:
            cfg = channels

        sections = []
        names = []
        for i in range(nodes):
            # Node sections
            sec = h.Section()

            for chan in cfg['node_channels']:
                sec.insert(chan)
            for var, val in cfg['node_vars'].items():
                setattr(sec, var, val)

            sections.append(sec)
            names.append('node')

            # Internode sections
            sec = h.Section()

            for chan in cfg['internode_channels']:
                sec.insert(chan)
            for var, val in cfg['internode_vars'].items():
                setattr(sec, var, val)

            sections.append(sec)
            names.append('internode')

        for var, val in cfg['global_vars'].items():
            setattr(h, var, val)

        # Terminal node
        sections[0].L = terminal_length
        sections[0].nseg = terminal_nseg

        for sec in sections:
            sec.v = cfg['vrest']
        self._vrest = cfg['vrest']

        # Connect sections
        for a, b in zip(sections[:-1], sections[1:]):
            b.connect(a)

        # IHC Synapse
        self._syn = h.Exp2Syn(sections[0](0.5))
        self._syn.tau1 = 0.000392876886274 * 1e3 / calc_tf(q10=2.4)
        self._syn.tau2 = 0.000392990435175 * 1e3 / calc_tf(q10=2.4)
        assert self._syn.tau1 < self._syn.tau2
        self._syn.e = 0
        self._con = h.NetCon(None, self._syn)
        self._con.delay = 0
        if weight is None:
            # Weight at 22C from an EPSC fit to Li et al. (2009)
            # (spiral_ganglion_notes/li2009_epsc_fit/run_epsc_fit.py);
            # Q10 from rate fit at 37C
            # (spiral_ganglion_notes/acoustic_q10_weight/find_best_weight.py)
            self._con.weight[0] = 1e6 * 5.5686753127e-10 * calc_tf(q10=1.25187741383)
        else:
            self._con.weight[0] = 1e6 * weight

        # Recording spikes from the last section
        last = sections[-1]
        self._probe = h.NetCon(last(0.5)._ref_v, None, 0, 0, 0, sec=last)
        self._spikes = h.Vector()
        self._probe.record(self._spikes)

        # Record voltages from the last section
        last_voltage = h.Vector()
        last_voltage.record(
            sections[-1](0.5)._ref_v
        )
        self._last_voltage = last_voltage

        # Recording voltages from all sections
        self._voltages = []
        if record_voltages:
            logger.info("Recording voltages is on")

            for sec in sections:
                for seg in sec:
                    vec = h.Vector()
                    vec.record(seg._ref_v)
                    self._voltages.append(vec)

                    # try:
                    #     vec.record(seg.na_schwarz1987._ref_m)
                    #     self._voltages.append(vec)
                    # except Exception:
                    #     pass

                # vec = h.Vector()
                # vec.record(sec(0.5)._ref_v)
                # self._voltages.append(vec)

        assert len(sections) == len(names)
        self.sections = sections
        self.section_names = names

        # Meta data for spike trains
        if meta is None:
            self._meta = {'type': 'anf'}
        else:
            self._meta = meta


def calc_conductivity_cm2(conductance, capacitance):
    cm = 0.9e-6                 # [F/cm2]
    area = capacitance / cm        # [cm2]

    conductivity = conductance / area  # [S/cm2]
    return conductivity


def plot_vectors(vectors, fs):
    """TODO: use single axis, move to thorns"""

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(
        vectors.shape[1],
        1,
        sharex=True,
        sharey=True
    )

    for v, a in zip(vectors.T, ax):
        lines = a.plot(v)
        a.patch.set_visible(False)
        a.set_frame_on(False)
        # a.get_major_ticks().set_visible(False)

        for line in lines:
            line.set_clip_on(False)

    a.set_frame_on(True)
    # a.ax.get_xaxis().set_visible(True)

    return fig, ax


def plot_geometry(objects):
    """Plot `objects` (e.g. `ANF_Axon`, `Electrode`) onto XYZ planes.

    """
    import matplotlib.pyplot as plt

    fig = plt.figure()
    xy = fig.add_subplot(221)
    zy = fig.add_subplot(222, sharey=xy)
    xz = fig.add_subplot(223, sharex=xy)

    for obj in objects:
        if isinstance(obj, Electrode):
            fmt = 'k*'
            x = obj.x
            y = obj.y
            z = obj.z

        elif isinstance(obj, ANF):
            pos = obj.get_positions()
            fmt = 'ko'
            x = pos['x']
            y = pos['y']
            z = pos['z']

        _plot_object(x, y, z, fmt, fig)

    # return fig


def _plot_object(x, y, z, fmt, fig):

    xy, zy, xz = fig.axes

    xy.plot(x, y, fmt)
    xy.set_title("XY")
    xy.set_xlabel("X [m]")
    xy.set_ylabel("Y [m]")
    xy.set_aspect('equal', 'datalim')

    zy.plot(z, y, fmt)
    zy.set_title("ZY")
    zy.set_xlabel("Z [m]")
    zy.set_ylabel("Y [m]")
    zy.set_aspect('equal', 'datalim')

    xz.plot(x, z, fmt)
    xz.set_title("XZ")
    xz.set_xlabel("X [m]")
    xz.set_ylabel("Z [m]")
    xz.set_aspect('equal', 'datalim')
