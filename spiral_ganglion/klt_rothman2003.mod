TITLE klt_rothman2003.mod  The low threshold conductance of cochlear nucleus neurons

COMMENT

Neuron implementation of Jason Rothman's measurements of VCN conductances.

This file implements the low threshold potassium current found in several brainstem
nuclei of the auditory system, including the spherical and globular bushy cells
(Manis and Marx, 1991; Rothman and Manis, 2003a,b) and octopus cells (Bal and
Oertel, 2000) of the ventral cochlear nucleus, principal cells of the medial
nucleus of the trapzoid body (Brew and Forsythe, 1995, Wang and Kaczmarek,
1997) and neurons of the medial superior olive. The current is likely mediated by
heteromultimers of Kv1.1 and Kv1.2 potassium channel subunits. The specific
implementation is described in Rothman and Manis, J. Neurophysiol. 2003, in the
appendix. Measurements were made from isolated neurons from adult guinea pig,
under reasonably stringent voltage clamp conditions. The measured current is
sensitive to the mamba snake toxin dendrotoxin-I.


Similar conductrances are found in the homologous neurons of the avian auditory
system (Reyes and Rubel; Zhang and Trussell; Rathouz and Trussell), and the
conductance described here, in the absence of more detailed kinetic measurements
, is probably suitable for use in modeling that system.


Original implementation by Paul B. Manis, April (JHU) and Sept, (UNC)1999.

File split implementation, February 28, 2004.

Contact: pmanis@med.unc.edu

ENDCOMMENT

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (nA) = (nanoamp)
}

NEURON {
    SUFFIX klt_rothman2003
    USEION k READ ek WRITE ik
    RANGE gkltbar, ik
    GLOBAL winf, zinf, wtau, ztau, vrest
}


PARAMETER {
    v (mV)
    celsius = 22 (degC)  : model is defined on measurements made at room temp in Baltimore
    dt (ms)
    ek = -77 (mV)
    gkltbar = 0.01592 (mho/cm2) <0,1e9>
    zss = 0.5   <0,1>   : steady state inactivation of glt
    vrest = -64. (mV)
}

STATE {
    w z
}

ASSIGNED {
    ik (mA/cm2)
    winf zinf
    wtau (ms) ztau (ms)
}



BREAKPOINT {
    SOLVE states METHOD cnexp

    ik = gkltbar*w*w*w*w*z*(v - ek)
}


INITIAL {
    rates(v)
    w = winf
    z = zinf
}


DERIVATIVE states {
    rates(v)

    w' = (winf-w)/wtau
    z' = (zinf-z)/ztau
}



PROCEDURE rates(v(mV)) {
    LOCAL vm, Tf
    TABLE winf, wtau, zinf, ztau DEPEND celsius FROM -150 TO 150 WITH 300

    vm = v + (-64 - vrest)

    UNITSOFF

    Tf = 3^((celsius - 22)/10) : if you don't like room temp, it can be changed!

    winf = (1 / (1 + exp(-(vm+48) / 6)))^0.25
    zinf = zss + ((1-zss) / (1 + exp((vm + 71) / 10)))

    wtau = (100 / (6*exp((vm+60) / 6) + 16*exp(-(vm+60) / 45))) + 1.5
    wtau = wtau / Tf
    ztau = (1000 / (exp((vm+60) / 20) + exp(-(vm+60) / 8))) + 50
    ztau = ztau / Tf
}





FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    }else{
        vtrap = x/(exp(x/y) - 1)
    }
}

UNITSON
