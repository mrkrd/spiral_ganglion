TITLE ih_manis.mod  VCN conductances

COMMENT
Ih for VCN neurons - average from several studies in auditory neurons


Implementation by Paul B. Manis, April (JHU) and Sept, (UNC)1999.
revised 2/28/04 pbm

pmanis@med.unc.edu

ENDCOMMENT

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (nA) = (nanoamp)
}

NEURON {
    SUFFIX ih_manis
    NONSPECIFIC_CURRENT i
    RANGE ghbar, ih
    GLOBAL rinf, rtau, vrest
}


PARAMETER {
    v (mV)
    celsius = 22 (degC)
    dt (ms)
    ghbar = 0.00318 (mho/cm2) <0,1e9>
    eh = -43 (mV)
    vrest = -64. (mV)
}

STATE {
    r
}

ASSIGNED {
    i (mA/cm2)
    rinf
    rtau (ms)
}



BREAKPOINT {
    SOLVE states METHOD cnexp

    i = ghbar*r*(v - eh)
}



INITIAL {
    rates(v)
    r = rinf
}

DERIVATIVE states {
    rates(v)

    r' = (rinf-r)/rtau
}



PROCEDURE rates(v) {
    LOCAL vm, Tf
    TABLE rinf, rtau DEPEND celsius FROM -200 TO 150 WITH 350

    vm = v + (-64 - vrest)

    UNITSOFF

    Tf = 3^((celsius - 22)/10)

    rinf = 1 / (1+exp((vm + 76) / 7))

    rtau = (100000 / (237*exp((vm+60) / 12) + 17*exp(-(vm+60) / 14))) + 25
    rtau = rtau / Tf

}


FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    }else{
        vtrap = x/(exp(x/y) - 1)
    }
}

UNITSON
