
NEURON {
    SUFFIX k_mino2004
    USEION k READ ek WRITE ik
    RANGE gkbar, ik
    GLOBAL ninf, ntau
    THREADSAFE
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}

PARAMETER {
    gkbar=.120 (mho/cm2) <0,1e9>
}

STATE {
    n
}

ASSIGNED {
    v (mV)
    celsius (degC)
    ek (mV)
    ik (mA/cm2)
    ninf
    ntau (ms)
}

INITIAL {
    rates(v)
    n = ninf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gkbar*n*n*n*n*(v - ek)
}

DERIVATIVE states {
    rates(v)
    n' = (ninf - n)/ntau
}

FUNCTION nalpha(v(mV)) (/ms) {
    LOCAL Tf
    Tf = 3^((celsius - 37(degC))/10(degC))

    nalpha  = 0.129 * Tf * expM1( (35.0-v), 10.0 )
}

FUNCTION nbeta(v(mV)) (/ms) {
    LOCAL Tf
    Tf = 3^((celsius - 37(degC))/10(degC))

    nbeta = 0.3236 * Tf * expM1( (v-35.0), 10.0 )
}

FUNCTION expM1(x,y) {
    if (fabs(x/y) < 1e-6) {
	expM1 = y*(1 - x/y/2)
    } else {
	expM1 = x/(exp(x/y) - 1)
    }
}

PROCEDURE rates(v(mV)) {LOCAL a, b
    TABLE ninf, ntau DEPEND celsius FROM -100 TO 100 WITH 200
    ntau = 1/(nalpha(v) + nbeta(v))
    ninf = nalpha(v)/(nalpha(v) + nbeta(v))
}
