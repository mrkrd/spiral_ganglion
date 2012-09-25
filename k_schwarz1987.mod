
NEURON {
    SUFFIX k_schwarz1987
    USEION k READ ek WRITE ik
    RANGE gkbar, ik
    GLOBAL ninf, ntau, vrest
    THREADSAFE
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}

PARAMETER {
    gkbar=.120 (mho/cm2) <0,1e9>
    vrest=-78. (mV)
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
    Tf = 3.0^((celsius - 20(degC))/10(degC))

    nalpha  = 0.02 * Tf * expM1( (35.0-v), 10.0 )
}

FUNCTION nbeta(v(mV)) (/ms) {
    LOCAL Tf
    Tf = 3.0^((celsius - 20(degC))/10(degC))

    nbeta = 0.05 * Tf * expM1( (v-10.0), 10.0 )
}

FUNCTION expM1(x,y) {
    if (fabs(x/y) < 1e-6) {
	expM1 = y*(1 - x/y/2)
    } else {
	expM1 = x/(exp(x/y) - 1)
    }
}

PROCEDURE rates(v(mV)) {
    LOCAL vm
    TABLE ninf, ntau DEPEND celsius FROM -100 TO 100 WITH 200

    vm = v - vrest

    ntau = 1/(nalpha(vm) + nbeta(vm))
    ninf = nalpha(vm)/(nalpha(vm) + nbeta(vm))
}
