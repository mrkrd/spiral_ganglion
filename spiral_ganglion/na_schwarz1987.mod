NEURON {
    SUFFIX na_schwarz1987
    USEION na READ ena WRITE ina
    RANGE gnabar, ina
    GLOBAL minf, hinf, mtau, htau, vrest
    THREADSAFE
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}

PARAMETER {
    gnabar = .120 (mho/cm2) <0,1e9>
    vrest = -78. (mV)
}

STATE {
    m h
}

ASSIGNED {
    v (mV)
    celsius (degC)
    ena (mV)
    ina (mA/cm2)
    minf hinf
    mtau (ms)
    htau (ms)
}

INITIAL {
    rates(v)
    m = minf
    h = hinf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ina = gnabar*m*m*m*h*(v - ena)
}

DERIVATIVE states {
    rates(v)
    m' = (minf - m)/mtau
    h' = (hinf - h)/htau
}

FUNCTION malpha(v(mV)) (/ms) {
    LOCAL Tf
    Tf = 2.2^((celsius - 20(degC))/10(degC))

    malpha  = 0.49 * Tf * expM1( (25.41-v), 6.06 )
}

FUNCTION mbeta(v(mV)) (/ms) {
    LOCAL Tf
    Tf = 2.2^((celsius - 20(degC))/10(degC))

    mbeta = 1.04 * Tf * expM1( (v-21.00), 9.41 )
}

FUNCTION halpha(v(mV)) (/ms) {
    LOCAL Tf
    Tf = 2.9^((celsius - 20(degC))/10(degC))

    halpha = 0.09 * Tf * expM1( (v+27.74), 9.06 )
}

FUNCTION hbeta(v(mV)) (/ms) {
    LOCAL Tf
    Tf = 2.9^((celsius - 20(degC))/10(degC))

    hbeta = 3.70 * Tf / (1 + exp((56.0-v)/12.50))
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
    TABLE minf, hinf, mtau, htau DEPEND celsius FROM -100 TO 100 WITH 200

    vm = v - vrest

    mtau = 1/(malpha(vm) + mbeta(vm))
    minf = malpha(vm)/(malpha(vm) + mbeta(vm))

    htau = 1/(halpha(vm) + hbeta(vm))
    hinf = halpha(vm)/(halpha(vm) + hbeta(vm))
}
