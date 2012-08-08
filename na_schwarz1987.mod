
NEURON {
    SUFFIX na_schwarz1987
    USEION na READ ena WRITE ina
    RANGE gnabar, ina
    GLOBAL minf, hinf, mtau, htau
    THREADSAFE
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}

PARAMETER {
    gnabar=.120 (mho/cm2) <0,1e9>
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
    LOCAL Tf, vm
    Tf = 2.2^((celsius - 20(degC))/10(degC))
    vm = v + 78

    malpha  = 0.49 * Tf * expM1( (25.41-vm), 6.06 )
}

FUNCTION mbeta(v(mV)) (/ms) {
    LOCAL Tf, vm
    Tf = 2.2^((celsius - 20(degC))/10(degC))
    vm = v + 78

    mbeta = 1.04 * Tf * expM1( (vm-21.00), 9.41 )
}

FUNCTION halpha(v(mV)) (/ms) {
    LOCAL Tf, vm
    Tf = 2.9^((celsius - 20(degC))/10(degC))
    vm = v + 78

    halpha = 0.09 * Tf * expM1( (vm+27.74), 9.06 )
}

FUNCTION hbeta(v(mV)) (/ms) {
    LOCAL Tf, vm
    Tf = 2.9^((celsius - 20(degC))/10(degC))
    vm = v + 78

    hbeta = 3.70 * Tf / (1 + exp((56.0-vm)/12.50))
}

FUNCTION expM1(x,y) {
    if (fabs(x/y) < 1e-6) {
	expM1 = y*(1 - x/y/2)
    } else {
	expM1 = x/(exp(x/y) - 1)
    }
}

PROCEDURE rates(v(mV)) {LOCAL a, b
    TABLE minf, hinf, mtau, htau DEPEND celsius FROM -100 TO 100 WITH 200
    mtau = 1/(malpha(v) + mbeta(v))
    minf = malpha(v)/(malpha(v) + mbeta(v))

    htau = 1/(halpha(v) + hbeta(v))
    hinf = halpha(v)/(halpha(v) + hbeta(v))
}
