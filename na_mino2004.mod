
NEURON {
    SUFFIX na_mino2004
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
    LOCAL Tf
    Tf = 3^((celsius - 37(degC))/10(degC))

    malpha  = 1.872 * Tf * expM1( (25.41-v), 6.06 )
    :malpha = 1.872 * Tf * (25.41-v) / (exp((25.41-v)/6.06) - 1)
    :: x = (25.41-v)
    :: y = 6.06
}

FUNCTION mbeta(v(mV)) (/ms) {
    LOCAL Tf
    Tf = 3^((celsius - 37(degC))/10(degC))

    mbeta = 0.4 * Tf * expM1( (v-21.001), 9.41 )
    :mbeta = 3.973 * Tf * (v-21.001) / (exp((v-21.001)/9.41) - 1)
}

FUNCTION halpha(v(mV)) (/ms) {
    LOCAL Tf, T10
    Tf = 3^((celsius - 37(degC))/10(degC))

    halpha = 0.549 * Tf * expM1( (v+27.74), 9.06 )
}

FUNCTION hbeta(v(mV)) (/ms) {
    LOCAL Tf
    Tf = 3^((celsius - 37(degC))/10(degC))

    hbeta = 22.57 * Tf / (1 + exp((56.0-v)/12.5))
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
