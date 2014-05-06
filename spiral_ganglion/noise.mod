NEURON {
        SUFFIX noise
        NONSPECIFIC_CURRENT i : not in charge balance eq. and opposite sign -> ELECTRODE_CURRENT i 
	RANGE new_seed,sd
}

PARAMETER {
           i_randn (milliamp/cm2)
           sd = 0
}

ASSIGNED {
          i (milliamp/cm2)
}

BREAKPOINT { 
            SOLVE new_randn
            i = i_randn
            }

PROCEDURE new_randn() { 
            i_randn = normrand(0,sd)
            }

PROCEDURE new_seed(seed) {		: procedure to set the seed
	set_seed(seed)
	:VERBATIM
	:  printf("Setting random generator with seed = %g\n", _lseed);
	:ENDVERBATIM
}
