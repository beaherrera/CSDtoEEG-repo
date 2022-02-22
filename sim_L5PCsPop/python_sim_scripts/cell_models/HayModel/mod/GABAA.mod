TITLE GABAA receptor


COMMENT
GABAA receptor conductance using a dual-exponential profile
presynaptic short-term plasticity based on Fuhrmann et al, 2002
Implemented by Beatriz Herrera April 8 based on the implementation by
Srikanth Ramaswamy, Blue Brain Project, March 2009 and Guy 2018.
Some modificaations of Ramaswamy implementation were made based on
Mäki-Marttunen et al., J of Neurosci Methods. 293, 264-283 (2018).
GUY: Removed  plasticity and depression

ENDCOMMENT


NEURON {

        POINT_PROCESS GABAA
        RANGE  tau_r_GABAA, tau_d_GABAA
        RANGE Use
        RANGE i,  i_GABAA,  g_GABAA, e, gmax
        NONSPECIFIC_CURRENT i
}

PARAMETER {

	   tau_r_GABAA = 0.2   (ms) : dual-exponential conductance profile
        tau_d_GABAA = 1.7     (ms) : IMPORTANT: tau_r < tau_d
        Use = 1.0   (1)   : Utilization of synaptic efficacy (just initial values! Use, Dep and Fac are overwritten by BlueBuilder assigned values)

        e = -80     (mV)  : GABAA reversal potential
    	:gmax = .001 (uS) :1nS weight conversion factor (from nS to uS)
    	u0 = 0 :initial value of u, which is the running value of Use
}

COMMENT
The Verbatim block is needed to generate random nos. from a uniform distribution between 0 and 1
for comparison with Pr to decide whether to activate the synapse or not
ENDCOMMENT




ASSIGNED {

        v (mV)
        i (nA)
	i_GABAA (nA)
	g_GABAA (uS)
	factor_GABAA

}

STATE {


	A_GABAA       : GABAA state variable to construct the dual-exponential profile - decays with conductance tau_r_GABAA
    B_GABAA       : GABAA state variable to construct the dual-exponential profile - decays with conductance tau_d_GABAA
}

INITIAL{

    LOCAL  tp_GABAA

	A_GABAA = 0
	B_GABAA = 0

	tp_GABAA = (tau_r_GABAA*tau_d_GABAA)/(tau_d_GABAA-tau_r_GABAA)*log(tau_d_GABAA/tau_r_GABAA) :time to peak of the conductance



	factor_GABAA = -exp(-tp_GABAA/tau_r_GABAA)+exp(-tp_GABAA/tau_d_GABAA) :GABAA Normalization factor - so that when t = tp_GABAA, gsyn = gpeak
    factor_GABAA = 1/factor_GABAA

}

BREAKPOINT {

    SOLVE state METHOD cnexp
    g_GABAA = (B_GABAA-A_GABAA) :compute time varying conductance as the difference of state variables B_GABAA and A_GABAA and mggate kinetics

	i_GABAA = g_GABAA*(v-e) :compute the GABAA driving force based on the time varying conductance, membrane potential, and GABAA reversal
	i =  i_GABAA
}

DERIVATIVE state{


	A_GABAA' = -A_GABAA/tau_r_GABAA
    B_GABAA' = -B_GABAA/tau_d_GABAA
}





NET_RECEIVE (weight, weight_GABAA){

	weight_GABAA = weight



	A_GABAA = A_GABAA + weight_GABAA*factor_GABAA
    B_GABAA = B_GABAA + weight_GABAA*factor_GABAA


}
