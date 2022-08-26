TITLE AMPA and NMDA receptor with presynaptic short-term plasticity


COMMENT
AMPA and NMDA receptor conductance using a dual-exponential profile
presynaptic short-term plasticity based on Fuhrmann et al. 2002
Implemented by Srikanth Ramaswamy, Blue Brain Project, July 2009
GUY: Removed  plasticity and depression

ENDCOMMENT


NEURON {

        POINT_PROCESS AMPA
        RANGE  tau_r_AMPA, tau_d_AMPA
        RANGE Use
        RANGE i,  i_AMPA,  g_AMPA, e, gmax
        NONSPECIFIC_CURRENT i
}

PARAMETER {

	   tau_r_AMPA = 0.3   (ms) : dual-exponential conductance profile
        tau_d_AMPA = 1.8     (ms) : IMPORTANT: tau_r < tau_d
        Use = 1.0   (1)   : Utilization of synaptic efficacy (just initial values! Use, Dep and Fac are overwritten by BlueBuilder assigned values)

        e = 0     (mV)  : AMPA and NMDA reversal potential
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
	i_AMPA (nA)
	g_AMPA (uS)
	factor_AMPA

}

STATE {


	A_AMPA       : AMPA state variable to construct the dual-exponential profile - decays with conductance tau_r_AMPA
    B_AMPA       : AMPA state variable to construct the dual-exponential profile - decays with conductance tau_d_AMPA
}

INITIAL{

    LOCAL  tp_AMPA

	A_AMPA = 0
	B_AMPA = 0

	tp_AMPA = (tau_r_AMPA*tau_d_AMPA)/(tau_d_AMPA-tau_r_AMPA)*log(tau_d_AMPA/tau_r_AMPA) :time to peak of the conductance



	factor_AMPA = -exp(-tp_AMPA/tau_r_AMPA)+exp(-tp_AMPA/tau_d_AMPA) :AMPA Normalization factor - so that when t = tp_AMPA, gsyn = gpeak
    factor_AMPA = 1/factor_AMPA

}

BREAKPOINT {

    SOLVE state METHOD cnexp
    g_AMPA = (B_AMPA-A_AMPA) :compute time varying conductance as the difference of state variables B_AMPA and A_AMPA and mggate kinetics

	i_AMPA = g_AMPA*(v-e) :compute the AMPA driving force based on the time varying conductance, membrane potential, and AMPA reversal
	i =  i_AMPA
}

DERIVATIVE state{


	A_AMPA' = -A_AMPA/tau_r_AMPA
    B_AMPA' = -B_AMPA/tau_d_AMPA
}





NET_RECEIVE (weight, weight_AMPA){

	weight_AMPA = weight



	A_AMPA = A_AMPA + weight_AMPA*factor_AMPA
    B_AMPA = B_AMPA + weight_AMPA*factor_AMPA


}
