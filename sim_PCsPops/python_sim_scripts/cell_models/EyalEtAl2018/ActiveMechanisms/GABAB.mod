TITLE GABAB receptor


COMMENT
GABAB receptor conductance using a dual-exponential profile
presynaptic short-term plasticity based on Fuhrmann et al, 2002
Implemented by Beatriz Herrera April 8 based on the implementation by
Srikanth Ramaswamy, Blue Brain Project, March 2009 and Guy 2018.
Some modificaations of Ramaswamy implementation were made based on Yang et al.
2016 Nat Commun. 7, 1–14 (2016).
GUY: Removed  plasticity and depression

ENDCOMMENT


NEURON {

        POINT_PROCESS GABAB
        RANGE  tau_r_GABAB, tau_d_GABAB
        RANGE Use
        RANGE i,  i_GABAB,  g_GABAB, e, gmax
        NONSPECIFIC_CURRENT i
}

PARAMETER {

	   tau_r_GABAB = 30   (ms) : dual-exponential conductance profile
        tau_d_GABAB = 80     (ms) : IMPORTANT: tau_r < tau_d
        Use = 1.0   (1)   : Utilization of synaptic efficacy (just initial values! Use, Dep and Fac are overwritten by BlueBuilder assigned values)

        e = 0     (mV)  : GABAB and NMDA reversal potential
    	:gmax = .001 (uS) :1nS weight conversion factor (from nS to uS)
    	u0 = 0 :initial value of u, which is the running value of Use
        fGABAB
}

COMMENT
The Verbatim block is needed to generate random nos. from a uniform distribution between 0 and 1
for comparison with Pr to decide whether to activate the synapse or not
ENDCOMMENT




ASSIGNED {

        v (mV)
        i (nA)
	i_GABAB (nA)
	g_GABAB (uS)
	factor_GABAB
    factor_GABAB1

}

STATE {


	A_GABAB       : GABAB state variable to construct the dual-exponential profile - decays with conductance tau_r_GABAB
    B_GABAB       : GABAB state variable to construct the dual-exponential profile - decays with conductance tau_d_GABAB
}

INITIAL{

    LOCAL  tp_GABAB, temp_GABAB



	A_GABAB = 0
	B_GABAB = 0

	temp_GABAB = tau_r_GABAB/tau_d_GABAB :time to peak of the conductance

	factor_GABAB1 = temp_GABAB^(1/(temp_GABAB-1)) :GABAB Normalization factor - such that the peak of g_GABAB is 1

    tp_GABAB = (tau_r_GABAB*tau_d_GABAB)/(tau_d_GABAB-tau_r_GABAB)*log(tau_d_GABAB/tau_r_GABAB) :time to peak of the conductance

	factor_GABAB = -exp(-tp_GABAB/tau_r_GABAB)+exp(-tp_GABAB/tau_d_GABAB) :GABAB Normalization factor - so that when t = tp_GABAB, gsyn = gpeak
    factor_GABAB = 1/factor_GABAB

}

BREAKPOINT {

    SOLVE state METHOD cnexp
    fGABAB = 33.33*(0.5-2/(1+exp((v + 98.73)/12.5))) :fGABAB kinetics - Shoemaker, Neurocomputing 2011
    g_GABAB = B_GABAB :time varying conductance is given by the state variables B_GABAB

	i_GABAB = g_GABAB*fGABAB :compute the GABAB driving force based on the time varying conductance, and fGABAB(V) kinetics
	i =  i_GABAB
}

DERIVATIVE state{


	A_GABAB' = -A_GABAB/tau_r_GABAB
    B_GABAB' = (-B_GABAB + factor_GABAB1*A_GABAB)/tau_d_GABAB
}





NET_RECEIVE (weight, weight_GABAB){

	weight_GABAB = weight



	A_GABAB = A_GABAB + weight_GABAB
    B_GABAB = B_GABAB + weight_GABAB


}
