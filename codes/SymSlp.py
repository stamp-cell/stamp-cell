from __future__ import division
import numpy as np
import math as m
import pylab
import pandas as pd
import scipy 
import scipy.stats
from math import *
from numpy import *
from scipy.integrate import ode
NumberofEquations = 20
NumberofICs = 4
NumberofParameters = 20

#################### Add something like if the answer true then continue
#####------ Initial Conditions:
AmInput = input('Enter A internal concentration (mM) ')
AnInput = input('Enter A external concentration (mM) ')
BmInput = input('Enter B internal concentration (mM) ')
BnInput = input('Enter B external concentration (mM) ')

#####------ Parameters:
K_Am = input('Enter K_Am (mM) ')  #A internal dissocication constant (mM) 
K_Bm = input('Enter K_Bm (mM) ') #B internal dissocication constant (mM) 
K_BAm = input('Enter K_BAm (mM) ') #A internal dissocication constant (mM) 
K_An = input('Enter K_An (mM) ') #A external dissocication constant (mM) 
K_Bn = input('Enter K_Bn (mM) ') #B external dissocication constant (mM) 
K_BAn = input('Enter K_BAn (mM) ') #B external dissocication constant (mM) 


g_E_N = input(' Enter g_E_N (per_s) , (E translocation rate constant from ext to int) ')
g_EA_N = input(' Enter g_EA_N (per_s) , (EA translocation rate constant from ext to int) ')
g_EB_N = input(' Enter g_EB_N (per_s) , (EB translocation rate constant from ext to int) ')
g_EAB_N = input(' Enter g_EAB_N (per_s) , (EAB translocation rate constant from ext to int) ')
g_E_M = input(' Enter g_E_M (per_s) , (E translocation rate constant from int to ext) ')
g_EB_M = input(' Enter g_EB_M (per_s) , (EB translocation rate constant from int to ext) ')
g_EA_M = input(' Enter g_EA_M (per_s) , (EA translocation rate constant from int to ext) ')
g_EAB_M = input(' Enter g_EAB_M (per_s) , (EAB translocation rate constant from int to ext) ')

#M side: 
def terms(Am,k_Am,Bm,k_Bm):
	alpha_m = Am/k_Am
	beta_m = Bm/k_Bm
	return(alpha_m,beta_m)
print terms(3,4,5,6)


def createLegends():
    legend_ICs = [""] * NumberofICs
    legend_equation = [""] * NumberofEquations
    legend_voi = ""
    legend_parameters = [""] * NumberofParameters
    legend_voi = "time in component environment (second)"
    legend_ICs[0] = "HCO3_int in component concentrations (mM)"
    legend_ICs[1] = "HCO3_ext in component concentrations (mM)"
    legend_ICs[2] = "Cl_int in component concentrations (mM)"
    legend_ICs[3] = "Cl_ext in component concentrations (mM)"
    legend_equation[14] = "J_AE1_HCO3 in component AE1 (mM_per_s)"
    legend_equation[15] = "J_HCO3_influx in component AE1 (mM_per_s)"
    legend_equation[19] = "J_Cl_influx in component AE1 (mM_per_s)"
    legend_parameters[4] = "K_B_N (mM)" # "K_HCO3_ext in component AE1 (mM)"
    legend_parameters[1] = "K_B_M (mM)" # "K_HCO3_int in component AE1 (mM)"
    legend_parameters[3] = "K_A_N (mM)" # "K_Cl_ext in component AE1 (mM)"
    legend_parameters[0] = "K_A_M (mM)" # "K_Cl_int in component AE1 (mM)"
    legend_parameters[8] = "g_EB_N (per_s)" # "P_HCO3_ext in component AE1 (per_s)"
    legend_parameters[11] = "g_EB_M (per_s)" # "P_HCO3_int in component AE1 (per_s)"
    legend_parameters[7] = "g_EA_N (per_s)" # "P_Cl_ext in component AE1 (per_s)"
    legend_parameters[12] = "g_EA_M (per_s)" # "P_Cl_int in component AE1 (per_s)"
    legend_equation[1] = "alpha_M" # "gamma_int in component AE1 (dimensionless)"
    legend_equation[3] = "beta_M" # "beta_int in component AE1 (dimensionless)"
    legend_equation[4] = "alpha_N" # "gamma_ext in component AE1 (dimensionless)"
    legend_equation[6] = "beta_N" # "beta_ext in component AE1 (dimensionless)"
    legend_equation[13] = "sigma in component AE1 (per_s)"
    legend_parameters[14] = "E_Tmax in component AE1 (mM)"
    legend_parameters[15] = "K_I in component AE1 (mM)"
    legend_equation[0] = "E_T (mM)"
    return (legend_ICs, legend_equation, legend_voi, legend_parameters)

#_____________________________________________________________________
def ic_parameters():
	parameters = [0.0] * NumberofParameters; ICs = [0.0] * NumberofICs;
	ICs[0] = float(AmInput) 	# A_M 	# int concentrations (mM)
	ICs[1] = float(AnInput)  # A_N 	# ext concentrations (mM)
	ICs[2] = float(BmInput) # B_M	# int concentrations (mM)
	ICs[3] = float(BnInput) # B_N	# ext concentrations (mM)
# Model Parameters:
	parameters[0] = float(K_Am) 
	parameters[1] = float(K_Bm) 
	parameters[2] = float(K_BAm)
	parameters[3] = float(K_An) 
	parameters[4] = float(K_Bn) 
	parameters[5] = float(K_BAn) 
	parameters[6] = float(g_E_N) 
	parameters[7] = float(g_EA_N)
	parameters[8] = float(g_EB_N) 
	parameters[9] = float(g_EAB_N)
	parameters[10] = float(g_E_M) 
	parameters[11] = float(g_EB_M)
	parameters[12] = float(g_EA_M)
	parameters[13] = float(g_EAB_M)
	parameters[14] = 1	#E_Tmax (mM)
	parameters[15] = 172	#K_I  (mM)
	parameters[16] = 0.00000
	parameters[17] = 1.00000
	parameters[18] = 0.00000
	parameters[19] = 0.00000
	return (ICs, parameters)

def computeRates(voi, ICs, parameters):
    rates = [0.0] * NumberofICs; equation = [0.0] * NumberofEquations
    rates[0] = parameters[16] #d/dt HCO3_int in component concentrations (mM)
    rates[1] = parameters[17] #d/dt HCO3_ext in component concentrations (mM)
    rates[2] = parameters[18] #d/dt Cl_int in component concentrations (mM)
    rates[3] = parameters[19] #d/dt Cl_ext in component concentrations (mM)
    return(rates)

def computeEquation(parameters, ICs, voi):
	equation = array([[0.0] * len(voi)] * NumberofEquations)
	ICs = array(ICs)
	voi = array(voi)
# 0) E_T in component AE1 (mM):
	equation[0] = parameters[14]/(1.00000+ICs[0]/parameters[15])
# 1) alpha_M (dimensionless) # gamma_int in component AE1 (dimensionless)    
	equation[1] = ICs[0]/parameters[0]
# 2) alphaprime_M (dimensionless)    
	equation[2] = ICs[0]/parameters[2]
# 3) beta_M (dimensionless) # beta_int in component AE1 (dimensionless):    
	equation[3]=  ICs[2]/parameters[1]
# 4) alpha_N (dimensionless) # gamma_ext in component AE1 (dimensionless):    
	equation[4] = ICs[1]/parameters[3]
# 5) alphaprime_n (dimensionless)    
	equation[5] = ICs[1]/parameters[5]
# 6) beta_N (dimensionless)# beta_ext in component AE1 (dimensionless):
	equation[6]= ICs[3]/parameters[5]
# 7) R_N = 1 + alpha_n + beta_n + alphaprime_n*beta_n
	equation[7] = (1.00000 + equation[4] + equation[6] + equation[5]*equation[6])
# 8) R_MM = g_E_M + g_EA_M*alpha_m + g_EB_M*beta_M + g_EAB_M*alpha_prime_M*beta_M
	equation[8] = (parameters[10] + parameters[12]*equation[1] + parameters[11]*equation[3] + parameters[13]*equation[2]*equation[3])
# 9) R_M = 1 + alpha_m + beta_m + alphaPrime_M*beta_M
	equation[9] = (1.00000 + equation[1] + equation[3] + equation[2]*equation[3])
# 10) R_NN = g_E_N + g_EA_N*alpha_N + g_EB_N*beta_N + g_EAB_N*alphaprime_N*beta_N
	equation[10] = (parameters[6] + parameters[7]*equation[4] + parameters[8]*equation[6] + parameters[9]*equation[5]*equation[6])
# 11) J Denominator: R_M*R_NN + R_N*R_MM
	equation[11] = (equation[7]*equation[8])+(equation[9]*equation[10])
# 12) J_A First_Numerator: = (alpha_M*g_EA_M+g_EAB_M*beta_M*alphaprime_M)(g_E_N+g_EB_N*beta_N)
	equation[12] = (equation[1]*parameters[12]+parameters[13]*equation[3]*equation[2])*(parameters[6]+parameters[8]*equation[6])
# 13) J_A Second_Numerator: = (alpha_N*g_EA_N+g_EAB_N*beta_N*alphaprime_N)(g_E_M+g_EB_M*beta_M)
	equation[13] = (equation[4]*parameters[7]+parameters[9]*equation[6]*equation[5])*(parameters[10]+parameters[11]*equation[3])
# 14) J_A Numerator:
	equation[14] = equation[12]-equation[13]
# 15) J_A(M, N):
	equation[15] = equation[0] * (equation[14]/equation[11])
# 16) J_B First_Numerator: = (beta_M*g_EB_M+g_EAB_M*beta_M*alphaprime_M)(g_E_N+g_EA_N*alpha_N)
	equation[16] = (equation[3]*parameters[11]+parameters[13]*equation[3]*equation[2])*(parameters[6]+parameters[12]*equation[4])
# 17) J_B Second_Numerator: = (beta_N*g_EB_N+g_EAB_N*beta_N*alphaprime_N)(g_E_M+g_EA_M*alpha_M)
	equation[17] = (equation[6]*parameters[8]+parameters[9]*equation[6]*equation[5])*(parameters[10]+parameters[12]*equation[1])
# 18) J_B Numerator:
	equation[18] = equation[16]-equation[17]
# 19) J_B(M, N):
	equation[19] = equation[0] * (equation[18]/equation[11])
	return equation
#_____________________________________________________________________
def solve_model():
    """Solve model with ODE solver"""
    # Initialise parameters and state variables
    (init_ICs, parameters) = ic_parameters()

    # Set timespan to solve over
    voi = linspace(0, 10, 500)

    # Construct ODE object to solve
    r = ode(computeRates)
    r.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=1)
    r.set_initial_value(init_ICs, voi[0])
    r.set_f_params(parameters)

    # Solve model
    ICs = array([[0.0] * len(voi)] * NumberofICs)
    ICs[:,0] = init_ICs
    for (i,t) in enumerate(voi[1:]):
        if r.successful():
            r.integrate(t)
            ICs[:,i+1] = r.y
        else:
            break

    # Compute equation variables
    equation = computeEquation(parameters, ICs, voi)
    return (voi, ICs, equation)

def plot_model(voi, ICs, equation):
    """Plot variables against variable of integration"""
    (legend_ICs, legend_equation, legend_voi, legend_parameters) = createLegends()
    pylab.figure(1)
    pylab.plot(voi,vstack((ICs,equation)).T)
    pylab.xlabel(legend_voi)
    pylab.legend(legend_ICs + legend_equation, loc='best')
    pylab.show()

if __name__ == "__main__":
    (voi, ICs, equation) = solve_model()
    plot_model(voi, ICs, equation)

###### output as csv file or table as well
print(equation)


