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
NumberofEquations = 26
NumberofICs = 6
NumberofParameters = 26

#################### Add something like if the answer true then continue
#####------ Initial Conditions:
AmInput = input('Enter [A]_M (mM) (A internal concentration) ')
AnInput = input('Enter [A]_N (mM) (A external concentration) ')
BmInput = input('Enter [B]_M (mM) (B internal concentration) ')
BnInput = input('Enter [B]_N (mM) (B external concentration) ')
CmInput = input('Enter [C]_M (mM) (C internal concentration) ')
CnInput = input('Enter [C]_N (mM) (C external concentration) ')

#####------ Parameters:
K_Am = input('Enter K_Am (mM) ')  #A internal dissocication constant (mM) 
K_ABm = input('Enter K_ABm (mM) ') #B internal dissocication constant (mM) 
K_ACm = input('Enter K_ACm (mM) ') #C internal dissocication constant (mM) 
K_An = input('Enter K_An (mM) ') #A external dissocication constant (mM) 
K_ABn = input('Enter K_ABn (mM) ') #B external dissocication constant (mM) 
K_ACn = input('Enter K_ACn (mM) ') #C external dissocication constant (mM) 


g_E_N = input('Enter g_E_N (per_s) , (E translocation rate constant from ext to int) ')
g_EAC_N = input('Enter g_EAC_N (per_s) , (EAC translocation rate constant from ext to int) ')
g_EAB_N = input('Enter g_EAB_N (per_s) , (EAB translocation rate constant from ext to int) ')

g_E_M = input(' Enter g_E_M (per_s) , (E translocation rate constant from int to ext) ')
g_EAC_M = input(' Enter g_EAC_M (per_s) , (EAB translocation rate constant from int to ext) ')
g_EAB_M = input(' Enter g_EAB_M (per_s) , (EAC translocation rate constant from int to ext) ')
#M side: 
#def terms(Am,k_Am,Bm,k_Bm):
#	alpha_m = Am/k_Am
#	beta_m = Bm/k_Bm
#	return(alpha_m,beta_m)
#print terms(3,4,5,6)


def createLegends():
    legend_ICs = [""] * NumberofICs
    legend_equation = [""] * NumberofEquations
    legend_voi = ""
    legend_parameters = [""] * NumberofParameters
    legend_voi = "time in component environment (second)"
    legend_ICs[0] = "None_int in component concentrations (mM)"
    legend_ICs[1] = "None_ext in component concentrations (mM)"
    legend_ICs[2] = "None_int in component concentrations (mM)"
    legend_ICs[3] = "None_ext in component concentrations (mM)"
    legend_equation[14] = "None_HCO3 in component AE1 (mM_per_s)"
    legend_equation[15] = "None_influx in component AE1 (mM_per_s)"
    legend_equation[19] = "None_influx in component AE1 (mM_per_s)"
    legend_parameters[3] = "None (mM)"
    legend_parameters[1] = "None (mM)"
    legend_parameters[3] = "None (mM)" 
    legend_parameters[0] = "None (mM)" 
    legend_parameters[7] = "None (per_s)" 
    legend_parameters[10] = "None (per_s)"
    legend_parameters[6] = "None (per_s)" 
    legend_parameters[11] = "None (per_s)"
    legend_equation[1] = "None" 
    legend_equation[3] = "None" 
    legend_equation[4] = "None"
    legend_equation[6] = "None"
    legend_equation[13] = "None"
    legend_parameters[12] = "E_Tmax in component AE1 (mM)"
    legend_parameters[13] = "K_I in component AE1 (mM)"
    legend_equation[0] = "E_T (mM)"
    return (legend_ICs, legend_equation, legend_voi, legend_parameters)

#_____________________________________________________________________
def ic_parameters():
	parameters = [0.0] * NumberofParameters; ICs = [0.0] * NumberofICs;
	ICs[0] = float(AmInput) # A_M 	# int concentrations (mM)
	ICs[1] = float(AnInput)  # A_N 	# ext concentrations (mM)
	ICs[2] = float(BmInput) # B_M	# int concentrations (mM)
	ICs[3] = float(BnInput) # B_N	# ext concentrations (mM)
	ICs[4] = float(CmInput) # C_M	# int concentrations (mM)
	ICs[5] = float(CnInput) # C_N	# ext concentrations (mM)
# Model Parameters:
	parameters[0] = float(K_Am) 
	parameters[1] = float(K_ABm) 
	parameters[2] = float(K_ACm)
	parameters[3] = float(K_An) 
	parameters[4] = float(K_ABn) 
	parameters[5] = float(K_ACn)
	parameters[6] = float(g_E_N) 
	parameters[7] = float(g_EAC_N) 
	parameters[8] = float(g_EAB_N)
	parameters[9] = float(g_E_M ) 
	parameters[10] = float(g_EAC_M)
	parameters[11] = float(g_EAB_M)
	parameters[12] = 1	#E_Tmax (mM)
	parameters[13] = 172	#K_I  (mM)
	parameters[14] = 0.00000
	parameters[15] = 1.00000
	parameters[16] = 0.00000
	parameters[17] = 0.00000
	parameters[18] = 0.00000
	parameters[19] = 0.00000
	return (ICs, parameters)

def computeRates(voi, ICs, parameters):
    rates = [0.0] * NumberofICs; equation = [0.0] * NumberofEquations
    rates[0] = parameters[14] #d/dt A_M in component concentrations (mM) int=M
    rates[1] = parameters[15] #d/dt A_N in component concentrations (mM) ext=N
    rates[2] = parameters[16] #d/dt B_M in component concentrations (mM)
    rates[3] = parameters[17] #d/dt B_N in component concentrations (mM)
    rates[4] = parameters[18] #d/dt C_M in component concentrations (mM)
    rates[5] = parameters[19] #d/dt C_N in component concentrations (mM)
    return(rates)

def computeEquation(parameters, ICs, voi):
	equation = array([[0.0] * len(voi)] * NumberofEquations)
	ICs = array(ICs)
	voi = array(voi)
# 0) E_T in component AE1 (mM):
	equation[0] = parameters[12]/(1.00000+ICs[0]/parameters[13])
# 1) alpha_M (dimensionless)     
	equation[1] = ICs[0]/parameters[0]
# 2) beta_M (dimensionless)     
	equation[2]=  ICs[2]/parameters[1]
# 3) gamma_M (dimensionless)    
	equation[3] = ICs[4]/parameters[2]
# 4) alpha_N (dimensionless) #     
	equation[4] = ICs[1]/parameters[3]
# 5) beta_N (dimensionless)     
	equation[5]=  ICs[3]/parameters[3]
# 6) gamma_N (dimensionless)    
	equation[6] = ICs[5]/parameters[4]
# 7) R_N = 1 + alpha_N + alpha_N*beta_N + alpha_N*gamma_N
	equation[7] = (1.0 + equation[4] + equation[4]*equation[5] + equation[4]*equation[6])
# 8) R_MM = g_E_M + g_EAB_M*alpha_M*beta_M + g_EAC_M*alpha_M*gamma_M
	equation[8] = (parameters[9] + parameters[11]*equation[1]*equation[2] + parameters[10]*equation[1]*equation[3])
# 9) R_M = 1 + alpha_M + alpha_M*beta_M + alpha_M*gamma_M
	equation[9] = (1.0 + equation[1] + equation[1]*equation[2] + equation[1]*equation[3])
# 10) R_NN = g_E_N + g_EAB_N*alpha_N*beta_N + g_EAC_N*alpha_N*gamma_N
	equation[10] = (parameters[6] + parameters[8]*equation[4]*equation[5] + parameters[7]*equation[4]*equation[6])
# 11) J Denominator: R_M*R_NN + R_N*R_MM
	equation[11] = (equation[9]*equation[10]+equation[7]*equation[8])
# 12) J_B First_Numerator: = R_NN*(g_EAB_M*alpha_M*beta_M)
	equation[12] = equation[10]*(parameters[8]*equation[1]*equation[2])
# 13) J_B Second_Numerator: = R_MM*(g_EAB_N*alpha_N*beta_N)
	equation[13] = equation[8]*(parameters[8]*equation[4]*equation[5])
# 14) J_B Numerator:
	equation[14] = (equation[12]-equation[13])
# 15) J_B(M, N):
	equation[15] = equation[0] * (equation[14]/equation[11])
# 16) J_C First_Numerator: = R_NN*(g_EAC_M*alpha_M*gamma_M)
	equation[16] = equation[10]*(parameters[10]*equation[1]*equation[3])
# 17) J_C Second_Numerator: = R_MM*(g_EAC_N*alpha_N*gamma_N)
	equation[17] = equation[8]*(parameters[7]*equation[4]*equation[6])
# 18) J_C Numerator:
	equation[18] = (equation[16]-equation[17])
# 19) J_C(M, N):
	equation[19] = equation[0] * (equation[18]/equation[11])
# 20) J_A(M, N): = J_B(M, N) + J_C(M, N) 
	equation[20] = (equation[15]+equation[19])
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


