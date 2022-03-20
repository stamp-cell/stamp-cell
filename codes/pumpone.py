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
NumberofEquations = 10
NumberofICs = 2
NumberofParameters = 10

#################### Add something like if the answer true then continue
#####------ Initial Conditions:
AmInput = input('Enter [A]_M concentration (mM) ') # A internal concentration (mM) ')
AnInput = input('Enter [A]_N concentration (mM) ') #A external concentration (mM) ')
#BmInput = input('Enter [B]_M concentration (mM) ') #B internal concentration (mM) ')
#BnInput = input('Enter [B]_N concentration (mM) ') #B external concentration (mM) ')

#####------ Parameters:
K_EAPm = input('Enter K_EAPm (mM) ') #('Enter A internal dissocication constant (mM) ')
K_EAPn = input('Enter K_EAPn (mM) ') #('Enter A external dissocication constant (mM) ')
K_A = input('Enter K_A (mM) ') #('Enter B external dissocication constant (mM) ')
g = input('Enter g (per_s) , (EAP translocation rate constant from ext to int) ')

#g_EB_N = input('Enter g_EB_N (per_s) , (EB translocation rate constant from ext to int) ')
#g_EB_M = input('Enter g_EB_M (per_s) , (EB translocation rate constant from int to ext) ')
#g_EA_N = input('Enter g_EA_N (per_s) , (EA translocation rate constant from ext to int) ')
#g_EA_M = input('Enter g_EA_M (per_s) , (EA translocation rate constant from int to ext) ')


def createLegends():
    legend_ICs = [""] * NumberofICs
    legend_equation = [""] * NumberofEquations
    legend_voi = ""
    legend_parameters = [""] * NumberofParameters
    legend_voi = "time in component environment (second)"
    legend_ICs[0] = "HCO3_int in component concentrations (mM)"
    legend_ICs[1] = "HCO3_ext in component concentrations (mM)"
    legend_equation[8] = "J_AE1_HCO3 in component AE1 (mM_per_s)"
    legend_parameters[1] = "K_B_M (mM)" # "K_HCO3_int in component AE1 (mM)"
    legend_parameters[0] = "K_A_M (mM)" # "K_Cl_int in component AE1 (mM)"
    legend_parameters[2] = "K_A (mM)"
    legend_parameters[3] = "g_EB_M (per_s)" # "P_HCO3_int in component AE1 (per_s)"
    legend_parameters[4] = "g_EA_N (per_s)" # "P_Cl_ext in component AE1 (per_s)"
    legend_parameters[5] = "g_EA_M (per_s)" # "P_Cl_int in component AE1 (per_s)"
    legend_equation[1] = "alpha_M" # "gamma_int in component AE1 (dimensionless)"
    legend_equation[2] = "beta_M" # "beta_int in component AE1 (dimensionless)"
    legend_equation[3] = "alpha_N" # "gamma_ext in component AE1 (dimensionless)"
    legend_equation[4] = "beta_N" # "beta_ext in component AE1 (dimensionless)"
    legend_equation[7] = "sigma in component AE1 (per_s)"
    legend_parameters[4]  = "E_Tmax in component AE1 (mM)"
    legend_parameters[5]   = "K_I in component AE1 (mM)"
    legend_equation[0] = "E_T (mM)"
    return (legend_ICs, legend_equation, legend_voi, legend_parameters)

#_____________________________________________________________________
def ic_parameters():
	parameters = [0.0] * NumberofParameters; ICs = [0.0] * NumberofICs;
	ICs[0] = float(AmInput) # A_M 	# (50.0) A_int concentrations (mM)
	ICs[1] = float(AnInput) # A_N 	# (0.0) A_ext  concentrations (mM)
	#ICs[2] = float(BmInput) # B_M	# (0.0) B_int  concentrations (mM)
	#ICs[3] = float(BnInput) # B_N	# (0.0) B_ext  concentrations (mM)
# Model Parameters:
	parameters[0] = float(K_EAPm) # K_A_ATP	#50 K"_c, K_Cl_int in component AE1 (per_s)
	parameters[1] = float(K_EAPn) # K_A_ATP	#50 K"_c, K_Cl_int in component AE1 (per_s)
	parameters[2] = float(K_A)  # K_A	# 198 K'_b, K_HCO3_ext in component AE1 (per_s)
	parameters[3] = float(g)  
	parameters[4] = 1	#E_Tmax (mM)
	parameters[5] = 172	#K_I  (mM)
	parameters[6] = 0.00000
	parameters[7] = 1.00000
	parameters[8] = 0.00000
	parameters[9] = 0.00000
	return (ICs, parameters)

def computeRates(voi, ICs, parameters):
    rates = [0.0] * NumberofICs; equation = [0.0] * NumberofEquations
    rates[0] = parameters[6]   #d/dt HCO3_int in component concentrations (mM)
    rates[1] = parameters[7]   #d/dt HCO3_ext in component concentrations (mM)
    #rates[2] = parameters[8]   #d/dt Cl_int in component concentrations (mM)
    #rates[3] = parameters[9]   #d/dt Cl_ext in component concentrations (mM)
    return(rates)

def computeEquation(parameters, ICs, voi):
	equation = array([[0.0] * len(voi)] * NumberofEquations)
	ICs = array(ICs)
	voi = array(voi)
# 0) E_T in component AE1 (mM):
	equation[0] = parameters[4]/(1.00000+ICs[0]/parameters[5])
# 1) J_Denom_First = (K_EAPm + K_EAPn + 2*K_EAPm*K_EAPn)*A_M*A_N  
	equation[1] = (parameters[0] + parameters[1] + 2*parameters[0]*parameters[1])*ICs[0]*ICs[1]
# 2) J_Denom_Second = K_A*(1+2*K_EAPm)*A_M+2*K_A^2
	equation[2] = parameters[2]*(1.00000 + 2.00000 * parameters[0])*ICs[0]
# 3) J_Denom_Third = K_A*(1+2*K_EAPn)*A_N
	equation[3] = parameters[2]*(1.00000 + 2.00000 * parameters[1])*ICs[1]
# 4) J_Denom_Fourth = 2*K_A^2
	equation[4] = (2.0000*parameters[2]**2)
# 5) J Denominator: 
	equation[5] = (equation[1] + equation[2] + equation[3] + equation[4])
# 6) J_A First_Numerator: = (K_EAPm*A_M)
	equation[6] = (parameters[0]*ICs[0])
# 7) J_A Second_Numerator: = (K_EAPn*A_N)
	equation[7] = (parameters[1]*ICs[1])
# 8) J_A Numerator:
	equation[8] = parameters[2]*parameters[3]*(equation[6]-equation[7])
# 9) J_A:
	equation[9] = equation[8]/equation[5]
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
print("J_A = " , equation)


