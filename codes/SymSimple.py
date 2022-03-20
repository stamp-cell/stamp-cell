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
NumberofEquations = 18
NumberofICs = 6
NumberofParameters = 14

#################### Modify for g_ and ^
#####------ Initial Conditions:
AmInput = input('Enter [A]_M (mM) (A internal concentration) ')
AnInput = input('Enter [A]_N (mM) (A external concentration) ')
BmInput = input('Enter [B]_M (mM) (B internal concentration) ')
BnInput = input('Enter [B]_N (mM) (B external concentration) ')
CmInput = input('Enter [C]_M (mM) (C internal concentration) ')
CnInput = input('Enter [C]_N (mM) (C external concentration) ')

#####------ Parameters:
K_A = input('Enter K_A (mM) ')  #A dissocication constant (mM) 
K_B = input('Enter K_B (mM) ') #B dissocication constant (mM) 
K_C = input('Enter K_C (mM) ') #C dissocication constant (mM) 
a = input('Enter A stoichiometry ratio ') 
b = input('Enter B stoichiometry ratio ') 
c = input('Enter C stoichiometry ratio ')  


g_sym = input('Enter g_symporter (per_s) , (E translocation rate constant from ext to int) ')
E_t = input(' Enter E_t ')
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
    legend_parameters[7] = "None (mM)"
    legend_parameters[1] = "None (mM)"
    legend_parameters[6] = "None (mM)" 
    legend_parameters[0] = "None (mM)" 
    legend_parameters[8] = "None (per_s)" 
    legend_parameters[11] = "None (per_s)"
    legend_parameters[7] = "None (per_s)" 
    legend_parameters[12] = "None (per_s)"
    legend_equation[1] = "None" 
    legend_equation[3] = "None" 
    legend_equation[4] = "None"
    legend_equation[6] = "None"
    legend_equation[13] = "None"
    legend_parameters[13] = "E_Tmax in component AE1 (mM)"
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
	parameters[0] = float(K_A) 
	parameters[1] = float(K_B) 
	parameters[2] = float(K_C)
	parameters[3] = float(a)
	parameters[4] = float(b)	
	parameters[5] = float(c)

	parameters[6] = float(g_sym)
	parameters[7] = float(E_t)	#E_Tmax (mM)

	parameters[8] = 0.00000
	parameters[9] = 1.00000
	parameters[10] = 0.00000
	parameters[11] = 0.00000
	parameters[12] = 0.00000
	parameters[13] = 0.00000
	return (ICs, parameters)

def computeRates(voi, ICs, parameters):
    rates = [0.0] * NumberofICs; equation = [0.0] * NumberofEquations
    rates[0] = parameters[8] #d/dt A_M in component concentrations (mM) int=M
    rates[1] = parameters[9] #d/dt A_N in component concentrations (mM) ext=N
    rates[2] = parameters[10] #d/dt B_M in component concentrations (mM)
    rates[3] = parameters[11] #d/dt B_N in component concentrations (mM)
    rates[4] = parameters[12] #d/dt C_M in component concentrations (mM)
    rates[5] = parameters[13] #d/dt C_N in component concentrations (mM)
    return(rates)

def computeEquation(parameters, ICs, voi):
	equation = array([[0.0] * len(voi)] * NumberofEquations)
	ICs = array(ICs)
	voi = array(voi)
# 0) P_Symp: E_t*g_symp
	equation[0] = parameters[7]*parameters[6]
# 1) alpha_M (dimensionless)     
	equation[1] = ICs[0]/parameters[0]
# 2) beta_M (dimensionless)     
	equation[2]=  ICs[2]/parameters[1]
# 3) gamma_M (dimensionless)    
	equation[3] = ICs[4]/parameters[2]
# 4) alpha_N (dimensionless) #     
	equation[4] = ICs[1]/parameters[0]
# 5) beta_N (dimensionless)     
	equation[5]=  ICs[3]/parameters[1]
# 6) gamma_N (dimensionless)    
	equation[6] = ICs[5]/parameters[2]
# 7) J_Symp First_Denominator: (1+alpha_N)
	equation[7] = (1+equation[4])
# 8) J_Symp Second_Denominator: (1+beta_N)
	equation[8] = (1+equation[5])
# 9) J_Symp Third_Denominator: (1+gamma_N)
	equation[9] = (1+equation[6])
# 10) J_Symp Denominator: ((1+alpha_N)^a)*((1+beta_N)^b)*((1+gamma_N)^c)
	equation[10] = ((power(equation[7], parameters[3]))*(power(equation[8], parameters[4]))*(power(equation[9], parameters[5])))
# 11) J_symp First_Numerator: = A_M^a * B_M^b * C_M*c
	equation[11] = ((power(ICs[0], parameters[3]))*(power(ICs[2], parameters[4]))*(power(ICs[4], parameters[5])))
# 12) J_symp Second_Numerator: = A_N^a * B_N^b * C_N*c
	equation[12] = ((power(ICs[1], parameters[3]))*(power(ICs[3], parameters[4]))*(power(ICs[5], parameters[5])))
# 13) J_symp Numerator:
	equation[13] = (equation[11]-equation[12])
# 14) J_symp(M, N):
	equation[14] = equation[0] * (equation[13]/equation[10])
# 15) J_A(M, N): = a*J_symp
	equation[15] = parameters[3]*equation[14]
# 16) J_B(M, N): b*J_symp
	equation[16] = parameters[4]*equation[14]
# 17) J_C(M, N): = c*J_symp
	equation[17] = parameters[5]*equation[14]
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
