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


# Size of variable arrays:
NumberofEquations = 2
NumberofICs = 1
NumberofParameters = 9


#####------ Initial Conditions:
AmInput = input('Enter [A]_M concentration (mM) ') # A internal concentration (mM) ')

#####------ Parameters:
AnInput = input('Enter [A]_N concentration (mM) ') #A external concentration (mM) ')
BmInput = input('Enter [B]_M (mM) (B internal concentration) ')
BnInput = input('Enter [B]_N (mM) (B external concentration) ')
J_A_Max = input('Enter J_A, pump_Max  (mM_per_s)')

a = input('Enter A stoichiometry ratio ') 
b = input('Enter B stoichiometry ratio ') 
#c = input('Enter C stoichiometry ratio ')  


def createLegends():
	legend_ICs = [""] * NumberofICs
	legend_equation = [""] * NumberofEquations
	legend_voi = ""
	legend_parameters = [""] * NumberofParameters
	legend_voi = "time in component environment (second)"
	legend_ICs[0] = "Na_int in component concentrations (mM)"
	legend_parameters[0] = "Na_ext in component concentrations (mM)"
	legend_parameters[1] = "K_int in component concentrations (mM)"
	legend_parameters[2] = "K_ext in component concentrations (mM)"
	legend_equation[0] = "J_NaK_Na in component NaK (mM_per_s)"
	legend_equation[1] = "J_NaK_K in component NaK (mM_per_s)"
	legend_parameters[4]  = "K_Na in component NaK (mM)"
	legend_parameters[5]  = "K_K in component NaK (mM)"
	legend_parameters[3]  = "J_NaK_Na_Max in component NaK (mM_per_s)"
	#legend_rates[0] = "d/dt Na_int in component concentrations (mM)"
	return (legend_ICs, legend_equation, legend_voi, legend_parameters)

def initConsts():
	parameters = [0.0] * NumberofParameters; ICs = [0.0] * NumberofICs;
	ICs[0] = float(AmInput) #0.0
	parameters[0] = float(AnInput) # 0.0 Na_ext in component concentrations (mM)
	parameters[1] = float(BmInput) # 39.0 K_int in component concentrations (mM)
	parameters[2] = float(BnInput) # 10.0 K_ext in component concentrations (mM)
	parameters[3] = float(J_A_Max) # 10.8 (mM_per_s)
	parameters[4] = 0.200000*(1.00000+parameters[1]/8.33000) #K_Na in component NaK (mM)
	parameters[5] = 0.100000*(1.00000+parameters[0]/18.5000) #K_K in component NaK (mM)
	parameters[6] = 1.00000 # d/dt Na_int in component concentrations (mM)
	parameters[7] = float(a)
	parameters[8] = float(b)	
	#parameters[T5] = float(c)

    #parameters[7] = 0.200000*parameters[5] 
	return (ICs, parameters)

def computeRates(voi, ICs, parameters):
	rates = [0.0] * NumberofICs; equation = [0.0] * NumberofEquations
	rates[0] = parameters[6] # d/dt Na_int in component concentrations (mM)
	return(rates)

def computeequation(parameters, ICs, voi):
	equation = array([[0.0] * len(voi)] * NumberofEquations)
	ICs = array(ICs)
	voi = array(voi)
	equation[0] = parameters[3] *(power(ICs[0]/(ICs[0]+parameters[4] ), parameters[7]))*(power(parameters[2]/(parameters[2]+parameters[5] ), parameters[8]))
	equation[1] = ((-parameters[8]/parameters[7])*equation[0])
	return equation

def solve_model():
	"""Solve model with ODE solver"""
	# Initialise parameters and state variables
	(init_ICs, parameters) = initConsts()
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
	equation = computeequation(parameters, ICs, voi)
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
