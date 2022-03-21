#  Parameters are changed based on :<<A Mathematical Model of the Proton Balance in the Outer Mantle Epithelium of Anodonta cygnea L>> Na-K

from __future__ import division
import os
import numpy as np
import math as m
import pylab
import pandas as pd
import csv
import scipy 
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
import json
import pytz
from pandas import Series, DataFrame
from dateutil.parser import parse # dateutil is capable of parsing almost any human-intelligible date representation
from math import *
from numpy import *
from scipy.integrate import ode
from collections import OrderedDict

print ('\n' + "Please enter number for following inputs:" + '\n') 
# Size of variable arrays:
NumberofEquations = 2
NumberofICs = 1
#NumberofParameters = 9
NumberofParameters = 7


#####------ Initial Conditions:
AmInput = input('Enter [A]_M concentration (mM) ') # A internal concentration (mM) ')

#####------ Parameters:
#AnInput = input('Enter [A]_N concentration (mM) ') #A external concentration (mM) ')
#BmInput = input('Enter [B]_M (mM) (B internal concentration) ')
J_A_Max = input('Enter J_A, pump_Max  (mM_per_s) ')
K_A_M = input('Enter K_A_M (mM) ') 
BnInput = input('Enter [B]_N (mM) (B external concentration) ')
K_B_N = input('Enter K_B_N (mM) ')

a = input('Enter A stoichiometry ratio ') 
b = input('Enter B stoichiometry ratio ') 
#c = input('Enter C stoichiometry ratio ')  


def createLegends():
	legend_ICs = [""] * NumberofICs
	legend_equation = [""] * NumberofEquations
	legend_time = ""
	legend_parameters = [""] * NumberofParameters
	legend_time = "time in component environment (second)"
	legend_ICs[0] = "[A]_M, A_int concentrations (mM)" #"Na_int in component concentrations (mM)"
	#legend_parameters[0] = "[A]_N, A_ext concentrations (mM)" #"Na_ext in component concentrations (mM)"
	legend_parameters[0] = "K_A_M (mM)" #K_Na in component NaK (mM)
	#legend_parameters[1] = "[B]_M, B_int concentrations (mM)" # "K_int in component concentrations (mM)"
	legend_parameters[1] = "K_B_N (mM)" #K_K in component NaK (mM)
	legend_parameters[2] = "[B]_N, B_ext concentrations (mM)" # "K_ext in component concentrations (mM)"
	legend_parameters[3]  = "J_A_Max in component NaK (mM_per_s)"
	#legend_parameters[4] = "K_A (mM)" #K_Na in component NaK (mM)
	#legend_parameters[5] = "K_B (mM)" #K_K in component NaK (mM)
	legend_parameters[4] = "A stoichiometric"
	legend_parameters[5] = "B stoichiometric"	
	legend_equation[0] = "J_pump_A (mM_per_s)" #"J_NaK_Na in component NaK (mM_per_s)"
	legend_equation[1] = "J_pump_B in component NaK (mM_per_s)" # "J_NaK_K in component NaK (mM_per_s)"
	return (legend_ICs, legend_equation, legend_time, legend_parameters)
def initConsts():
	parameters = [0.0] * NumberofParameters; ICs = [0.0] * NumberofICs;
	ICs[0] = float(AmInput) #0.0
	#parameters[0] = float(AnInput) # 0.0 Na_ext in component concentrations (mM)
	parameters[0] = float(K_A_M) 
	#parameters[1] = float(BmInput) # 39.0 K_int in component concentrations (mM)
	parameters[1] = float(K_B_N)
	parameters[2] = float(BnInput) # 10.0 K_ext in component concentrations (mM)
	parameters[3] = float(J_A_Max) # 10.8 (mM_per_s)
	#parameters[4] = 0.200000*(1.00000+parameters[1]/8.33000) #K_Na in component NaK (mM)
	#parameters[5] = 0.100000*(1.00000+parameters[0]/18.5000) #K_K in component NaK (mM)
	parameters[4] = float(a)
	parameters[5] = float(b)	
	parameters[6] = 1.00000 # d/dt Na_int in component concentrations (mM)
	#parameters[T5] = float(c)

    #parameters[4] = 0.200000*parameters[5] 
	return (ICs, parameters)

def computeRates(time, ICs, parameters):
	rates = [0.0] * NumberofICs; equation = [0.0] * NumberofEquations
	rates[0] = parameters[6] # d/dt Na_int in component concentrations (mM)
	return(rates)

def computeequation(parameters, ICs, time):
	equation = array([[0.0] * len(time)] * NumberofEquations)
	ICs = array(ICs)
	time = array(time)
	equation[0] = parameters[3] *(power(ICs[0]/(ICs[0]+parameters[0] ), parameters[4]))*(power(parameters[2]/(parameters[2]+parameters[1] ), parameters[5]))
	equation[1] = ((-parameters[5]/parameters[4])*equation[0])
	return equation

def solve_model():
	"""Solve model with ODE solver"""
	# Initialise parameters and state variables
	(init_ICs, parameters) = initConsts()
	# Set timespan to solve over
	time = linspace(0, 10, 500)

    # Construct ODE object to solve
	r = ode(computeRates)
	r.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=1)
	r.set_initial_value(init_ICs, time[0])
	r.set_f_params(parameters)

	# Solve model
	ICs = array([[0.0] * len(time)] * NumberofICs)
	ICs[:,0] = init_ICs
	for (i,t) in enumerate(time[1:]):
		if r.successful():
			r.integrate(t)
			ICs[:,i+1] = r.y
		else:
			break
	# Compute equation variables
	equation = computeequation(parameters, ICs, time)
	return (time, ICs, equation)
if __name__ == "__main__":
    (time, ICs, equation) = solve_model()

#def plot_model(time, ICs, equation):
#    """Plot variables against variable of integration"""
#    (legend_ICs, legend_equation, legend_time, legend_parameters) = createLegends()
#    pylab.figure(1)
#    pylab.plot(time,vstack((ICs,equation)).T)
#    pylab.xlabel(legend_time)
#    pylab.legend(legend_ICs + legend_equation, loc='best')
#    pylab.show()

#if __name__ == "__main__":
#    (time, ICs, equation) = solve_model()
#    plot_model(time, ICs, equation)

#_______________________________________________Saving in csv file_______________________________________________________#
# Time:
time = pd.Series(time)
#print time

# ICs: 
AM = pd.Series(ICs[0])	 	# (50.0) A_int concentrations (mM)
#Equations:
J_pump_A = pd.Series(equation[0]) #"J_NaK_Na in component NaK (mM_per_s)"
J_pump_B = pd.Series(equation[1]) # "J_NaK_K in component NaK (mM_per_s)"
df = pd.DataFrame({'time':time, 'A_int(mM)':AM, 'J_BMN':J_pump_B, 'J_AMN':J_pump_A})

#print("J_A = " , equation[0])
#np.savetxt("foo.csv", Et, delimiter=",")
###### output as csv file or table as well
column_order = ['time','A_int(mM)','J_AMN', 'J_BMN']
df.to_csv('unordered' + '.csv', sep = ',')
s3 = input('Do you want to save data in a separate csv file?'+'\n'+' if yes, enter YES ')
if (str(s3) == 'YES'):
	out3 = input('Choose filename: (do not include extension) ')
#	df.to_csv(str(out3) + '.csv', sep = ',')
	df[column_order].to_csv(str(out3) + '_pumptwo' +'.csv', sep = ',',index=False)
	print('The new csv file has been saved under the name of ' + str(out3) + '_pumptwo' +'.csv')
else:
	print('Ok, Continue')
#_______________________________________________Plot_______________________________________________________#

# to save the plot:
out6 = input('Choose a name to save the plot, with no extension: ') 
def plot_model(time, ICs, equation):
    """Plot variables against variable of integration"""
    (legend_ICs, legend_equation, legend_time, legend_parameters) = createLegends()
    #plt.figure(1)
    fig = plt.figure(figsize=(17,5))
    plt.plot(time,vstack((ICs,equation)).T)
    plt.xlabel(legend_time)
    plt.legend(legend_ICs + legend_equation, loc='best')
    plt.savefig(str(out6) + '_pumptwo' + '.png')
    plt.savefig(str(out6) + '_pumptwo',format = 'eps', dpi=1000, bbox_inches='tight')
    plt.show()
if __name__ == "__main__":
	plot_model(time, ICs, equation)
	print ('\n',  'The plot has been saved. ' ,'\n')
	#print  'Rerun the script if new study is desired' , '\n'


#________________________________________Post Processing:
if (str(s3) == 'YES'):
	post = input('Do you want to use results from different transporters?'+'\n'+' if yes please enter YES ')
	if (str(post) == 'YES'):
		postprocessing = open("postprocessing.txt","a")#append mode 
		postprocessing.write(str(out3) + '_SymSlp' +'.csv' + "\n") 
		postprocessing.close() 
		#print ("Run the postprocessing.py to get the net result.")
		#print "The results has been saved in postprocessing.txt file.", '\n' , "Run the postprocessing.py to get the net result."
		print ( ' running stamp_general.py script ... ')
		os.system('python3 stamp_general.py')
	else:
		print('End of simulation. Run the postprocessing.py to get the net result.' +'\n' + "Exitting")

else:
	rerun = input('Do you want to start a new simulation?'+'\n'+' if yes please enter YES  ') 

#rerun = input('If a new study is desired,please enter YES ') 

if (str(rerun) == 'YES'):
	print ( ' running stamp_general.py script ... ')
	os.system('python3 stamp_general.py')
else:
	print('End of simulation. Exitting.')






