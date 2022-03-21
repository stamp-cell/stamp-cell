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

# Size of variable arrays:
NumberofEquations = 1
NumberofICs = 1
NumberofParameters = 3

print ('\n' + "Please enter number for following inputs:" + '\n') 
#####------ Initial Conditions:
AmInput = input('Enter [A]_M concentration (mM) ') # A internal concentration (mM) ')
#####------ Parameters:
V_max = input('Enter V_max  (mM_per_s) ')
K_M = input('Enter K_M (mM) ') 
def createLegends():
	legend_ICs = [""] * NumberofICs
	legend_equation = [""] * NumberofEquations
	legend_time = ""
	legend_parameters = [""] * NumberofParameters
	legend_time = "time in component environment (second)"
	legend_ICs[0] = "[A]_M, A_int in component concentrations (mM)" #"Na_int in component concentrations (mM)"
	#legend_parameters[0] = "Na_ext in component concentrations (mM)"
	legend_parameters[0] = "V_max (mM_per_s)" #"J_NaK_Na_Max in component NaK (mM_per_s)"
	legend_parameters[2] = "K_M (mM)"
	legend_equation[0] = "J_A,pump"
	#legend_rates[0] = "d/dt Na_int in component concentrations (mM)"
	return (legend_ICs, legend_equation, legend_time, legend_parameters)

def initConsts():
	parameters = [0.0] * NumberofParameters; ICs = [0.0] * NumberofICs;
	ICs[0] = float(AmInput) #0.0
	parameters[0] = float(V_max) # 10.8 (mM_per_s)
	parameters[1] = 1.00000 # d/dt Na_int in component concentrations (mM)
	parameters[2] = float(K_M)
	return (ICs, parameters)

def computeRates(time, ICs, parameters):
	rates = [0.0] * NumberofICs; equation = [0.0] * NumberofEquations
	rates[0] = parameters[1] # d/dt Na_int in component concentrations (mM)
	return(rates)

def computeequation(parameters, ICs, time):
	equation = array([[0.0] * len(time)] * NumberofEquations)
	ICs = array(ICs)
	time = array(time)
# 0) J_A,pump = J_A,pump,Max *(A_M/(A_M + K_M)) 
	equation[0] = parameters[0]*(ICs[0]/(ICs[0]+parameters[2]))
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

#_______________________________________________Saving in csv file_______________________________________________________#
# Time:
time = pd.Series(time)
#print time
# ICs: 
AM = pd.Series(ICs[0])	 	# (50.0) A_int concentrations (mM)
#Equations:
J_A_pump = pd.Series(equation[0])

df = pd.DataFrame({'time':time,'J_AMN':J_A_pump, 'A_int(mM)':AM})
#print("J_A = " , equation[0])
#np.savetxt("foo.csv", Et, delimiter=",")
###### output as csv file or table as well
column_order = ['time', 'A_int(mM)','J_AMN']
df.to_csv('unordered' + '.csv', sep = ',')
s3 = input('Do you want to save data in a separate csv file?'+'\n'+' if yes, enter YES ')
if (str(s3) == 'YES'):
	out3 = input('Choose filename: (do not include extension) ')
#	df.to_csv(str(out3) + '.csv', sep = ',')
	df[column_order].to_csv(str(out3) + '_uniMM' +'.csv', sep = ',',index=False)
	print('The new csv file has been saved under the name of ' + str(out3) + '_uniMM' +'.csv')
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
    plt.savefig(str(out6) + '_uniMM'+ '.png')
    plt.savefig(str(out6) + '_uniMM',format = 'eps', dpi=1000, bbox_inches='tight')
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







