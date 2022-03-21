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

NumberofEquations = 10
NumberofICs = 2
NumberofParameters = 10

print ('\n' + "Please enter number for following inputs:" + '\n') 
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
	legend_time = ""
	legend_parameters = [""] * NumberofParameters
	legend_time = "time in component environment (second)"
	legend_ICs[0] = "[A]_M, ,A_int in component concentrations (mM)"
	legend_ICs[1] = "[A]_N, A_ext in component concentrations (mM)"
	legend_parameters[0] = "K_EAPm (mM)" 
	legend_parameters[1] = "K_EAPn (mM)" 
	legend_parameters[2] = "K_A (mM)"
	legend_parameters[3] = "g (per_s)" 
	legend_parameters[4]  = "E_Tmax  (mM)"
	legend_parameters[5]   = "K_I (mM)"
	legend_equation[0] = "E_T (mM)"
	legend_equation[1] = "J_Denom_First" #(K_EAPm + K_EAPn + 2*K_EAPm*K_EAPn)*A_M*A_N 
	legend_equation[2] = "J_Denom_Second" #K_A*(1+2*K_EAPm)*A_M+2*K_A^2
	legend_equation[3] = "J_Denom_Third" #K_A*(1+2*K_EAPn)*A_N
	legend_equation[4] = "J_Denom_Fourth" # 2*K_A^2
	legend_equation[5] = "JDenominator"
	legend_equation[6] = "J_A First_Numerator" # (K_EAPm*A_M)
	legend_equation[7] = "J_A Second_Numerator" # (K_EAPn*A_N)
	legend_equation[8] = "J_A Numerator"
	legend_equation[9] = "J_A (mM_per_s)"
	return (legend_ICs, legend_equation, legend_time, legend_parameters)

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

def computeRates(time, ICs, parameters):
    rates = [0.0] * NumberofICs; equation = [0.0] * NumberofEquations
    rates[0] = parameters[6]   #d/dt HCO3_int in component concentrations (mM)
    rates[1] = parameters[7]   #d/dt HCO3_ext in component concentrations (mM)
    #rates[2] = parameters[8]   #d/dt Cl_int in component concentrations (mM)
    #rates[3] = parameters[9]   #d/dt Cl_ext in component concentrations (mM)
    return(rates)

def computeEquation(parameters, ICs, time):
	equation = array([[0.0] * len(time)] * NumberofEquations)
	ICs = array(ICs)
	time = array(time)
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
    equation = computeEquation(parameters, ICs, time)
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
Et = pd.Series(equation[0]) 
J_Denom_First = pd.Series(equation[1])
J_Denom_Second = pd.Series(equation[2])
J_Denom_Third =	pd.Series(equation[3])
J_Denom_Fourth = pd.Series(equation[4])
JDenominator = pd.Series(equation[5])
J_AFirst_Numerator = pd.Series(equation[6])
J_ASecond_Numerator = pd.Series(equation[7])
J_ANumerator = pd.Series(equation[8])
J_AMN  = pd.Series(equation[9])





df = pd.DataFrame({'time':time,'J_AMN':J_AMN ,'A_int(mM)':AM,'Et':Et,'J_Denom_First':J_Denom_First,'J_Denom_Second':J_Denom_Second,'J_Denom_Third':J_Denom_Third,
'J_Denom_Fourth':J_Denom_Fourth,'JDenominator':JDenominator, 'J_AFirst_Numerator':J_AFirst_Numerator,'J_ASecond_Numerator':J_ASecond_Numerator,
'J_ANumerator':J_ANumerator})
#print("J_A = " , equation[0])
#np.savetxt("foo.csv", Et, delimiter=",")
###### output as csv file or table as well
column_order = ['time','J_AMN','A_int(mM)','Et','J_Denom_First','J_Denom_Second','J_Denom_Third',
'J_Denom_Fourth','JDenominator','J_AFirst_Numerator','J_ASecond_Numerator','J_ANumerator']
df.to_csv('unordered' + '.csv', sep = ',')
s3 = input('Do you want to save data in a separate csv file?'+'\n'+' if yes, enter YES ')
if (str(s3) == 'YES'):
	out3 = input('Choose filename: (do not include extension) ')
#	df.to_csv(str(out3) + '.csv', sep = ',')
	df[column_order].to_csv(str(out3) + '_pumpone' +'.csv', sep = ',',index=False)
	print('The new csv file has been saved under the name of ' + str(out3) + '_pumpone' +'.csv')
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
    plt.savefig(str(out6) + '_pumpone' + '.png')
    plt.savefig(str(out6) + '_pumpone',format = 'eps', dpi=1000, bbox_inches='tight')
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







