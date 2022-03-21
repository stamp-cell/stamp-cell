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

NumberofEquations = 12
NumberofICs = 2
NumberofParameters = 12


print ('\n' + "Please enter number for following inputs:" + '\n') 
#################### Add something like if the answer true then continue
#####------ Initial Conditions:
AmInput = input('Enter A internal concentration (mM) ')
AnInput = input('Enter A external concentration (mM) ')

#####------ Parameters:
K_Am = input('Enter K_Am (mM) ')  #A internal dissocication constant (mM) 
K_An = input('Enter K_An (mM) ')  #A external dissocication constant (mM) 

g_E_N = input('Enter g_E_N (per_s) , (E translocation rate constant from ext to int) ')
g_EA_N = input('Enter g_EA_N (per_s) , (EA translocation rate constant from ext to int) ')
g_E_M = input('Enter g_E_M (per_s) , (E translocation rate constant from int to ext) ')
g_EA_M = input('Enter g_EA_M (per_s) , (EA translocation rate constant from int to ext) ')

#M side: 
#def terms(Am,k_Am,Bm,k_Bm):
#	alpha_m = Am/k_Am
#	beta_m = Bm/k_Bm
#	return(alpha_m,beta_m)
#print terms(3,4,5,6)


def createLegends():
	legend_ICs = [""] * NumberofICs
	legend_equation = [""] * NumberofEquations
	legend_time = ""
	legend_parameters = [""] * NumberofParameters
	legend_time = "time in component environment (second)"
	legend_ICs[0] = "[A]_M, A_int concentrations (mM)"
	legend_ICs[1] = "[A]_N, A_ext concentrations (mM)"
	legend_parameters[0] = "K_A_M (mM)" # "K_Cl_int in component AE1 (mM)"
	legend_parameters[1] = "K_A_N (mM)" # "K_Cl_ext in component AE1 (mM)"
	legend_parameters[2] = "g_E_N (per_s)" 
	legend_parameters[3] = "g_EA_N (per_s)"
	legend_parameters[4] = "g_E_M (per_s)" 
	legend_parameters[5] = "g_EA_M (per_s)"
	legend_parameters[8] = "E_Tmax (mM)"
	legend_parameters[9] = "K_I  (mM)"
	legend_equation[0] = "E_T (mM)"
	legend_equation[1] = "alpha_M (dimensionless)"
	legend_equation[2] = "alpha_N (dimensionless)"
	legend_equation[3] = "R_N" # 1 + alpha_n
	legend_equation[4] = "R_MM" # g_E_M + g_EA_M*alpha_M
	legend_equation[5] = "R_M" # 1 + alpha_m
	legend_equation[6] = "R_NN" # g_E_N + g_EA_N*alpha_N
	legend_equation[7] = "J Denominator" # R_N*R_MM+R_M*R_NN
	legend_equation[8] =  "J_A First_Numerator"
	legend_equation[9] =  "J_A Second_Numerator"
	legend_equation[10] = "J_A Numerator"
	legend_equation[11] = "J_A(M, N), J_A_influx (mM_per_s)"
	return (legend_ICs, legend_equation, legend_time, legend_parameters)


#_____________________________________________________________________
def ic_parameters():
	parameters = [0.0] * NumberofParameters; ICs = [0.0] * NumberofICs;
	ICs[0] = float(AmInput) 	# A_M 	# int concentrations (mM)
	ICs[1] = float(AnInput)  # A_N 	# ext concentrations (mM)
# Model Parameters:
	parameters[0] = float(K_Am) 
	parameters[1] = float(K_An) 
	parameters[2] = float(g_E_N) 
	parameters[3] = float(g_EA_N)
	parameters[4] = float(g_E_M) 
	parameters[5] = float(g_EA_M)
	parameters[8] = 1	#E_Tmax (mM)
	parameters[9] = 172	#K_I  (mM)
	parameters[10] = 0.00000
	parameters[11] = 1.00000
	#parameters[12] = 0.00000
	#parameters[13] = 0.00000
	return (ICs, parameters)

def computeRates(time, ICs, parameters):
    rates = [0.0] * NumberofICs; equation = [0.0] * NumberofEquations
    rates[0] = parameters[10] #d/dt HCO3_int in component concentrations (mM)
    rates[1] = parameters[11] #d/dt HCO3_ext in component concentrations (mM)
    #rates[2] = parameters[12] #d/dt Cl_int in component concentrations (mM)
    #rates[3] = parameters[13] #d/dt Cl_ext in component concentrations (mM)
    return(rates)

def computeEquation(parameters, ICs, time):
	equation = array([[0.0] * len(time)] * NumberofEquations)
	ICs = array(ICs)
	time = array(time)
# 0) E_T in component AE1 (mM):
	equation[0] = parameters[8]/(1.00000+ICs[0]/parameters[9])
# 1) alpha_M (dimensionless) # gamma_int in component AE1 (dimensionless)    
	equation[1] = ICs[0]/parameters[0]
# 2) alpha_N (dimensionless) # gamma_ext in component AE1 (dimensionless):    
	equation[2] = ICs[1]/parameters[1]
# 3) R_N = 1 + alpha_n 
	equation[3] = (1.00000 + equation[2])
# 4) R_MM = g_E_M + g_EA_M*alpha_M
	equation[4] = (parameters[4] + parameters[5]*equation[1])
# 5) R_M = 1 + alpha_m 
	equation[5] = (1.00000 + equation[1])
# 6) R_NN = g_E_N + g_EA_N*alpha_N
	equation[6] = (parameters[2] + parameters[3]*equation[2])
# 7) J Denominator: R_N*R_MM+R_M*R_NN
	equation[7] = (equation[3]*equation[4])+(equation[5]*equation[6])
# 8) J_A First_Numerator: = (g_EA_M*alpha_M)*g_E_N
	equation[8] =  (parameters[5]*equation[1])*parameters[2]
# 9) J_A Second_Numerator: = (g_EA_N*alpha_N)*g_E_M
	equation[9] =  (parameters[3]*equation[2])*parameters[4]
# 10) J_A Numerator:
	equation[10] = equation[8]-equation[9]
# ) J_A(M, N):
	equation[11] = equation[0] * (equation[10]/equation[7])
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


if __name__ == "__main__":
    (time, ICs, equation) = solve_model()

#_______________________________________________Saving in csv file_______________________________________________________#
# Time:
time = pd.Series(time)
#print time
# ICs: 
AM = pd.Series(ICs[0])	 	# (50.0) A_int concentrations (mM)
AN = pd.Series(ICs[1])	 	# (0.0) A_ext  concentrations (mM)
#Equations:
Et = pd.Series(equation[0])
alphaM = pd.Series(equation[1])
alphaN = pd.Series(equation[2])
RN = pd.Series(equation[3])
RMM = pd.Series(equation[4])
RM = pd.Series(equation[5])
RNN = pd.Series(equation[6])
JDenominator = pd.Series(equation[7]) #R_M*R_NN + R_N*R_MM
J_AFirst_Numerator = pd.Series(equation[8])
J_ASecond_Numerator = pd.Series(equation[9])
J_ANumerator = pd.Series(equation[10])
J_AMN = pd.Series(equation[11])

df = pd.DataFrame({'time':time, 'A_int(mM)':AM, 'A_ext(mM)':AN,'Et':Et,'alphaM':alphaM,'alphaN':alphaN,
'RN':RN,'RMM':RMM,'RM':RM,'RNN':RNN,'JDenominator':JDenominator,'J_AFirst_Numerator':J_AFirst_Numerator\
,'J_ASecond_Numerator':J_ASecond_Numerator,'J_ANumerator':J_ANumerator,'J_AMN':J_AMN})
#print("J_A = " , equation[0])
#np.savetxt("foo.csv", Et, delimiter=",")
###### output as csv file or table as well
column_order = ['time', 'A_int(mM)', 'A_ext(mM)','J_AMN','Et',\
'alphaM','alphaN','RN','RMM','RM'\
,'RNN','JDenominator','J_AFirst_Numerator'\
,'J_ASecond_Numerator','J_ANumerator']
df.to_csv('unordered' + '.csv', sep = ',')
s3 = input('Do you want to save data in a separate csv file?'+'\n'+' if yes, enter YES ')
if (str(s3) == 'YES'):
	out3 = input('Choose filename: (do not include extension) ')
#	df.to_csv(str(out3) + '.csv', sep = ',')
	df[column_order].to_csv(str(out3) + '_unidet' +'.csv', sep = ',',index=False)
	print('The new csv file has been saved under the name of ' + str(out3) + '_unidet' +'.csv')
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
    plt.savefig(str(out6) + '_unidet'+ '.png')
    plt.savefig(str(out6) + '_unidet',format = 'eps', dpi=1000, bbox_inches='tight')
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




