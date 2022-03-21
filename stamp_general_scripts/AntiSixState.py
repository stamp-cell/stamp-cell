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
NumberofEquations = 15
NumberofICs = 4
NumberofParameters = 14

print ('\n' + "Please enter number for following inputs:" + '\n') 
#################### Add something like if the answer true then continue
#####------ Initial Conditions:
AmInput = input('Enter [A]_M concentration (mM) ') # A internal concentration (mM) ')
AnInput = input('Enter [A]_N concentration (mM) ') #A external concentration (mM) ')
BmInput = input('Enter [B]_M concentration (mM) ') #B internal concentration (mM) ')
BnInput = input('Enter [B]_N concentration (mM) ') #B external concentration (mM) ')

#####------ Parameters:
K_Am = input('Enter K_Am (mM) ') #('Enter A internal dissocication constant (mM) ')
K_An = input('Enter K_An (mM) ') #('Enter A external dissocication constant (mM) ')
K_Bm = input('Enter K_Bm (mM) ') #('Enter B internal dissocication constant (mM) ')
K_Bn = input('Enter K_Bn (mM) ') #('Enter B external dissocication constant (mM) ')


g_EB_N = input('Enter g_EB_N (per_s) , (EB translocation rate constant from ext to int) ')
g_EB_M = input('Enter g_EB_M (per_s) , (EB translocation rate constant from int to ext) ')
g_EA_N = input('Enter g_EA_N (per_s) , (EA translocation rate constant from ext to int) ')
g_EA_M = input('Enter g_EA_M (per_s) , (EA translocation rate constant from int to ext) ')


def createLegends():
	legend_ICs = [""] * NumberofICs
	legend_equation = [""] * NumberofEquations
	legend_time = ""
	legend_parameters = [""] * NumberofParameters
	legend_time = "time in component environment (second)"
	legend_ICs[0] = "HCO3_int in component concentrations (mM)" #A_int concentrations (mM)
	legend_ICs[1] = "HCO3_ext in component concentrations (mM)" #A_ext  concentrations (mM)
	legend_ICs[2] = "Cl_int in component concentrations (mM)" #B_int  concentrations (mM)
	legend_ICs[3] = "Cl_ext in component concentrations (mM)" #B_ext  concentrations (mM)
	legend_equation[0] = "E_T (mM)"
	legend_equation[1] = "alpha_M" # "gamma_int in component AE1 (dimensionless)"
	legend_equation[2] = "beta_M" # "beta_int in component AE1 (dimensionless)"
	legend_equation[3] = "alpha_N" # "gamma_ext in component AE1 (dimensionless)"
	legend_equation[4] = "beta_N" # "beta_ext in component AE1 (dimensionless)"
	legend_equation[5] = "R_N" # 
	legend_equation[6] = "R_MM" # 
	legend_equation[7] = "R_M" # 
	legend_equation[8] = "R_NN" # 
	legend_equation[9] = "R_M*R_NN + R_N*R_MM" # 
	legend_equation[10] = "(g_EA_M*alpha_M)*(g_EB_N*beta_N)" # 
	legend_equation[11] = " (g_EA_N*alpha_N)*(g_EB_M*beta_M)" # 
	legend_equation[12] = "(g_EA_M*alpha_M)*(g_EB_N*beta_N)-(g_EA_N*alpha_N)*(g_EB_M*beta_M)"#J_A Numerator
	legend_equation[13] = "J_A_influx in component AE1 (mM_per_s)" #J_A(M, N)
	legend_equation[14] = "J_B_influx in component AE1 (mM_per_s)" #J_B(M, N): = -J_A(M, N)
	legend_parameters[0] = "K_A_M (mM)" # "K_Cl_int in component AE1 (mM)"
	legend_parameters[1] = "K_B_M (mM)" # "K_HCO3_int in component AE1 (mM)"
	legend_parameters[2] = "K_A_N (mM)" # "K_Cl_ext in component AE1 (mM)"
	legend_parameters[3] = "K_B_N (mM)" # "K_HCO3_ext in component AE1 (mM)"
	legend_parameters[4] = "g_EB_N (per_s)" # "P_HCO3_ext in component AE1 (per_s)"
	legend_parameters[5] = "g_EB_M (per_s)" # "P_HCO3_int in component AE1 (per_s)"
	legend_parameters[6] = "g_EA_N (per_s)" # "P_Cl_ext in component AE1 (per_s)"
	legend_parameters[7] = "g_EA_M (per_s)" # "P_Cl_int in component AE1 (per_s)"
	legend_parameters[8] = "E_Tmax in component AE1 (mM)"
	legend_parameters[9] = "K_I in component AE1 (mM)"
	return (legend_ICs, legend_equation, legend_time, legend_parameters)

#_____________________________________________________________________
def ic_parameters():
	parameters = [0.0] * NumberofParameters; ICs = [0.0] * NumberofICs;
	ICs[0] = float(AmInput) # A_M 	# (50.0) A_int concentrations (mM)
	ICs[1] = float(AnInput) # A_N 	# (0.0) A_ext  concentrations (mM)
	ICs[2] = float(BmInput) # B_M	# (0.0) B_int  concentrations (mM)
	ICs[3] = float(BnInput) # B_N	# (0.0) B_ext  concentrations (mM)
# Model Parameters:
	parameters[0] = float(K_Am) # K_Am	#50 K"_c, K_Cl_int in component AE1 (per_s)
	parameters[1] = float(K_Bm) # K_Bm	#198 K"_b, K_HCO3_int in component AE1 (per_s)
	parameters[2] = float(K_An) # K_An	#50 K'_c, K_Cl_ext in component AE1 (per_s)
	parameters[3] = float(K_Bn)  # K_Bn	# 198 K'_b, K_HCO3_ext in component AE1 (per_s)
	parameters[4] = float(g_EB_N) # g_EB_N # 1247 P'_b, P_HCO3_ext in component AE1 (per_s)
	parameters[5] = float(g_EB_M) # g_EB_M	# 135 P"_b, P_HCO3_int in component AE1 (per_s)
	parameters[6] = float(g_EA_N) # g_EA_N	# 562 P'_c, P_Cl_ext in component AE1 (per_s)
	parameters[7] = float(g_EA_M) # g_EA_M	#61 P"_c, P_Cl_int in component AE1 (per_s)
	parameters[8] = 1	#E_Tmax (mM)
	parameters[9] = 172	#K_I  (mM)
	parameters[10] = 0.00000
	parameters[11] = 1.00000
	parameters[12] = 0.00000
	parameters[13] = 0.00000
	return (ICs, parameters)

def computeRates(time, ICs, parameters):
    rates = [0.0] * NumberofICs; equation = [0.0] * NumberofEquations
    rates[0] = parameters[10] #d/dt HCO3_int in component concentrations (mM)
    rates[1] = parameters[11] #d/dt HCO3_ext in component concentrations (mM)
    rates[2] = parameters[12] #d/dt Cl_int in component concentrations (mM)
    rates[3] = parameters[13] #d/dt Cl_ext in component concentrations (mM)
    return(rates)

def computeEquation(parameters, ICs, time):
	equation = array([[0.0] * len(time)] * NumberofEquations)
	ICs = array(ICs)
	time = array(time)
# 0) E_T in component AE1 (mM):
	equation[0] = parameters[8]/(1.00000+ICs[0]/parameters[9])
# 1) alpha_M (dimensionless) # gamma_int in component AE1 (dimensionless)    
	equation[1] = ICs[0]/parameters[0]
# 2) beta_M (dimensionless) # beta_int in component AE1 (dimensionless):    
	equation[2] = ICs[2]/parameters[1]
# 3) alpha_N (dimensionless) # gamma_ext in component AE1 (dimensionless):    
	equation[3] = ICs[1]/parameters[2]
# 4) beta_N (dimensionless)# beta_ext in component AE1 (dimensionless):
	equation[4] = ICs[3]/parameters[3]
# 5) R_N = equation[5] = 1 + alpha_n + beta_n  
	equation[5] = (1.00000 + equation[3] + equation[4])
# 6) R_MM = equation[6] = g_EA_M*alpha_m + g_EB_M*beta_M
	equation[6] = (parameters[7]*equation[1] + parameters[5]*equation[2])
# 7) R_M = equation[7] = 1 + alpha_m + beta_m
	equation[7] = (1.00000 + equation[1] + equation[2])
# 8) R_NN = equation [21] = g_EB_N*beta_N + g_EA_N*alpha_N
	equation[8] = (parameters[4]*equation[4]+parameters[6]*equation[3])
# 9) J Denominator: R_M*R_NN + R_N*R_MM
	equation[9] = (equation[7]*equation[8])+(equation[5]*equation[6])
# 10) J_A First_Numerator: = (g_EA_M*alpha_M)*(g_EB_N*beta_N)
	equation[10] = (parameters[7]*equation[1])*(parameters[4]*equation[4])
# 11) J_A Second_Numerator: = (g_EA_N*alpha_N)*(g_EB_M*beta_M)
	equation[11] = (parameters[6]*equation[3])*(parameters[5]*equation[2])
# 12) J_A Numerator:
	equation[12] = equation[10]-equation[11]
# 13) J_A(M, N):
	equation[13] = equation[0] * (equation[12]/equation[9])
# 14) J_B(M, N): = -J_A(M, N)
	equation[14] = -equation[13]
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
AN = pd.Series(ICs[1])	 	# (0.0) A_ext  concentrations (mM)
BM = pd.Series(ICs[2])		# (0.0) B_int  concentrations (mM)
BN = pd.Series(ICs[3])		# (0.0) B_ext  concentrations (mM)
#Equations:
Et = pd.Series(equation[0])
alphaM = pd.Series(equation[1])
betaM = pd.Series(equation[2])
alphaN = pd.Series(equation[3])
betaN = pd.Series(equation[4])
RN = pd.Series(equation[5])
RMM = pd.Series(equation[6])
RM = pd.Series(equation[7])
RNN = pd.Series(equation[8])
JDenominator = pd.Series(equation[9]) #R_M*R_NN + R_N*R_MM
J_AFirst_Numerator = pd.Series(equation[10])
J_ASecond_Numerator = pd.Series(equation[11])
J_ANumerator = pd.Series(equation[12])
J_AMN = pd.Series(equation[13])
J_BMN = pd.Series(equation[14])
df = pd.DataFrame({'time':time, 'A_int(mM)':AM, 'A_ext(mM)':AN, 'B_int(mM)':BM,'B_ext(mM)' :BN,\
'Et':Et,'alphaM':alphaM,'betaM':betaM,'alphaN':alphaN,'betaN':betaN,'RN':RN,'RMM':RMM,'RM':RM\
,'RNN':RNN,'JDenominator':JDenominator,'J_AFirst_Numerator':J_AFirst_Numerator\
,'J_ASecond_Numerator':J_ASecond_Numerator,'J_ANumerator':J_ANumerator,'J_AMN':J_AMN,'J_BMN':J_BMN})
#print("J_A = " , equation[0])
#np.savetxt("foo.csv", Et, delimiter=",")
###### output as csv file or table as well
column_order = ['time', 'A_int(mM)', 'A_ext(mM)','B_int(mM)','B_ext(mM)','J_AMN','J_BMN',\
'Et','alphaM','betaM','alphaN','betaN','RN','RMM','RM'\
,'RNN','JDenominator','J_AFirst_Numerator'\
,'J_ASecond_Numerator','J_ANumerator']
df.to_csv('unordered' + '.csv', sep = ',')

s3 = input('Do you want to save data in a separate csv file?'+'\n'+' if yes, enter YES ')
if (str(s3) == 'YES'):
	out3 = input('Choose filename: (do not include extension) ')
#	df.to_csv(str(out3) + '.csv', sep = ',')
	df[column_order].to_csv(str(out3) + '_AntiSixState' +'.csv', sep = ',',index=False)
	print('The new csv file has been saved under the name of ' + str(out3) + '_AntiSixState' +'.csv')
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
    plt.savefig(str(out6) + '_AntiSixState'  + '.png')
    plt.savefig(str(out6) + '_AntiSixState' ,format = 'eps', dpi=1000, bbox_inches='tight')
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







