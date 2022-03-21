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
NumberofEquations = 20
NumberofICs = 6
NumberofParameters = 20

print ('\n' + "Please enter number for following inputs:" + '\n') 
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
K_ABCm = input('Enter K_ABCm (mM) ') #C internal dissocication constant (mM) 
K_ABCBm = input('Enter K_ABCBm (mM) ') #B internal dissocication constant (mM) 
K_ABCBn = input('Enter K_ABCBn (mM) ') #A external dissocication constant (mM) 
K_BCBn = input('Enter K_BCBn (mM) ') #B external dissocication constant (mM) 
K_CBn = input('Enter K_CBn (mM) ') #C external dissocication constant (mM) 
K_Bn = input('Enter K_Bn (mM) ') #B external dissocication constant (mM) 

g_E_N = input('Enter g_E_N (per_s) , (E translocation rate constant from ext to int) ')
g_EABCB_N = input('Enter g_EABCB_N (per_s) , (EABCB translocation rate constant from ext to int) ')
g_E_M = input('Enter g_E_M (per_s) , (E translocation rate constant from int to ext) ')
g_EABCB_M = input('Enter g_EABCB_M (per_s) , (EABCB translocation rate constant from int to ext) ')
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
	legend_ICs[2] = "[B]_M, B_int concentrations (mM)"
	legend_ICs[3] = "[B]_N, B_ext concentrations (mM)"
	legend_ICs[4] = "[C]_M, C_int concentrations (mM)"
	legend_ICs[5] = "[C]_N, C_ext concentrations (mM)"
	legend_parameters[0] = "K_Am (mM)" 
	legend_parameters[1] = "K_ABm (mM)" 
	legend_parameters[2] = "K_ABCm (mM)"
	legend_parameters[3] = "K_ABCBm (mM)"
	legend_parameters[4] = "K_ABCBn (mM)" 
	legend_parameters[5] = "K_BCBn (mM)" 
	legend_parameters[6] = "K_CBn (mM)"
	legend_parameters[7] = "K_Bn (mM)"
	legend_parameters[8] = "g_E_N (per_s)" 
	legend_parameters[9] = "g_EABCB_N (per_s" 
	legend_parameters[10] = "g_E_M (per_s)"
	legend_parameters[11] = "g_EABCB_M (per_s)"
	legend_parameters[12] = "E_Tmax (mM)"
	legend_parameters[13] = "K_I  (mM)"
	legend_equation[0] = "E_T (mM)"
	legend_equation[1] = "alpha_M (dimensionless)"
	legend_equation[2]=  "beta_M (dimensionless)"
	legend_equation[3] = "gamma_M (dimensionless)"
	legend_equation[4] = "betazegond_M (dimensionless)"
	legend_equation[5] = "alpha_N (dimensionless)"
	legend_equation[6]=  "beta_N (dimensionless)"
	legend_equation[7] = "gamma_N (dimensionless)"
	legend_equation[8] = "betaazegond_N (dimensionless)"    
	legend_equation[9] = "R_N" # 1 + betazegond_N + gamma_N*betazegond_N + beta_N*gamma_N*betazegond_N + alpha_N*beta_N*gamma_N*betazegond_N
	legend_equation[10] = "R_MM" # g_E_M + g_EABCB_M*alpha_M*beta_M*gamma_M*betazegond_M
	legend_equation[11] = "R_M" # 1 + alpha_M + alpha_M*beta_M + alpha_M*beta_M*gamma_M + alpha_M*beta_M*gamma_M*betazegond_M
	legend_equation[12] = "R_NN" # g_E_N + g_EABCB_N*alpha_N*beta_N*gamma_N*betazegond_N
	legend_equation[13] = "J Denominator" # R_M*R_NN + R_N*R_MM
	legend_equation[14] = "J_A First_Numerator" # (g_EABCB_M*alpha_M*beta_M*gamma_M*betazegond_M*g_E_N)
	legend_equation[15] = "J_A Second_Numerator" # (g_EABCB_N*alpha_N*beta_N*gamma_N*betazegond_N*g_E_M)
	legend_equation[16] = "J_ANumerator"
	legend_equation[17] = "J_A(M, N), A_influx (mM_per_s)"
	legend_equation[18] = "J_B(M, N), B_influx (mM_per_s)"  
	legend_equation[19] = "J_C(M, N), C_influx (mM_per_s)" 
	return (legend_ICs, legend_equation, legend_time, legend_parameters)

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
	parameters[2] = float(K_ABCm)
	parameters[3] = float(K_ABCBm)
	parameters[4] = float(K_ABCBn) 
	parameters[5] = float(K_BCBn) 
	parameters[6] = float(K_CBn)
	parameters[7] = float(K_Bn)
	parameters[8] = float(g_E_N) 
	parameters[9] = float(g_EABCB_N) 
	parameters[10] = float(g_E_M)
	parameters[11] = float(g_EABCB_M)
	parameters[12] = 1	#E_Tmax (mM)
	parameters[13] = 172	#K_I  (mM)
	parameters[14] = 0.00000
	parameters[15] = 1.00000
	parameters[16] = 0.00000
	parameters[17] = 0.00000
	parameters[18] = 0.00000
	parameters[19] = 0.00000
	return (ICs, parameters)

def computeRates(time, ICs, parameters):
	rates = [0.0] * NumberofICs; equation = [0.0] * NumberofEquations
	rates[0] = parameters[14] #d/dt A_M in component concentrations (mM) int=M
	rates[1] = parameters[15] #d/dt A_N in component concentrations (mM) ext=N
	rates[2] = parameters[16] #d/dt B_M in component concentrations (mM)
	rates[3] = parameters[17] #d/dt B_N in component concentrations (mM)
	rates[4] = parameters[18] #d/dt C_M in component concentrations (mM)
	rates[5] = parameters[19] #d/dt C_N in component concentrations (mM)
	return(rates)

def computeEquation(parameters, ICs, time):
	equation = array([[0.0] * len(time)] * NumberofEquations)
	ICs = array(ICs)
	time = array(time)
# 0) E_T in component AE1 (mM):
	equation[0] = parameters[12]/(1.00000+ICs[0]/parameters[13])
# 1) alpha_M (dimensionless)     
	equation[1] = ICs[0]/parameters[0]
# 2) beta_M (dimensionless)     
	equation[2]=  ICs[2]/parameters[1]
# 3) gamma_M (dimensionless)    
	equation[3] = ICs[4]/parameters[2]
# 4) betazegond_M (dimensionless)    
	equation[4] = ICs[2]/parameters[3]
# 5) alpha_N (dimensionless) #     
	equation[5] = ICs[1]/parameters[4]
# 6) beta_N (dimensionless)     
	equation[6]=  ICs[3]/parameters[5]
# 7) gamma_N (dimensionless)    
	equation[7] = ICs[5]/parameters[6]
# 8) betaazegond_N (dimensionless)    
	equation[8] = ICs[3]/parameters[7]
# 9) R_N = 1 + betazegond_N + gamma_N*betazegond_N + beta_N*gamma_N*betazegond_N + alpha_N*beta_N*gamma_N*betazegond_N
	equation[9] = (1 + equation[8] + equation[7]*equation[8] + equation[6]*equation[7]*equation[8] + equation[5]*equation[6]*equation[7]*equation[8])  
# 10) R_MM = g_E_M + g_EABCB_M*alpha_M*beta_M*gamma_M*betazegond_M
	equation[10] = (parameters[10] + parameters[11]*equation[1]*equation[2]*equation[3]*equation[4])
# 11) R_M = 1 + alpha_M + alpha_M*beta_M + alpha_M*beta_M*gamma_M + alpha_M*beta_M*gamma_M*betazegond_M
	equation[11] = (1 + equation[1] + equation[1]*equation[2] + equation[1]*equation[2]*equation[3] + equation[1]*equation[2]*equation[3]*equation[4])
# 12) R_NN = g_E_N + g_EABCB_N*alpha_N*beta_N*gamma_N*betazegond_N
	equation[12] = (parameters[8] + parameters[9]*equation[5]*equation[6]*equation[7]*equation[8])
# 13) J Denominator: R_M*R_NN + R_N*R_MM
	equation[13] = (equation[11]*equation[12] + equation[9]*equation[10])
# 14) J_A First_Numerator: = (g_EABCB_M*alpha_M*beta_M*gamma_M*betazegond_M*g_E_N)
	equation[14] = (parameters[11]*equation[1]*equation[2]*equation[3]*equation[4]*parameters[8])
# 15) J_A Second_Numerator: = (g_EABCB_N*alpha_N*beta_N*gamma_N*betazegond_N*g_E_M)
	equation[15] = (parameters[9]*equation[5]*equation[6]*equation[7]*equation[8]*parameters[10])
# 16) J_A Numerator:
	equation[16] = (equation[15]-equation[14])
# 17) J_A(M, N):
	equation[17] = equation[0] * (equation[16]/equation[13])
# 18) J_B(M, N):
	equation[18] = 2*equation[17] 
# 19) J_C(M, N):
	equation[19] = equation[17] 
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
AN = pd.Series(ICs[1])	 	# (0.0) A_ext  concentrations (mM)
BM = pd.Series(ICs[2])		# (0.0) B_int  concentrations (mM)
BN = pd.Series(ICs[3])		# (0.0) B_ext  concentrations (mM)
CM = pd.Series(ICs[4])		# C_int concentrations (mM)
CN = pd.Series(ICs[5])		# C_ext concentrations (mM)

#Equations:
Et = pd.Series(equation[0])
alphaM = pd.Series(equation[1])
betaM = pd.Series(equation[2])
gammaM = pd.Series(equation[3])
betazegondM = pd.Series(equation[4])
alphaN = pd.Series(equation[5])
betaN = pd.Series(equation[6])
gammaN = pd.Series(equation[7])
betazegondN = pd.Series(equation[8])
RN = pd.Series(equation[9])
RMM = pd.Series(equation[10])
RM = pd.Series(equation[11])
RNN = pd.Series(equation[12])
JDenominator = pd.Series(equation[13]) #R_M*R_NN + R_N*R_MM
J_AFirst_Numerator = pd.Series(equation[14])
J_ASecond_Numerator = pd.Series(equation[15])
J_ANumerator = pd.Series(equation[16])
J_AMN = pd.Series(equation[17])
J_BMN = pd.Series(equation[18])
J_CMN = pd.Series(equation[19])

df = pd.DataFrame({'time':time, 'A_int(mM)':AM, 'A_ext(mM)':AN, 'B_int(mM)':BM,'B_ext(mM)' :BN,\
'C_int(mM)':CM, 'C_ext(mM)':CN, 'Et':Et,'J_AMN':J_AMN, 'J_BMN':J_BMN,'J_CMN':J_CMN,'alphaM':alphaM,'alphaN':alphaN,'betaM':betaM,'betaN':betaN, 'gammaM':gammaM,'gammaN':gammaN,'betazegondM':betazegondM,'betazegondN':betazegondN,
'R_N':RN,'R_MM':RMM ,'R_M':RM, 'R_NN':RNN, 'JDenominator':JDenominator, 'J_AFirst_Numerator':J_AFirst_Numerator,\
'J_ASecond_Numerator':J_ASecond_Numerator, 'J_ANumerator':J_ANumerator})

#print("J_A = " , equation[0])
#np.savetxt("foo.csv", Et, delimiter=",")
###### output as csv file or table as well
column_order = ['time', 'A_int(mM)', 'A_ext(mM)','B_int(mM)','B_ext(mM)','C_int(mM)', 'C_ext(mM)',\
'J_AMN', 'J_BMN','J_CMN','Et','alphaM','betaM',
'alphaN','betaN','betazegondM','betazegondN','gammaM','gammaN','R_N','R_MM','R_M','R_NN',\
'JDenominator','J_AFirst_Numerator','J_ASecond_Numerator','J_ANumerator']
df.to_csv('unordered' + '.csv', sep = ',')
s3 = input('Do you want to save data in a separate csv file?'+'\n'+' if yes, enter YES ')
if (str(s3) == 'YES'):
	out3 = input('Choose filename: (do not include extension) ')
#	df.to_csv(str(out3) + '.csv', sep = ',')
	df[column_order].to_csv(str(out3) + '_SymABCB' +'.csv', sep = ',',index=False)
	print('The new csv file has been saved under the name of ' + str(out3) + '_SymABCB' +'.csv')
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
    plt.savefig(str(out6) + '_SymABCB' + '.png')
    plt.savefig(str(out6) + '_SymABCB',format = 'eps', dpi=1000, bbox_inches='tight')
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







