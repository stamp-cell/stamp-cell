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
NumberofEquations = 24
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
K_Bm = input('Enter K_Bm (mM) ') #B internal dissocication constant (mM) 
K_Cm = input('Enter K_Cm (mM) ') #C internal dissocication constant (mM) 
K_An = input('Enter K_An (mM) ') #A external dissocication constant (mM) 
K_Bn = input('Enter K_Bn (mM) ') #B external dissocication constant (mM) 
K_Cn = input('Enter K_Cn (mM) ') #C external dissocication constant (mM) 


g_EA_N = input('Enter g_EA_N (per_s) , (E translocation rate constant from ext to int) ')
g_EB_N = input('Enter g_EB_N (per_s) , (EB translocation rate constant from ext to int) ')
g_EC_N = input('Enter g_EC_N (per_s) , (EC translocation rate constant from ext to int) ')

g_EA_M = input('Enter g_EA_M (per_s) , (EA translocation rate constant from int to ext) ')
g_EB_M = input('Enter g_EB_M (per_s) , (EB translocation rate constant from int to ext) ')
g_EC_M = input('Enter g_EC_M (per_s) , (EC translocation rate constant from int to ext) ')
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
    legend_ICs[0] = "A_M (A_int in component concentrations (mM))"
    legend_ICs[1] = "A_N (A_ext in component concentrations (mM))"
    legend_ICs[2] = "B_M (B_int in component concentrations (mM))"
    legend_ICs[3] = "B_N (B_ext in component concentrations (mM))"
    legend_ICs[4] = "C_M (C_int in component concentrations (mM))"
    legend_ICs[5] = "C_N (C_ext in component concentrations (mM))"
    legend_parameters[0] = "K_Am (mM)"
    legend_parameters[1] = "K_Bm (mM)"
    legend_parameters[2] = "K_Cm (mM)"
    legend_parameters[3] = "K_An (mM)" 
    legend_parameters[4] = "K_Bn (mM)" 
    legend_parameters[5] = "K_Cn (mM)" 
    legend_parameters[6] = "g_EA_N (per_s)" 
    legend_parameters[7] = "g_EC_N (per_s)" 
    legend_parameters[8] = "g_EB_N (per_s)" 
    legend_parameters[9] = "g_EA_M (per_s)" 
    legend_parameters[10] = "g_EC_M (per_s)"
    legend_parameters[11] = "g_EB_M (per_s)"
    legend_parameters[12] = "E_Tmax in component AE1 (mM)"
    legend_parameters[13] = "K_I in component AE1 (mM)"
    legend_equation[0] = "E_T (mM)"
    legend_equation[1] = "alpha_M" 
    legend_equation[2] = "beta_M" 
    legend_equation[3] = "gamma_M" 
    legend_equation[4] = "alpha_N"
    legend_equation[5] = "beta_N"
    legend_equation[6] = "gamma_N"
    legend_equation[7] = "R_N"
    legend_equation[8] = "R_MM"
    legend_equation[9] = "R_M"
    legend_equation[10] = "R_NN"
    legend_equation[11] = "JDenominator"
    legend_equation[12] = "J_AFirst_Numerato"
    legend_equation[13] = "J_ASecond_Numerator"
    legend_equation[14] = "J_A Numerator"
    legend_equation[15] = "J_A(M, N) A influx in component AE1 (mM_per_s)"
    legend_equation[16] = "J_B First_Numerator"
    legend_equation[17] = "J_B Second_Numerator" # g_EB_N*beta_N*(g_EA_M*alpha_M + g_EC_M*gama_M)
    legend_equation[18] = "J_B Numerator"
    legend_equation[19] = "J_B(M, N) B_influx in component AE1 (mM_per_s)"
    legend_equation[20] = "J_C First_Numerator" # g_EC_M*gama_M*(g_EA_N*alpha_N + g_EB_N*beta_N)
    legend_equation[21] = "J_C Second_Numerator"
    legend_equation[22] = "J_C Numerator C_influx in component AE1 (mM_per_s)"
    legend_equation[23] = "J_C(M, N) C_influx in component AE1 (mM_per_s)"
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
	parameters[1] = float(K_Bm) 
	parameters[2] = float(K_Cm)
	parameters[3] = float(K_An) 
	parameters[4] = float(K_Bn) 
	parameters[5] = float(K_Cn)
	parameters[6] = float(g_EA_N) 
	parameters[7] = float(g_EC_N) 
	parameters[8] = float(g_EB_N)
	parameters[9] = float(g_EA_M) 
	parameters[10] = float(g_EC_M)
	parameters[11] = float(g_EB_M)
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
# 4) alpha_N (dimensionless) #     
	equation[4] = ICs[1]/parameters[3]
# 5) beta_N (dimensionless)     
	equation[5]=  ICs[3]/parameters[3]
# 6) gamma_N (dimensionless)    
	equation[6] = ICs[5]/parameters[4]
# 7) R_N = 1 + alpha_N + beta_N + gamma_N
	equation[7] = (1.0 + equation[4] + equation[5] + equation[6])
# 8) R_MM = g_EA_M*alpha_M + g_EB_M*beta_M + g_EC_M*gamma_M
	equation[8] = (parameters[9]*equation[1] + parameters[10]*equation[2] + parameters[11]*equation[3])
# 9) R_M = 1 + alpha_M + beta_M + gamma_M
	equation[9] = (1.0 + equation[1] + equation[2] + equation[3])
# 10) R_NN = g_EA_N*alpha_N + g_EB_N*beta_N + g_EC_N*gamma_N
	equation[10] = (parameters[6]*equation[4] + parameters[8]*equation[5] + parameters[7]*equation[6])
# 11) J Denominator: R_M*R_NN + R_N*R_MM
	equation[11] = (equation[9]*equation[10]+equation[7]*equation[8])
# 12) J_A First_Numerator: = g_EA_M*alpha_M*(g_EB_N*beta_N + g_EC_N*gama_N)
	equation[12] =  parameters[9]*equation[1]*(parameters[8]*equation[5]+parameters[7]*equation[6]) 
# 13) J_A Second_Numerator: = g_EA_N*alpha_N*(g_EB_M*beta_M + g_EC_M*gama_M)
	equation[13] = parameters[6]*equation[4]*(parameters[11]*equation[2]+parameters[10]*equation[3])
# 14) J_A Numerator:
	equation[14] = (equation[12]-equation[13])
# 15) J_A(M, N):
	equation[15] = equation[0] * (equation[14]/equation[11])
# 16) J_B First_Numerator: = g_EB_M*beta_M*(g_EA_N*alpha_N + g_EC_N*gama_N)
	equation[16] = parameters[11]*equation[2]*(parameters[6]*equation[4]+parameters[7]*equation[4])
# 17) J_B Second_Numerator: = g_EB_N*beta_N*(g_EA_M*alpha_M + g_EC_M*gama_M)
	equation[17] = parameters[8]*equation[5]*(parameters[9]*equation[1]+parameters[10]*equation[3])
# 18) J_B Numerator:
	equation[18] = (equation[16]-equation[17])
# 19) J_B(M, N):
	equation[19] = equation[0] * (equation[18]/equation[11])
# 20) J_C First_Numerator: = g_EC_M*gama_M*(g_EA_N*alpha_N + g_EB_N*beta_N)
	equation[20] = parameters[10]*equation[3]*(parameters[6]*equation[4]+parameters[8]*equation[5])
# 21) J_C Second_Numerator: = g_EC_N*gama_N*(g_EA_M*alpha_M + g_EB_M*beta_M)
	equation[21] = parameters[7]*equation[6]*(parameters[9]*equation[1]+parameters[11]*equation[2])
# 22) J_C Numerator:
	equation[22] = (equation[20]-equation[21])
# 23) J_C(M, N):
	equation[23] = equation[0] * (equation[22]/equation[11])
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
alphaN = pd.Series(equation[4])
betaN = pd.Series(equation[5])
gammaN = pd.Series(equation[6])
RN = pd.Series(equation[7])
RMM = pd.Series(equation[8])
RM = pd.Series(equation[9])
RNN = pd.Series(equation[10])
JDenominator = pd.Series(equation[11]) #R_M*R_NN + R_N*R_MM
J_AFirst_Numerator = pd.Series(equation[12])
J_ASecond_Numerator = pd.Series(equation[13])
J_ANumerator = pd.Series(equation[14])
J_BFirst_Numerator = pd.Series(equation[16])
J_BSecond_Numerator = pd.Series(equation[17])
J_BNumerator = pd.Series(equation[18])
J_CFirst_Numerator = pd.Series(equation[20])
J_CSecond_Numerator = pd.Series(equation[21])
J_CNumerator = pd.Series(equation[22])
J_AMN = pd.Series(equation[15])
J_BMN = pd.Series(equation[19])
J_CMN = pd.Series(equation[23])


df = pd.DataFrame({'time':time, 'A_int(mM)':AM, 'A_ext(mM)':AN, 'B_int(mM)':BM,'B_ext(mM)' :BN,\
'C_int(mM)':CM, 'C_ext(mM)':CN, 'Et':Et,'J_AMN':J_AMN, 'J_BMN':J_BMN,'J_CMN':J_CMN,'alphaM':alphaM,'alphaN':alphaN,'betaM':betaM,'betaN':betaN,'gammaM':gammaM,'gammaN':gammaN,\
'R_N':RN,'R_MM':RMM , 'R_M':RM, 'R_NN':RNN, 'JDenominator':JDenominator, 'J_AFirst_Numerator':J_AFirst_Numerator,\
'J_ASecond_Numerator':J_ASecond_Numerator, 'J_ANumerator':J_ANumerator, 'J_BFirst_Numerator':J_BFirst_Numerator,\
'J_BSecond_Numerator':J_BSecond_Numerator, 'J_BNumerator':J_BNumerator, 'J_CFirst_Numerator':J_CFirst_Numerator,\
'J_CSecond_Numerator':J_CSecond_Numerator, 'J_CNumerator':J_CNumerator})

#print("J_A = " , equation[0])
#np.savetxt("foo.csv", Et, delimiter=",")
###### output as csv file or table as well
column_order = ['time', 'A_int(mM)', 'A_ext(mM)','B_int(mM)','B_ext(mM)','C_int(mM)', 'C_ext(mM)',\
'J_AMN', 'J_BMN','J_CMN','Et','alphaM','betaM',
'alphaN','betaN','gammaM','gammaN','R_N','R_MM','R_M','R_NN',\
'JDenominator','J_AFirst_Numerator','J_ASecond_Numerator','J_ANumerator','J_BFirst_Numerator','J_BSecond_Numerator','J_BNumerator',\
'J_CFirst_Numerator','J_CSecond_Numerator','J_CNumerator']
df.to_csv('unordered' + '.csv', sep = ',')

s3 = input('Do you want to save data in a separate csv file?'+'\n'+' if yes, enter YES ')
if (str(s3) == 'YES'):
	out3 = input('Choose filename: (do not include extension) ')
#	df.to_csv(str(out3) + '.csv', sep = ',')
	df[column_order].to_csv(str(out3) + '_AntiCompete' +'.csv', sep = ',',index=False)
	print('The new csv file has been saved under the name of ' + str(out3) + '_AntiCompete' +'.csv')
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
    plt.savefig(str(out6) + '_AntiCompete' + '.png')
    plt.savefig(str(out6) + '_AntiCompete',format = 'eps', dpi=1000, bbox_inches='tight')
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










