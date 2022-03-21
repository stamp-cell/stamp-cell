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
NumberofEquations = 19
NumberofICs = 6
NumberofParameters = 14

print ('\n' + "Please enter number for following inputs:" + '\n') 
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


g_anti = input('Enter g_antiporter (per_s) , (E translocation rate constant from ext to int) ')
E_t = input('Enter [E]_t (mM) ')
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
	legend_ICs[0] = "[A]_M (A_int concentrations (mM))"
	legend_ICs[1] = "[A]_N (A_ext concentrations (mM))"
	legend_ICs[2] = "[B]_M (B_int concentrations (mM))"
	legend_ICs[3] = "[B]_N (B_ext concentrations (mM))"
	legend_ICs[4] = "[C]_M (C_int concentrations (mM))"
	legend_ICs[5] = "[C]_N (C_ext concentrations (mM))"
	legend_parameters[0] = "K_A" 
	legend_parameters[1] = "K_B" 
	legend_parameters[2] = "K_C"
	legend_parameters[3] = "A stoichiometric"
	legend_parameters[4] = "B stoichiometric"	
	legend_parameters[5] = "C stoichiometric"
	legend_parameters[6] = "g_anti"
	legend_parameters[7] =" E_t_Tmax (mM)"
	legend_equation[0] = "p_anti" # E_t*g_anti
	legend_equation[1] = "alpha_M" # (dimensionless)
	legend_equation[2]= "beta_M" # (dimensionless)
	legend_equation[3] = "gamma_M" # (dimensionless)
	legend_equation[4] = "alpha_N" # (dimensionless)
	legend_equation[5]=  "beta_N" # (dimensionless)
	legend_equation[6] = "gamma_N" # (dimensionless)
	legend_equation[7] = "R_N" # 1 + (alpha_N)^a + (beta_N)^b + (gamma_N)^c
	legend_equation[8] = "R_MM" # (alpha_M)^a + (beta_M)^b + (gamma_M)^c
	legend_equation[9] = "R_M" # 1 + (alpha_M)^a + (beta_M)^b + (gamma_M)^c
	legend_equation[10] = "R_NN" # (alpha_N)^a + (beta_N)^b + (gamma_N)^c
	legend_equation[11] = "J Denominator" #R_M*R_NN + R_N*R_MM
	legend_equation[12] = "J_anti_First_Numerator"
	legend_equation[13] = "J_anti_Second_Numerator"
	legend_equation[14] = "J_anti_Numerator"
	legend_equation[15] = "J_antiporter(M, N) antiporter influx in component AE1 (mM_per_s)"
	legend_equation[16] = "J_A(M, N) A_influx in component AE1 (mM_per_s)"
	legend_equation[17] = "J_B(M, N) B_influx in component AE1 (mM_per_s)"
	legend_equation[18] = "J_C(M, N) C_influx in component AE1 (mM_per_s)"

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
	parameters[0] = float(K_A) 
	parameters[1] = float(K_B) 
	parameters[2] = float(K_C)
	parameters[3] = float(a)
	parameters[4] = float(b)	
	parameters[5] = float(c)
	parameters[6] = float(g_anti)
	parameters[7] = float(E_t)	#E_Tmax (mM)
	parameters[8] = 0.00000
	parameters[9] = 1.00000
	parameters[10] = 0.00000
	parameters[11] = 0.00000
	parameters[12] = 0.00000
	parameters[13] = 0.00000
	return (ICs, parameters)

def computeRates(time, ICs, parameters):
    rates = [0.0] * NumberofICs; equation = [0.0] * NumberofEquations
    rates[0] = parameters[8] #d/dt A_M in component concentrations (mM) int=M
    rates[1] = parameters[9] #d/dt A_N in component concentrations (mM) ext=N
    rates[2] = parameters[10] #d/dt B_M in component concentrations (mM)
    rates[3] = parameters[11] #d/dt B_N in component concentrations (mM)
    rates[4] = parameters[12] #d/dt C_M in component concentrations (mM)
    rates[5] = parameters[13] #d/dt C_N in component concentrations (mM)
    return(rates)

def computeEquation(parameters, ICs, time):
	equation = array([[0.0] * len(time)] * NumberofEquations)
	ICs = array(ICs)
	time = array(time)
# 0) p_anti: E_t*g_anti
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
# 7) R_N = 1 + (alpha_N)^a + (beta_N)^b + (gamma_N)^c
	equation[7] = (1.000000 + power(equation[4], parameters[3]) + power(equation[5], parameters[4]) + power(equation[6], parameters[5]))
# 8) R_MM = (alpha_M)^a + (beta_M)^b + (gamma_M)^c
	equation[8] = (power(equation[1], parameters[3]) + power(equation[2], parameters[4]) + power(equation[3], parameters[5]))
# 9) R_M = 1 + (alpha_M)^a + (beta_M)^b + (gamma_M)^c
	equation[9] = (1.000000 + power(equation[1], parameters[3]) + power(equation[2], parameters[4]) + power(equation[3], parameters[5]))
# 10) R_NN = (alpha_N)^a + (beta_N)^b + (gamma_N)^c
	equation[10] = (power(equation[4], parameters[3]) + power(equation[5], parameters[4]) + power(equation[6], parameters[5]))
# 11) J Denominator: R_M*R_NN + R_N*R_MM
	equation[11] = (equation[9]*equation[10]+equation[7]*equation[8])
# 12) J_symp First_Numerator: = alpha_M^a * beta_N^b * gamma_N^c
	equation[12] = ((power(equation[1], parameters[3]))*(power(equation[5], parameters[4]))*(power(equation[6], parameters[5])))
# 13) J_symp Second_Numerator: = alpha_N^a * beta_M^b * gamma_M^c
	equation[13] = ((power(equation[4], parameters[3]))*(power(equation[2], parameters[4]))*(power(equation[3], parameters[5])))
# 14) J_symp Numerator:
	equation[14] = (equation[12]-equation[13])
# 15) J_symp(M, N):
	equation[15] = equation[0] * (equation[14]/equation[11])
# 16) J_A(M, N): = a*J_symp
	equation[16] = parameters[3]*equation[15]
# 17) J_B(M, N): b*J_symp
	equation[17] = parameters[4]*equation[15]
# 18) J_C(M, N): = c*J_symp
	equation[18] = parameters[5]*equation[15]
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
p_anti = pd.Series(equation[0]) # E_t*g_anti
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
J_anti_First_Numerator = pd.Series(equation[12])
J_anti_Second_Numerator = pd.Series(equation[13])
J_anti_Numerator = pd.Series(equation[14])
J_antiporterMN = pd.Series(equation[15])
J_AMN = pd.Series(equation[16])
J_BMN = pd.Series(equation[17])
J_CMN = pd.Series(equation[18])


df = pd.DataFrame({'time':time, 'A_int(mM)':AM, 'A_ext(mM)':AN, 'B_int(mM)':BM,'B_ext(mM)' :BN,\
'C_int(mM)':CM, 'C_ext(mM)':CN, 'p_anti':p_anti, 'J_antiporterMN':J_antiporterMN,\
'J_AMN':J_AMN, 'J_BMN':J_BMN,'J_CMN':J_CMN,'alphaM':alphaM,'alphaN':alphaN,'betaM':betaM,'betaN':betaN,'gammaM':gammaM,'gammaN':gammaN,\
'R_N':RN,'R_MM':RMM , 'R_M':RM, 'R_NN':RNN, 'JDenominator':JDenominator, 'J_anti_First_Numerator':J_anti_First_Numerator,\
'J_anti_Second_Numerator':J_anti_Second_Numerator, 'J_anti_Numerator':J_anti_Numerator})

#print("J_A = " , equation[0])
#np.savetxt("foo.csv", Et, delimiter=",")
###### output as csv file or table as well
column_order = ['time', 'A_int(mM)', 'A_ext(mM)','B_int(mM)','B_ext(mM)','C_int(mM)', 'C_ext(mM)',\
'J_antiporterMN','J_AMN', 'J_BMN','J_CMN','p_anti','alphaM','betaM',
'alphaN','betaN','gammaM','gammaN','R_N','R_MM','R_M','R_NN',\
'JDenominator','J_anti_First_Numerator','J_anti_Second_Numerator','J_anti_Numerator']
df.to_csv('unordered' + '.csv', sep = ',')

s3 = input('Do you want to save data in a separate csv file?'+'\n'+' if yes, enter YES ')
if (str(s3) == 'YES'):
	out3 = input('Choose filename: (do not include extension) ')
#	df.to_csv(str(out3) + '.csv', sep = ',')
	df[column_order].to_csv(str(out3) + '_AntiSimple' +'.csv', sep = ',',index=False)
	print('The new csv file has been saved under the name of ' + str(out3) + '_AntiSimple' +'.csv')
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
    plt.savefig(str(out6) + '_AntiSimple' + '.png')
    plt.savefig(str(out6) + '_AntiSimple',format = 'eps', dpi=1000, bbox_inches='tight')
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



