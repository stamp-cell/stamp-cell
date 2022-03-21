from __future__ import division
import numpy as np
import math 
import pylab
import pandas as pd
import csv
import scipy as sp
import pylab as plt
from scipy.integrate import odeint
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
import subprocess
import sys
import re
#plt.rc('text', usetex=True)
############################################ Constants: ###################################################

T = 310.0
R = 8314 # j/mol/k #R = 1.987 # calories/mol/K
F = 96484 # coulombs/mol




#____________________________________________________
##### Or

# Ions:
ions = ["Ca", "Na", "K"]

z_i = {
'Ca' : 2 , 'Na': 1, 'K': 1
 }

# print z_i['Ca']
for i in ions:
	print(z_i[i])

P = {
'Ca' : 5.4*(10**(-4)) , 'Na': 6.75*(10**(-7)), 'K': 1.93*(10**(-7))
 } #cm/s

print(P['Ca'])
for i in ions:
	print (P[i])


gamma_i_in = {
'Ca' : 1.0 , 'Na': 0.75, 'K': 0.75
 } #cm/s

print(gamma_i_in['Ca'])
for i in ions:
	print(gamma_i_in[i])


gamma_i_out = {
'Ca' : 0.341 , 'Na': 0.75, 'K': 0.75
 } #cm/s

print(gamma_i_out['Ca'])
for i in ions:
	print(gamma_i_out[i])


c_in = { 'Ca_rest':10.0, 'Na':10.0, 'K':145.0} #mmol/liter
c_o = { 'Ca':140.0, 'Na':140.0, 'K':5.4} #mmol/liter
for i in ions:
	print(c_o[i])

#______________________________________________________
# Cell Geomwtry:

# Dimensions:  
l = 100	# length micro meter
r = 11.0	# radius in micro meter
V_cell_m = (math.pi)*(r**2)*l # CellVolume-micro meter cubed
V_cell = V_cell_m*(10**(-9)) # CellVolume-micro liter
# print V_cell
AGeo = 2.0*(math.pi)*(r**2)+2.0*(math.pi)*r*l # GeometricMembraneArea micro meter^2
AGeo_cm = AGeo*1e-8 # cm2
# print AGeo_cm
R_CG = 2.0 # R_CG=ratio of the A_cap to A_Geo
ACap = R_CG*AGeo_cm # Capacitive membrane area 
# print ACap 
V_myo = V_cell*0.68 # Myoplasm Volume microliter
V_mit = V_cell*0.26 # Mitochondria Volume microliter
V_SR = V_cell*0.06 # SR volume microliter
V_NSR = V_cell*0.0552 # A network SR (NSR) volume microliter
V_JSR = V_cell*0.0048 # A junctional SR (JSR) volume microliter
V_cleft = (V_cell/0.88)*0.12 # Cleft Volume microliter

#________________________________________________________

class ChannelsTransporters():#Variables:
	"""Full Mathematical Frame Work implemented in Python"""
	# Cell Geomwtry:
	# Dimensions:  
	l = 100.00	# length micro meter
	r = 11.0	# radius in micro meter
	#V_cell_m = (math.pi)*(r**2)*l # CellVolume-micro meter cubed
	#V_cell = V_cell_m*(10**(-9)) # CellVolume-micro liter
	# print V_cell
	#AGeo = 2.0*(math.pi)*(r**2)+2.0*(math.pi)*r*l # GeometricMembraneArea micro meter^2
	#AGeo_cm = AGeo*1e-8 # cm2
	# print AGeo_cm
	R_CG = 2.0 # R_CG=ratio of the A_cap to A_Geo
	#ACap = R_CG*AGeo_cm # Capacitive membrane area 
	# print ACap 
	#V_myo = V_cell*0.68 # Myoplasm Volume microliter
	#V_mit = V_cell*0.26 # Mitochondria Volume microliter
	#V_SR = V_cell*0.06 # SR volume microliter
	#V_NSR = V_cell*0.0552 # A network SR (NSR) volume microliter
	#V_JSR = V_cell*0.0048 # A junctional SR (JSR) volume microliter
	#V_cleft = (V_cell/0.88)*0.12 # Cleft Volume microliter

	# Standard Ionic Concentrations:
	K_o = 5.4 #float(K_o_input) 
	K_i = 145.0 # float(K_i_input)
	Na_o = 140.0 # float(Na_o_input)
	Na_i = 10.0 # float(Na_i_input)
	Ca_o = 140.0 # float(Ca_o_input)
	Ca_i_rest = 10.0 # float(Ca_i_input) # 0.00017948816388306 mmolar CellMl

	# Library:
	z_iCa = 2.0
	z_iNa = 1.0
	z_iK = 1.0

	P_Ca = 5.4*(10**(-4))
	P_Na = 6.75*(10**(-7))
	P_K = 1.93*(10**(-7)) #cm/s

	gamma_i_in_Ca = 1.0
	gamma_i_in_Na = 0.75
	gamma_i_in_K = 0.75
	
	gamma_i_out_Ca = 0.341
	gamma_i_out_Na = 0.75
	gamma_i_out_K = 0.75
 
	# Constants:
	T = 310.0
	R = 8314 # j/mol/k # 1.987 calories/mol/K
	F = 96484 # coulombs/mol


	# Half saturation concentration of ion channels micromol/L
	KmCa = 0.6 # float(KmCa_input) # half-saturation concentration of channel Ca, micromol/L 

	############ Time-dependet K current I_K:
	P_NaK = 0.01833


	############ Plateau K+ current: I_Kp
	g_Kp = 0.0183 #millisiemens/,uF;


	############ NaCa Exchanger
	K_NaCa = 2000.0 # microA/microF
	K_mNa = 87.5 # mmol/L
	K_mCa = 1.38 # mmol/L
	k_sat = 0.1 
	eta = 0.35

	############ NaK ATPase
	I_bar_NaK = 1.5 # microA/microF
	K_mNai = 10.0 # mmol/liter
	K_mKo = 1.5 # mmol/liter
	
	############ NonSpecific Ca2+ activated Current (I_nsCa):
	P_nsCa = 1.75*(10**(-7)) # cm/s
	K_mnsCa = 1.2 # micromol/liter
	############ Sacrolemmal Ca2+ pump (I_pCa):
	K_mpCa  = 0.5 # micromol/liter
	I_bar_pCa = 1.15 # microA/microF
	
	############ Ca background leakage current (I_Ca_b):
	g_Cab = 0.003016 # millisiemens/microF
	

	# Na background leakage current (I_Na_ohm):
	g_Nab = 0.00141 # millisiemens/microF

	#  Troponin (TRPN) and calmodulin (CMDN);
	TRPN_ave = 70.0 # micromol/L
	CMDN_ave = 50.0 # micromol/L
	KmTRPN = 0.5 # micromol/L
	KmCMDN = 2.38 # micromol/L

	# Ca2+induced Ca2+ release of JSR
	#situ = str(raw_input('VoltageClamp or ActionPotential '))
	#g_re1_bar = 18.0 # ms-1 for voltage clamp situation and =60 for action potential simulation
	deltaCa_ith = 0.18 # (micrmol/liter) threshold for external triggering of Ca2+ release from JSR, gmol/L
	K_mre1 = 0.8 # micromol/liter
	tau_on = 2.0 # milliseconds
	tau_off = 2.0 # milliseconds



	# Ca2+ release of JSR undet Ca overload conditions
	tau_on = 2.0 # milliseconds
	tau_off = 2.0 # milliseconds
	g_re1_bar = 4.0 # ms-1 
	CSQN = 0.7 # or higher

	# Ca buffer in JSR and CSQN
	CSQN_ave = 10 # mmol/liter
	K_mCSQN = 0.8 # mmol/liter


	# Ca uptake and leakage of NSR: I_up and I_leak
	I_up_ave = 0.005 # mmol/liter per milliseond
	K_mup = 0.92 # micromol/liter
	Ca_NSR_ave = 15.0 # mmol/liter

	# Transduction of Ca ions from NSR to JSR: Itr
	tau_tr =  180.0 # milliseconds
	"""Full EP17 Model implemented in Python"""
	C_m  =   1.0 #membrane capacitance, in uF/cm^2
	g_Na = 16.0 # maximum conductance of Na Channel (millisiemens/microF)
	#g_Na = 120.0 # Sodium (Na) maximum conductances, in mS/cm^2"""

	#g_K  =  36.0
	#"""Postassium (K) maximum conductances, in mS/cm^2"""
	g_Ki  = 0.75
	g_L  =   0.3
	"""Leak maximum conductances, in mS/cm^2"""

	E_Na  = (R*T/F)*np.log(Na_o/Na_i) # ENa, from 54.4 to 70 mV)
	#print E_Na
	#E_Na =  50.0
	#"""Sodium (Na) Nernst reversal potentials, in mV"""

	E_K  = -77.0
	"""Postassium (K) Nernst reversal potentials, in mV"""

	E_L  = -54.387
	"""Leak Nernst reversal potentials, in mV"""

	t = sp.arange(0.0, 600.0, 0.01)
	""" The time to integrate over """

	#______________________________________________________
	# Cell Geomwtry:

	# Dimensions:  
	#l = 100.00	# length micro meter
	#r = 11.0	# radius in micro meter
	def V_cell_m (self, r,l):
		return (math.pi)*(self.r**2)*self.l # CellVolume-micro meter cubed
	
	def V_cell (self, r,l):
		return self.V_cell_m(r,l)*(10**(-9)) # CellVolume-micro liter
	#print V_cell
	def AGeo (self,r,l):
		return 2.0*(math.pi)*(self.r**2)+2.0*(math.pi)*self.r*self.l # GeometricMembraneArea micro meter^2
	def AGeo_cm (self,r,l):
		return self.AGeo(r,l)*1e-8 # cm2
	# print AGeo_cm
	# R_CG = 2.0 # R_CG=ratio of the A_cap to A_Geo
	def ACap(self, r,l): 
		return self.R_CG*self.AGeo_cm(r,l) # Capacitive membrane area 
	# print ACap 
	
	def V_myo (self,r,l):
		return self.V_cell(r,l)*0.68 # Myoplasm Volume microliter

	def V_mit(self, r,l):
		return self.V_cell(r,l)*0.26 # Mitochondria Volume microliter

	def V_SR(self,r, l):
		return self.V_cell(r,l)*0.06 # SR volume microliter

	def V_NSR(self, r,l):
		self.V_cell(r,l)*0.0552 # A network SR (NSR) volume microliter

	def V_JSR(self, r,l):
		return self.V_cell(r,l)*0.0048 # A junctional SR (JSR) volume microliter

	def V_cleft(self, r,l):
		return (self.V_cell(r,l)/0.88)*0.12 # Cleft Volume microliter
################################################	Channels	################################################
	###########################	Calcium Dependent Chloride Channels (CaCC)
	def I_bar_Cacc (self, V):
		IbarCacc = self.P_Cl*((V*(self.F**2))/(self.R*self.T))* \
(self.Cl_o*sp.exp(-V*self.F/(self.R*self.T))-self.Cl_i)/(sp.exp(-V*self.F/(self.R*self.T))-1.0)
		return IbarCaCacc
	def f_Cl_Cacc_arreola (self, Ca_i): # Changed by SZ from f_Ca to f_Ca_l
		return 1.0/(1.0 + self.K2_CaccArreola(V)*((Ca_i/self.K2_CaccArreola(V))**2)+(Ca_i/self.K2_CaccArreola(V))+1) # Ca-dependent inactivation gate of the L-type Ca2+ channel

	def K2_CaccArreola(self, V):
		return self.k1*sp.exp((-1*self.gamma1*self.F)/(self.R*self.T))

	def K1_CaccArreola(self, V):
		return self.k2*sp.exp((-1*self.gamma2*self.F)/(self.R*self.T))

	def I_Cl_Cacc (self, V, d, f, Ca_i):
		ICa = d*f*self.f_Ca_l(Ca_i)*self.I_bar_Ca(V)  # f*f_Ca = inactivation process
		return ICa


	###########################	Sodium Channels
	############ Voltage Gated Sodium Channels (Na_v) Fast Sodium current (i_Na_v):
	# T = 310.0
	# R = 1.987 # calories/mol/K
	# F = 96484 # coulombs/mol
	# g_Na = 16.0 # maximum conductance of Na Channel (millisiemens/microF)
	def I_Na_v(self, V, m, h, j): # SZ 
		"""
		Membrane current (in uA/cm^2)
		Sodium (Na = element name)

		|  :param V:
		|  :param m:
		|  :param h:
		|  :param j:
		|  :return:
		"""
		return self.g_Na * m**3 * h * j * (V - self.E_Na) # m: activation gate, h: fast inactivation gate, j: slow inactivation gate


	def alpha_j(self, V):
		" fast_sodium_current_j_gate (Slow Inactivation) (per_millisecond)"
		if V < -40.0:
			alphaj = ((-127140.0*sp.exp(0.244400*V)-3.47400e-05*sp.exp(-0.0439100*V))*(V+37.7800))/(1.00000+sp.exp(0.311000*(V+79.2300)))
		else:
			alphaj = 0.0
		return alphaj


	def beta_j(self, V):
		" fast_sodium_current_j_gate (Slow Inactivation) (per_millisecond)"
		if V < -40.0:
			betaj = (0.121200*sp.exp(-0.0105200*V))/(1.00000+sp.exp(-0.137800*(V+40.1400)))	
		else: 
			betaj = (0.300000*sp.exp(-2.53500e-07*V))/(1.00000+sp.exp(-0.100000*(V+32.0000)))
		return  betaj


	def alpha_m(self, V):
		"""Channel gating kinetics. Functions of membrane voltage"""
		return 0.32*(V+47.13)/(1.0 -sp.exp(-0.1*(V+47.13))) # modified from ref 44


	def beta_m(self, V):
		"""Channel gating kinetics. Functions of membrane voltage"""
		return 0.08*sp.exp(-V/11.0)


	def alpha_h(self,V):
		" fast_sodium_current_j_gate (Slow Inactivation) (per_millisecond)"
		if V < -40.0:
			alphah = 0.135*sp.exp((V+8.0)/-6.8)
		else:
			alphah = 0.0
		return alphah


	def beta_h(self,V):
		" fast_sodium_current_j_gate (Slow Inactivation) (per_millisecond)"
		if V < -40.0:
			betah = 3.56*sp.exp(0.079*V) + 3.1*(10**5)*sp.exp(0.35*V)
		else:
			betah = 1.0/(0.13*(1.0 + sp.exp((V+10.66)/-11.1)))
		return betah

	def alpha_n(self, V):
		"""Channel gating kinetics. Functions of membrane voltage"""
		return 0.01*(V+55.0)/(1.0 - sp.exp(-(V+55.0) / 10.0))

	def beta_n(self, V):
		"""Channel gating kinetics. Functions of membrane voltage"""
		return 0.125*sp.exp(-(V+65) / 80.0)
	

	############ Sodium Channel Regular, Na background current (i_Na_ohm):
	E_NaN = float(E_Na)
		# g_Nab = 0.00141 # millisiemens/microF
	def I_Na_ohm (self, V):
		return self.g_Nab*(V-self.E_NaN)

	def I_Na_ghk (self, V):
		IbarCaNa = self.P_Na*((V*(self.F**2))/(self.R*self.T))* \
(self.Na_i*sp.exp(-V*self.F/(self.R*self.T))-self.Na_o)/(sp.exp(V*self.F/(self.R*self.T))-1.0)
		return IbarCaNa
	################################################################ 
	###########################	Calcium Channels
	############ Calcium SOC:
	def f_Ca_soc (self, Ca_i):
		return 1.0/(1.0 + (Ca_sr/self.KmCa)**eta_soc) # Ca-dependent inactivation gate of the L-type Ca2+ channel
	
	############ Calcium Ltype:

	def f_Ca_l (self, Ca_i): # Changed by SZ from f_Ca to f_Ca_l
		return 1.0/(1.0 + (Ca_i/self.KmCa)**2) # Ca-dependent inactivation gate of the L-type Ca2+ channel

	def d_inft (self, V):
	#activation gate and (voltage-dependent) of the L-type Ca channel
		return 1.0/(1.0+(sp.exp(-(V+10.0))/6.24)) 
	def tau_d (self, V):
		return self.d_inft(V)*(1.0 - sp.exp(-(V+ 10.0)/6.24))/(0.035*(V+ 10.0))


	def f_inft (self, V):
		return 1.0/(1.0 + sp.exp((V+35.06)/8.6))+0.6/(1 + sp.exp((50.0-V)/20.0))


	def tau_f (self, V): 
		return 1.0/( 0.0197*sp.exp(-(0.0337*(V+ 10.0))**2.0)+ 0.02)


	def alpha_d (self,V):
		return self.d_inft(V)/self.tau_d(V)
		

	def beta_d (self,V):
		return (1.0-self.d_inft(V))/self.tau_d(V)


	def alpha_f (self,V):
		return self.f_inft(V)/self.tau_f(V)


	def beta_f (self,V):
		return (1.0-self.f_inft(V))/self.tau_f(V)

	def E_si(self, Ca_i):
		return 7.70000-13.0287*sp.log(Ca_i/1.00000)

	def I_Ca_l(self, V, d, f, Ca_i):
		return 0.0900000*d*f*(V-self.E_si(Ca_i))

	#def I_si(self, V, d, f, Ca_i): -> SZ: I_si chenged to I_Ca_l
	#	return 0.0900000*d*f*(V-self.E_si(Ca_i))

	############ Calcium Ltype Permeable to Sodium and Potassium:
	def I_bar_Ca (self, V): 
		IbarCa = self.P_Ca*self.z_iCa*((V*(self.F**2))/(self.R*self.T))* \
(self.gamma_i_in_Ca*self.Ca_i_rest*sp.exp(self.z_iCa*V*self.F/(self.R*self.T)) - self.gamma_i_out_Ca*self.Ca_o)/(sp.exp((self.z_iCa*V*self.F)/(self.R*self.T))-1.0)
		return IbarCa

	def I_bar_CaNa (self, V):
		IbarCaNa = self.P_Na*self.z_iNa*((V*(self.F**2))/(self.R*self.T))* \
(self.gamma_i_in_Na*self.Na_i*sp.exp(self.z_iNa*V*self.F/(self.R*self.T))-self.gamma_i_out_Na*self.Na_o)/(sp.exp(self.z_iNa*V*self.F/(self.R*self.T))-1.0)
		return IbarCaNa

	def I_bar_CaK (self, V):
		IbarCaK = self.P_K*self.z_iK*((V*(self.F**2))/(self.R*self.T))* \
(self.gamma_i_in_K*self.K_i*sp.exp(self.z_iK*V*self.F/(self.R*self.T))-self.gamma_i_out_K*self.K_o)/(sp.exp(self.z_iK*V*self.F/(self.R*self.T))-1.0)
		return IbarCaK



	def I_Ca_lcanak (self, V, d, f, Ca_i):
		ICa = d*f*self.f_Ca_l(Ca_i)*self.I_bar_Ca(V)  # f*f_Ca = inactivation process
		return ICa
	#def I_Ca (self, V, d, f, Ca_i):
	#	ICa = d*f*self.f_Ca(Ca_i)*self.I_bar_Ca(V)  # f*f_Ca = inactivation process
	#	return ICa


	def I_CaK (self, V, d, f, Ca_i):
		ICaK = d*f*self.f_Ca_l(Ca_i)*self.I_bar_CaK(V)
		return ICaK

	def I_CaNa (self, V, d, f, Ca_i):
		ICaNa = d*f*self.f_Ca_l(Ca_i)*self.I_bar_CaNa(V)
		return ICaNa

	def I_Cal_total (self, V, d, f, Ca_i):
		ICat = self.I_Ca_lcanak(V, d, f, Ca_i) + self.I_CaK (V, d, f, Ca_i) + self.I_CaNa (V, d, f, Ca_i) 
		return ICat

	############ Calcium T type:
	def I_Ca_t(self, V, m, h, j): # SZ 
		"""
		Membrane current (in uA/cm^2)
		Sodium (Na = element name)

		|  :param V:
		|  :param m:
		|  :param h:
		|  :param j:
		|  :return:
		"""
		return self.g_Ca_t * m**2 * h * \
	((self.z_iCa**2)*V*(self.F**2)/(self.R*self.T))*(self.Ca_o*sp.exp(-self.z_iCa*V*self.F/(self.R*self.T))-self.Ca_i)/(sp.exp(-self.z_iCa*V*self.F/(self.R*self.T))-1.0)
 # m: activation gate, h: fast inactivation gate, j: slow inactivation gate

############ Time-dependet K current I_K:
# P_NaK = 0.01833

	g_K = 0.282*math.sqrt(K_o/5.4) # maximum conductance of K Channel (millisiemens/microF) (K_o dependense condutance)
	E_K  = (R*T/F)*sp.log((K_o + P_NaK*Na_o)/(K_i + P_NaK*Na_i)) 
	#def X_i (self, V):
	#	return 1.0/(1.0 + sp.exp((V-56.26)/32.1))

	def I_K (self, V, x):
		return self.g_K*(1.0/(1.0 + sp.exp((V-56.26)/32.1)))*(x**2)*(V-self.E_K)

	def alpha_x (self, V):
		return 7.19*(10.0**(-5))*(V+30)/(1.0-sp.exp(-0.148* (V+30)))
	
	def beta_x (self, V):
		return 1.31 * (10**(-4))*(V+30)/(-1.0+sp.exp(0.0687*(V+30)))



############ Time-independet K current I_Ki:
	#g_Ki = 0.75*math.sqrt(K_o/5.4) # g_bar_K = 0.75, maximum conductance of K Channel (millisiemens/microF)
	E_Ki  = (R*T/F)*sp.log((K_o)/(K_i)) 

	def K_inf (self, V):
	    return self.alpha_Ki(V)/(self.alpha_Ki(V)+self.beta_Ki(V))

	def I_Ki (self, V):
		return self.g_Ki*self.K_inf(V)*(V-self.E_Ki)

	def alpha_Ki(self, V):
		return 1.02/(1.0 + sp.exp(0.2385*(V- self.E_Ki-59.215)))

	def beta_Ki(self, V): 
		return (0.49124*sp.exp(0.08032*(V-self.E_Ki+5.476))+ sp.exp(0.06175*(V- self.E_Ki-594.31)))/(1+sp.exp(-0.5143*(V-self.E_Ki +4.753)))


############ Plateau K+ current: I_Kp
# g_Kp = 0.0183 #millisiemens/,uF;
	E_Kp = (R*T/F)*np.log((K_o)/(K_i))
	def Kp(self, V):
		return 1.0/(1.0 +sp.exp((7.488-V)/5.98))

	def I_Kp(self, V):
		return self.g_Kp*self.Kp(V)*(V-self.E_Kp)
	
	def I_pCa(self, Ca_i):
		return self.I_bar_pCa*(Ca_i/(self.K_mpCa + Ca_i))

################################################################ NaCa Exchanger i_NaCa: 

#K_NaCa = 2000.0 # microA/microF
#K_mNa = 87.5 # mmol/L
#K_mCa = 1.38 # mmol/L
#k_sat = 0.1
#eta = 0.35

	def I_NaCa(self, V): #, Na_o, Ca_o, Na_i, Na_o):
		return self.K_NaCa*(1.0/((self.K_mNa)**3+(self.Na_o)**3))*(1.0/(self.K_mCa+(self.Ca_o)**3))*(1.0/(1.0+(self.k_sat*sp.exp((self.eta-1.0)*V*(self.F/(self.R*self.T)))))*((sp.exp(self.eta*V*(self.F/(self.R*self.T)))*(self.Na_i**3)*self.Ca_o) - (sp.exp((self.eta-1.0)*V*(self.F/(self.R*self.T)))*(self.Na_o**3)*self.Ca_i_rest)))

################################################################ NaK Pump i_NaK:

# I_bar_NaK = 1.5 # microA/microF
# K_mNai = 10.0 # mmol/liter
# K_mKo = 1.5 # mmol/liter
	#def sigma (self, Na_o):
	#sigma = (1.0/7.0)*(sp.exp(Na_o/67.3)-1.0)
	
	def f_NaK (self, V):
		return 1.0/(1.0+0.1245*sp.exp(-0.1*((V*self.F)/(self.R*self.T)))+0.0365*((1.0/7.0)*(sp.exp(self.Na_o/67.3)-1.0))*sp.exp(-(V*self.F)/(self.R*self.T)))

	def I_NaK (self, V):
		return self.I_bar_NaK*self.f_NaK(V)*(1.0/(1.0+(self.K_mNai/self.Na_i)**(1.5)))*(self.K_o/(self.K_o+ self.K_mKo))

################################################################ NonSpecific Ca2+ activated Current (I_nsCa):
# P_nsCa = 1.75*(10**(-7)) # cm/s
# K_mnsCa = 1.2 # micromol/liter
	E_nsCa = ((R*T)/F)*sp.log((K_o + Na_o)/(K_i+Na_i))

	def I_bar_nsNa (self, V):
		IbarnsNa = self.P_nsCa*self.z_iNa*((V*(self.F**2))/(self.R*self.T))* \
(self.gamma_i_in_Na*self.Na_i*sp.exp(self.z_iNa*V*self.F/(self.R*self.T))-self.gamma_i_out_Na*self.Na_o)/(sp.exp(self.z_iNa*V*self.F/(self.R*self.T))-1.0)
		return IbarnsNa

	def I_bar_nsK (self, V):
		IbarnsK = self.P_nsCa*self.z_iK*((V*(self.F**2))/(self.R*self.T))* \
(self.gamma_i_in_K*self.K_i*sp.exp(self.z_iK*V*self.F/(self.R*self.T))-self.gamma_i_out_K*self.K_o)/(sp.exp(self.z_iK*V*self.F/(self.R*self.T))-1.0)
		return IbarnsK

	def I_nsK (self,V,Ca_i):
		return self.I_bar_nsK(V)*(1.0/(1.0 + (self.K_mnsCa/Ca_i))**3)

	def I_nsNa (self,V,Ca_i):
		return self.I_bar_nsNa(V)*(1.0/(1.0 + (self.K_mnsCa/Ca_i))**3)

	def I_nsCa (self, V, Ca_i):
		return self.I_nsK(V, Ca_i) + self.I_nsNa(V, Ca_i)

################################################################ Sacrolemmal Ca2+ pump (I_pCa):

#K_mpCa  = 0.5 # micromol/liter
#I_bar_pCa = 1.15 # microA/microF
	def I_pCa (self, Ca_i):
		return self.I_bar_pCa*(Ca_i/(self.K_mpCa + Ca_i))

################################################################ Ca background current (I_Ca_b):

# g_Cab = 0.003016 # millisiemens/microF

	E_CaN = ((R*T)/(2.0*F))*sp.log(Ca_o/Ca_i_rest)
	def I_Ca_b (self, V):
		return self.g_Cab*(V-self.E_CaN)

################################################################ Total time independent current (I_v):
	#def I_v (self, I_Ki, I_Kp, I_pCa, I_Na_ohm, I_Ca_b, I_NaK):
	def I_v (self, V, Ca_i):
		return self.I_Ki(V) + self.I_Kp(V) + self.I_pCa(Ca_i) + self.I_Na_ohm(V) + self.I_Ca_b(V) + self.I_NaK(V)#, Na_o, Nai, K_o)

# Total Membrane current (I_i):

	def I_i (self, V, d, f, Ca_i, x, m, h, j): #, I_Ca_b, I_Ca, I_NaCa, I_pCa, I_nsCa, I_Kp, I_Ki, I_K, I_NaK, I_Na, I_Na_ohm):
		return self.I_Ca_b(V) + self.I_Ca_lcanak(V, d, f, Ca_i)  + self.I_pCa(Ca_i) + self.I_nsCa(V, Ca_i) + self.I_Ki(V) + \
self.I_K(V, x) + self.I_NaK(V) +self.I_Na_v(V, m, h, j) + self.I_Na_ohm(V) + self.I_NaCa (V)



	#	"""
	#	Membrane current (in uA/cm^2)
	#	Potassium (K = element name)
	#
	#	|  :param V:
	#	|  :param h:
	#	|  :return:
	#	"""
	#	return self.g_K  * n**4 * (V - self.E_K)
    #  Leak
	def I_L(self, V):
		"""
		Membrane current (in uA/cm^2)
		Leak

		|  :param V:
		|  :param h:
		|  :return:
		"""
		return self.g_L * (V - self.E_L)

	def I_inj(self, t):
		"""
		External Current

		|  :param t: time
		|  :return: step up to 10 uA/cm^2 at t>100
		|           step down to 0 uA/cm^2 at t>200
		|           step up to 35 uA/cm^2 at t>300
		|           step down to 0 uA/cm^2 at t>400
		"""
#param period=500, iStim_mag=8.0, iStim_beg=20.0,  iStim_dur=5.0

		return 10*(t>100) - 10*(t>200) + 35*(t>300) - 35*(t>400)

		#return 10*(t>100) - 10*(t>105) + 25*(t>625) - 25*(t>1025)
	@staticmethod
	def dALLdt(X, t, self):
		"""
		Integrate

		|  :param X:
		|  :param t:
		|  :return: calculate membrane potential & activation variables
		"""
		V, m, h, j, n , d, f, Ca_i, x, Na_i = X
		#dVdt = (self.I_inj(t) - self.I_Na_v(V, m, h, j) - self.I_K(V, n) - self.I_L(V)) / self.C_m
		#dVdt = -1.0000*(self.I_Na_v(V, m, h, j) + self.I_Ca_l(V, d, f, Ca_i) + self.I_K(V, x) + self.I_Ki (V)+ self.IK_p(V)+self.I_pCa(Ca_i)+self.I_Na_ohm(V)+self.I_Ca_b(V)+self.I_NaK(V, Na_o, Na_i, K_o)) / self.C_m
		dVdt = 1.000*(self.I_inj(t) - self.I_i(V, d, f, Ca_i, x, m, h, j) - self.I_Ca_l(V, d, f, Ca_i))/self.C_m
		dmdt = self.alpha_m(V)*(1.0-m) - self.beta_m(V)*m # Na_v
		dhdt = self.alpha_h(V)*(1.0-h) - self.beta_h(V)*h # Na_v
		djdt = self.alpha_j(V)*(1.0-j) - self.beta_j(V)*j # Na_v
		dndt = self.alpha_n(V)*(1.0-n) - self.beta_n(V)*n
		dddt =  self.alpha_d(V)*(1.00000-d)-self.beta_d(V)*d
		dfdt = self.alpha_f(V)*(1.00000-d)-self.beta_f(V)*f
		dCaidt = -1.0000*(self.I_Ca_b(V)+self.I_Ca_lcanak(V, d, f, Ca_i)+self.I_pCa(Ca_i)+self.I_nsCa(V, Ca_i)+self.I_Ca_lcanak(V, d, f, Ca_i)+self.I_CaK(V, d, f, Ca_i)+self.I_CaNa(V, d, f, Ca_i))#*self.AGeo_cm(r,l)/(self.V_cell(r,l)*self.z_iCa*self.F)		
		#dBdt = -1.00000*(self.I_B*self.A_)
		dxdt =  self.alpha_x(V)*(1.00000-x)-self.beta_x(V)*x
		dNaidt = -1.0000*(self.I_Na_ohm(V)+self.I_Na_ghk(V)+self.I_Na_v(V, m, h, j)) # test by SZ 20191227
		return dVdt, dmdt, dhdt, djdt, dndt, dddt, dfdt, dCaidt, dxdt, dNaidt


	def Main(self):
		"""
		Main demo for the Channels Transporters model
		"""
    #m = 0.00171338077730188 m CellML
    #h = 0.982660523699656 h CellML
    #j = 0.989108212766685  j Cel ML
    #d = 0.00302126301779861 slow_inward_current_d_gate (dimensionless)  = CellML  states[5] = 0.00302126301779861
    #f = 0.999967936476325 slow_inward_current_f_gate (dimensionless) (Ref1ofEP17CELLML)"
    #Cai = 0.00017948816388306 intracellular_calcium_concentration (millimolar)

		X = odeint(self.dALLdt, [-84.3801107371, 0.05, 0.6, 0.98, 0.32, 0.00302126301779861, 0.999967936476325, 0.00017948816388306, 0.0417603108167287, 0.00017948816388306], self.t, args=(self,)) #    Ref1ofEP17 or CellML states[0] = -84.3801107371
		V = X[:,0]
		m = X[:,1]
		h = X[:,2]
		j = X[:,3]
		n = X[:,4]
		d = X[:,5]
		f = X[:,6]
		Ca_i = X[:,7]
		x = X[:,8]
		Na_ii = X[:,9]
		inav = self.I_Na_v(V, m, h, j)
		inaohm = self.I_Na_ohm(V)
		inaghk = self.I_Na_ghk(V)
		ical = self.I_Ca_l(V, d, f, Ca_i) # isi changed to ical by SZ
		ik = self.I_K (V, x)# self.I_K(V, n)
		iv = self.I_v(V, Ca_i)
		il = self.I_L(V)
		iki = self.I_Ki(V)
		ikp = self.I_Kp(V)
		ipCa = self.I_pCa(Ca_i)
		iCab = self.I_Ca_b(V)
		iCalcanak = self.I_Ca_lcanak(V, d, f, Ca_i) # iCa changed to ICalcanak by SZ
		#inaca = self.I_NaCa(V)
		inak = self.I_NaK(V)
		finft = self.f_inft(V)
		icaltotal = self.I_Cal_total(V, d, f, Ca_i)
		voi = array(np.arange(0.0, 600.0, 2.0))	
		print (inaohm) 
		#print self.K_i
		print (Na_ii)

		fig, axs = plt.subplots(5, 1, figsize=(9, 9), sharex=True)

		fig.suptitle('Major Ioninc Currents for Cardiac Ventricular Action Potential')
		plt.plot(self.t, Ca_i, 'b', label='$Ca_{t}$')
		plt.plot(self.t, Na_ii, 'r', label='$Na_{t}$')
		#plt.subplot(7,1,1)
		#plt.title('Major Ioninc Currents for Cardiac Ventricular Action Potential')
		axs[0].plot(self.t, V, 'k', label = 'Action Potential')
		#print V
		axs[0].set_xlim(0, 300)
		axs[0].set_ylabel('V (mV)')
		axs[0].set_ylim(-100.0,80.0)
		#plt.grid()
		#plt.legend()

		#plt.subplot(7,1,3)
		axs[1].plot(self.t, inaohm, 'c', label='$I_{Na}$')
		axs[1].set_xlim(0, 300)
		axs[1].set_ylabel('$I_{Na}$ \n ($\\mu{A}/cm^2$)')
		axs[1].set_ylim(-10.0,10.0)
		#plt.grid()
		plt.legend()
		plt.show()
if __name__ == '__main__':
    runner = ChannelsTransporters()
    runner.Main()


