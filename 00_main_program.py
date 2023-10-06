# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 15:37:58 2022

@author: hbusson
"""
#initial mass in cylinder : ideal gas law
#initial mass receivers : ideal gas law as we have tir pir for both


#libraries needed


import numpy as np
from scipy.interpolate import CubicSpline
from System_01a_IO_PID_fuelflow import *
from System_02a_IO_Propeller import *
from System_03a_IO_Powermanagement import *
from System_04a_IO_Engine0d import *
from System_05a_IO_EngineMVEM import *
from Init_system import *
from math import *



class Current:
    
    def __init__(self, NE_Current, NTC_Current, Vs_Current):


        #initial conditions
        self.NE_Current = NE_Current
        self.NTC_Current = NTC_Current
        self.Vs_Current = Vs_Current
        
    def update(self,NE, NTC, Vs):
        
        self.NE_Current = NE
        self.NTC_Current = NTC
        self.Vs_Current = Vs
        




PID=PID_sys(name='PID')
baldi0d=InCylinder(name='baldi0d')
Prop = Propeller_sys('Prop')
EngineRot = Rotationspeed('EngineRot')
PowerMag = Powermanagement('PowerMag')

############################### MODEL PARAMETERS #######################################################
############################### MODEL PARAMETERS #######################################################
############################### MODEL PARAMETERS #######################################################
############################### MODEL PARAMETERS #######################################################
############################### MODEL PARAMETERS #######################################################


time_t=0
total_time=0.005 

deltat2= 0.005


 
N_cycle = int(total_time // deltat2+1)


################ PARAM PID Propeller BEGIN #######################################




#ship displacement volume [m3]
sh_displ=112404



#density of sea water [kg/m^3]
rho_sw = 1026


################## to be determined ###################################

#mass of ship [kg]

mass_of_ship =260000000

#add virtual mass [kg]
m_hydro = rho_sw*sh_displ

#diameter [m]
D_p =10

#thrust
thrust = 0.2

# dVs_dt = c_1 * V_s^2 [kg/m]
c_1=26250

# the number of propeller blades, zp
z_p = 6
# the disk area coefficient, AE/Ao
A_e_A_o = 0.82
#the pitch to diameter ratio, p/Dp
p_D_p = 0.93846154


#ship wake fraction, w, which is considered constant taking values in the range from 0.20 
#to 0.45 for ships with a single propeller
omega = 0.3


########################wageningen coefficients imports ########################################


import pandas as pd

df = pd.read_excel (r'C:/Users/hbusson/Documents/TNTM/PythonScripts/specific ship code/champs_elysees/wageningen.xlsx')


n1_nump = df['n1'].to_numpy()
ct_nump = df['ct'].to_numpy()
s1_nump = df['s'].to_numpy()
t1_nump = df['t'].to_numpy()
u1_nump = df['u'].to_numpy()
v1_nump = df['v'].to_numpy()


n2_nump = df['n2'].to_numpy()
cq_nump = df['cq'].to_numpy()
s2_nump = df['s2'].to_numpy()
t2_nump = df['t2'].to_numpy()
u2_nump = df['u2'].to_numpy()
v2_nump = df['v2'].to_numpy()





################ PARAM Propeller END #######################################

################ PARAM PID Controller BEGIN #######################################


zc=12	#number of engine cylinders	


############### to be determined ###############################################
kp	= 0.003#constants for PID controller	3,00E-07
kd  = 0	 #constants for PID controller	4,00E-07
ki  =	 0.0015 #constants for PID controller	0
xr_o=0.01427586	 #original injeciton time (at the beginning of the simulation)	0,015






#


steadystate= pd.read_excel(r'C:/Users/hbusson/Documents/TNTM/ship data/champs_elysees/steady_state1.xlsx')

#order of speed 

Nord = 2*math.pi/60*steadystate['me_ne'].to_numpy()















################ PARAM PID Controller END #######################################


################ PARAM 0d cycle BEGIN #######################################


deltat=0.0001

#number of revolution per cycle: 
rev_cy=1


Vc=	164671 #clearance volume [cm3]
#lower heating value 
lhv=	42707

#wall surface [units ?]
Sw	=1.102343986
#stroke
stroke	=3.468
l=	3468
ar=	1734
#perfect gas constant
R=	8.314462


##### to be confirmed ############################

Tw=	300+273.15  # [K]
p_ivc=	4.1 #pressure at inlet valve closing [bar]
T_ivc=	307 #temperature at the inlet valve closing [K]
ma_ini=	7462 #initial masse in the cylinder [g]
m_ivc = ma_ini 

#compression ratio
rc	=15

#reference temperature
Tref=	298.15


#stochiometric coefficient
Cstoich=	14.8

#perfect gas constant divided by molar mass of exhaust, air and fuel gas respectively
re=	 0.1341 # J K-1 g-1
ra=	 0.287 # J K-1 g-1
rf=	 0.0489086 # J K-1 g-1




#Wiebe parameters parameters estimated previously 
acomb =  0.42
bcomb=  0.49
am=0.59
bm=-0.21
cm = -0.96
delta_m = 0.27


#parameters taken in simulating combustion p167
a_id=0.39
b_id=0.105
c_id=3.12
mWiebe=0.1







#Engine reference point that was used to calibrate the Wiebe curve
#air fuel ratio for reference point
lamb_ref = 2.435696566131349
#ignition delay for reference point
phi_id_ref = 1
#mass at inlet valve closing for reference point [g]
m_ivc_ref = 7462
# engine speed for reference point [rpm]
N_ref = 80
# wiebe curve form coefficient for reference point
mWieberef =   0.1
#combustion duration for reference point 
delta_phi_comb_ref= 0.3667076110839844*60/(2*math.pi)

SMFR= 14500
topen = 0.0005
tclose= 0.0005



################ PARAM 0d cycle END #######################################


################ PARAM MVEM Model BEGIN #######################################

#pressure air [bar]
p_a =1
# air Temperature [K]
T_a = 11+273.15

#heat specific capacity ratio for air
gamma_a= 1.401114206
#heat specific capacity ratio for exhaust
gamma_e = 1.375


#coefficient for pressure after turbine pt = pa + k_ep* me^2 [bar g^-2 s^2]
k_ep=4.20*pow(10,-12)



# shaft efficiency
eta_sh =    0.99
#combustion efficiency
eta_comb=0.99
#turbine efficiency
eta_t=0.84

#compressor efficiency
eta_c=0.8


#initial pressures and temperatures in inlet receiver IR and exhaust rceiver ER
#pressure exhaust [bar]
p_ER_ini =3.8
#temperature exhuast [K]
T_ER_ini =643.15
#temperature inlet [K]
T_IR_ini =298.15
#pressure inlet [bar]
p_IR_ini =2.4

# heat capacities, air and exhaust 
cpe = 1.100
cpa =  1.006

# coeffcient for pressure drop in air cooler delta_pac = kpac*pow(mdot_c,2)  [bar*s^2/g^2]
kpac=3.80*pow(10,-12)
# coeffcient for pressure drop in air filter delta_pac = kpac*pow(mdot_c,2)  [bar*s^2/g^2]
kpaf=0

# cv for inlet and exhaust
cvir=0.718
cver=0.8
#volume of inlet receiver [l or dm^3]
V_ir=20.000
# volume for exhaust receiver [l or dm^3]
V_er=20.000




#friction mean effective pressure is considered function of  the engine crankshaft speed, these are the coefficients
 # pressure friction coefficients as a function of engine speed p_f = kf0+ kf1*NE + kf2*NE^2
kf0 = 2.5645
kf1 = -0.0532
kf2=  0.0005



#air cooler effectiveness is assumed to be a function of the air mass flow rate, these are the coefficients
#  epsilon = kAC0+  kAC1*mdot_a+ kAC2*mdot_a^2

k_ac0 = 0.95
k_ac1 = -5.00*pow(10,-6)
k_ac2 = 7*pow(10,-11)




# diametre cylindre [m]
bore=0.92
#number of cylinders 
zc=12
# volume displacement 
V_D = zc*math.pi* pow(bore,2)*stroke/4


# Different Intertia: turbocharger, engine, shaft, propeller
 #turbocharger inertia
ITC = 600
#engine intertia
I_e =0   
#shaft inertia
I_sh = 0
#propeller inertia
I_p = 500000

#mass calculations
#coefficient discharge
cd = 1
# Area
A_eq=1



R_a = 287

# the temperature of the air cooler coolant medium [K]
T_w = 32.32+273.15



################ PARAM MVEM Model END #######################################


################ PARAM Optimisation energy demand BEGIN #######################################


Max_power = 30
Aux_Power_demand = np.random.randint(int(Max_power), size=(N_cycle))

#indexes : ME1 = Main Engine 1,ME2 = Main Engine 2, AE = Auxiliary Engine1, AE2 = Auxiliary Engine1

#maximum power for each engine [W]
PME1_max=63840000
PME2_max=1
PAE1_max=1
PAE2_max=1

 #minimum power for each engine [W]
PME1_min=0.0
PME2_min=0.0
PAE1_min=0.0
PAE2_min=0.0

 #lower heating values
LHV1=100.0
LHV2=0.0
LHV3=0.0
LHV4=0.0
LHV5=0.0
 
#maximum continuous rating
MCR1=63840000
MCR2=0
MCR3=100
MCR4=0
MCR5=0
 

#coefficients for fonction for fuel cons
an=1.0
an2=1.0
an3=1.0
an4=1.0
an5=1.0
 
bn=1.0
bn2=1.0
bn3=1.0
bn4=1.0
bn5=1.0
 
cn=1.0
cn2=1.0
cn3=1.0
cn4=1.0
cn5=1.0
 
#efficiencies
nu1ME1=0.8
nu1ME2=0.8
nu1AE1=0.8
nu1AE2=0.8

 
nu2ME1=0.8
nu2ME2=0.8
nu2AE1=0.8
nu2AE2=0.8


 
nuelME1=0.8
nuelME2=0.8
nuelAE1=0.8
nuelAE2=0.8

 
Q_ab_max=30.0
Q_ab_min=3

################ PARAM Optimisation energy demand END #######################################

############################### MODEL PARAMETERS #######################################################
############################### MODEL PARAMETERS #######################################################
############################### MODEL PARAMETERS #######################################################
############################### MODEL PARAMETERS #######################################################
############################### MODEL PARAMETERS #######################################################
    
############################### Initialisations #######################################################
############################### Initialisations #######################################################

#initial engine speed 80 rpm translated into   rad per seconds       
NE_ini =49*2*math.pi/60

#initial value of turbocharger speed [unit?]
NTC_ini = 7522
#initial value of engine rotation speed
NE_Current = NE_ini
#initial pressure and temperature in exhaust receiver
T_exh = T_ER_ini
p_exh = p_ER_ini
#initial pressure and temperature in exhaust receiver
T_IR = T_IR_ini
p_IR = p_IR_ini

#initial pressure ratio compressor
prc = 2.4
#initial pressure ratio turbine
prt = 3.8

#Initial_speed

Vs_ini = 12.15


#model time of the simulation
total_time = 0.035



############################### Initialisations #######################################################
############################### Initialisations #######################################################






##############################Initialisation#######################################

C_Current= Current(NE_ini, NTC_ini, Vs_ini)


NE_tot = np.array([NE_ini])
NTC_tot = np.array([NTC_ini])

#Putting the different parameters in every single sub models
Initialisations(  
                  #PID controller model #
                  PID, rev_cy, kp, kd, ki, 
                  xr_o, time_t, T_a, p_a, k_ac0,
                  k_ac1, k_ac2, T_w,  V_D,  cd,
                  A_eq, R_a, gamma_a, eta_c,
                  
                  #Engine 0d Baldi #
                  baldi0d,deltat, zc,
                  delta_phi_comb_ref, lamb_ref, N_ref, m_ivc_ref,phi_id_ref,
                   Vc,
                  mWieberef, lhv,  Tw, p_ivc, T_ivc, ma_ini, mWiebe,
                  a_id,b_id,c_id,stroke,rc,acomb, bcomb, am,bm,cm, 
                  Tref, delta_m, Sw, ra,re,rf, l,ar, Cstoich,R, kf0, kf1, kf2,bore,
                  SMFR, topen, tclose,
                  #Propeller model #
                  Prop,  z_p,A_e_A_o, p_D_p, D_p, omega,  
                  rho_sw, sh_displ, m_hydro, c_1, I_p,pi,cq_nump,ct_nump,s2_nump,t2_nump,u2_nump,v2_nump,s1_nump, 
                  t1_nump,u1_nump,v1_nump,n1_nump,n2_nump,thrust,
                  #Optimisation energy model #
                  PowerMag,  PME1_max,PME2_max, PAE1_max,  PAE2_max,
                  PME1_min,PME2_min, PAE1_min,  PAE2_min,
                  LHV1, LHV2, LHV3, LHV4, LHV5, 
                  MCR1, MCR2, MCR3, MCR4, MCR5, 
                  an,  an2,  an3,  an4,  an5,
                  bn,bn2,bn3,bn4,bn5,
                  cn, cn2, cn3, cn4, cn5, 
                  nu1ME1,  nu1ME2, nu1AE1, nu1AE2, 
                  nu2ME1,nu2ME2,nu2AE1,nu2AE2,
                  nuelME1,   nuelME2,   nuelAE1,   nuelAE2, 
                  Q_ab_max, Q_ab_min,
                  #MVEM model #
                  EngineRot,  k_ep, gamma_e, eta_sh, eta_t, I_e, I_sh,
                   ITC, cpa, cpe, eta_comb, kpac, kpaf, cvir, cver, V_ir, V_er )







count = 0

#storage of pressure, temperature, time and mass flow
pstor = np.array(p_IR_ini)
Tstor = np.array(T_IR_ini)
mestor = np.array(90919)
tstor = np.array([0])




while time_t <total_time:
  
########################################### Sub model PID Controller ###########################################
  

#Inputs for PiD fuel flow #
    #current time
    PID.var_in.time_t=time_t
    #Vector of engine speeds (present and past values)
    PID.var_in.NE= NE_tot
    #Current engine speed
    PID.var_in.NE_Current= C_Current.NE_Current

    #vector of desired speed
    PID.var_in.Nord = Nord
    #counting
    PID.var_in.count = count
    
    #running PID fuel flow model found in System_01a_IO_Compressor: "compute" function is called
    PID.run_once()



########################################### Sub model Propeller ###########################################

 

#Inputs for Propeller #

    #Current time of simulation
    Prop.var_in.time_t = time_t
    # Current speed of ship 
    Prop.var_in.Vs = C_Current.Vs_Current
    # Current engine speed
    Prop.var_in.NE=C_Current.NE_Current

#running Propeller model found in System_02a_IO_Propeller: "compute" function is called
    Prop.run_once()


########################################### Sub model Engine 0d cycle ###########################################


#Inputs for Engine 0d cycle #
 
    #Temperature inlet 
    baldi0d.var_in.T_IR=T_IR
    #Pressure inlet 
    baldi0d.var_in.p_IR=p_IR
    #Current engine speed 
    baldi0d.var_in.NE_Current=C_Current.NE_Current
    
    baldi0d.var_in.mf_cy=PID.var_out.mf_cy
    
    baldi0d.var_in.t_inj = PID.var_out.t_inj

#running 0d model found in System_04a_IO_Engine0d: "compute" function is called    
    baldi0d.run_once()
    
    
    



########################################## Sub model Engine MVEM ############################################


#Inputs for Engine MVEM cycle #

    #Power requested requested Torque 
    EngineRot.var_in.Q_p =  Prop.var_out.Power_Motor/C_Current.NE_Current

    EngineRot.var_in.deltat= deltat
 
 
    
    #Temperature inlet receiver
    EngineRot.var_in.T_IR= T_IR
    #pressure inlet receiver
    EngineRot.var_in.p_IR= p_IR
    #pressure exhaust receiver
    EngineRot.var_in.p_exh= p_exh
    #temperature exhaust receiver
    EngineRot.var_in.T_exh= T_exh
    
    
    #pressure ratio for turbine
    EngineRot.var_in.prt= prt
    #pressure ratio for compressor
    EngineRot.var_in.prc= prc
    # fuel flow 
    EngineRot.var_in.mf_dot = PID.var_out.dmf
    # current engine speed
    EngineRot.var_in.NE_Current= C_Current.NE_Current
    # current turbo charger speed
    EngineRot.var_in.NTC=C_Current.NTC_Current
 



      

    EngineRot.var_in.W_dot= baldi0d.ex.W_dot
    EngineRot.var_in.Qw_dot= baldi0d.ex.Qw_dot

#running Engine Rotation found in System_05a_IO_EngineMVEM: "compute" function is called     
    EngineRot.run_once()
      
    Prop.run_once()

########################################## Sub model Next time step ############################################
    
    #taking the different derivatives computed: output of previous sub models
    #derivative of TC speed output of engine MVEM sub model
    dNTC = EngineRot.var_out.dNTC_dt
    #derivative of engine speed output of engine MVEM sub model
    dNE = EngineRot.var_out.dNE_dt

    # updating current engine and turbo charger speed and the ship speed
    NE_Current = C_Current.NE_Current*60/(2*math.pi) + dNE*deltat2
    NTC_Current = C_Current.NTC_Current + dNTC*deltat2
    
    NE_Current = 2*math.pi/60*NE_Current
  
    C_Current.update(NE_Current, NTC_Current, Vs_ini)
    
    p_IR = EngineRot.var_out.p_IR
    p_exh = EngineRot.var_out.p_exh
    T_IR = EngineRot.var_out.T_IR
    T_exh = EngineRot.var_out.T_exh
    
    prc = EngineRot.var_out.prc
    prt = EngineRot.var_out.prt
    
    # time t
    
    time_t = time_t +deltat2
    tstor = np.append(tstor, time_t)
    pstor = np.append(pstor, p_IR)
    Tstor = np.append(Tstor, T_IR)

    
   

    
    # putting actual engine speed at the end of the vector of speeds
    NE_tot=np.append(NE_tot,NE_Current)
    NTC_tot=np.append(NTC_tot,NTC_Current)
    
    
    count = count+1
    
    print(9.5493 *NE_tot)
    
    print(NTC_tot)

    
  





