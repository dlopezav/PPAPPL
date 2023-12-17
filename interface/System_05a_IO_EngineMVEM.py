# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 15:01:32 2022

@author: hbusson
"""


from cosapp.base import System, Port
import numpy as np
from System_05b_F_EngineMVEM import *

#File a: determine the input and output necessary for the submodel to work
     #input Q_p, Qe, mdot_a,mdot_e, T_exh, Ne, x_r, NTC, T_c
     #output dNe_dt, dNC_dt

############################ INPUTS ####################################################
#Input port to determine what are the necessary input produced previously from other submodels

class PortInputMVeng(Port):
    """Handles the input cosapp Port for the engine's MVEM model.

    The attributes are initiated with cosapp's setup method.

    Attributes
    ----------
    Q_p: float
        Torque required by the propeller 
 
    p_IR: float
        Pressure at the inlet receiver
    T_IR: float
        Temperature at the inlet receiver
    p_exh: float
        Pressure of the exhaust gas
    T_exh: float
        Temperature of the exhaust gaz

    NTC: float
        Turbo charger speed
    mf_dot: float
        Fuel mass flow
 
    W_dot: float
        The derivative of the work
    Qw_dot: float
        The derivation of the heat transfer
    NE_Current: float
        Engine speed
    deltat: float
        The timestep for the simulation
    pstor : vector
        storage of pressure
    Tstor : vector
        storage of temperature
    mestor : vector
        storage of 
    """
    def setup(self):
        """Cosapp setup of the Port

        Populates the input Port of the engine's MVEM model.
        """
        #engine and propeller torque
        self.add_variable('Q_p')
    
     
 
        self.add_variable('p_IR')
        self.add_variable('T_IR')
        self.add_variable('p_exh')
        self.add_variable('T_exh')
        #current speed of engine and turbo charger
        self.add_variable('NTC')
        #fuel flow
        self.add_variable('mf_dot')

        self.add_variable('W_dot')
        self.add_variable('Qw_dot')
        self.add_variable('NE_Current')
        self.add_variable('deltat')
        
        
        self.add_variable('pstor')
        self.add_variable('Tstor')
        self.add_variable('mestor')
        
        self.add_variable('prt')
        self.add_variable('prc')
        

############################ OUTPUTS ####################################################
#Output port to determine what are the output produced by this model        
class PortDerivOut(Port):
    """Handles the input cosapp Port for the engine's MVEM model.

    The attributes are initiated with cosapp's setup method.

    Attributes
    ----------
    dNE_dt: float
        Engine rotation speed acceleration
    dNTC_dt: float
        Turbo charger rotation speed acceleration
    T_IR: float
        Temperature at the inlet receiver
    p_IR: float
        Pressure at the inlet receiver
    p_exh: float
        Pressure of the exhaust gas
    T_exh: float
        Temperature of the exhaust gas
    pstor : vector
        storage of pressure
    Tstor : vector
        storage of temperature
    mestor : vector
        storage of 
    prt : vector
        pressure ratio turbine
    prc : vector
        pressure ratio compressor
    """
    def setup(self):
        """Cosapp setup of the Port

        Populates the output Port of the engine's MVEM model.
        """
        #two output variables
        #engine speed derivative
        self.add_variable('dNE_dt')
        #turbo charger speed derivative
        self.add_variable('dNTC_dt')
        
        self.add_variable('T_IR')
        self.add_variable('p_IR')
        self.add_variable('p_exh')
        self.add_variable('T_exh')
        
        self.add_variable('pstor')
        self.add_variable('Tstor')
        self.add_variable('mestor')
        
        self.add_variable('prt')
        self.add_variable('prc')


  
#MVEM class model, include both input and output port, additional parameters variables        

class Rotationspeed(System):
    """The engine model implementation as a cosapp System.

    The attributes are initiated with cosapp's setup method.

    Attributes
    ----------

    var_in: PortInputMVeng
        The input Port of the system.
    var_out: PortDerivOut
        The output Port of the system.

    p_a: float
        Air pressure
    T_a: float
        Air temperature
    R_a: float
        R divided by air's molar mass
    T_w: float
        Water temperature
    k_ep: float
        Pressure turbine
    k_ac0: float
        Air cooler efficiency's polynome's coefficient
    k_ac1: float
        Air cooler efficiency's polynome's coefficient
    k_ac2: float
        Air cooler efficiency's polynome's coefficient

    kf0: float
        friction pressure as a function of NE polynome's coefficient 
    kf1: float
       friction pressure as a function of NE polynome's coefficient 
    kf2: float
        friction pressure as a function of NE polynome's coefficient 

    kpac: float
        This is a value
    kpaf: float
        This is a value
    # ratio of specific heat of ambiant and exhaust
    gamma_a: float
        Ratio for the specific heat of ambiant air
    gamma_e: float
        Ratio for the specific heat of exhaust gas
    eta_t: float
        Tubine efficiency
    eta_sh: float
        Shaft efficiency
    eta_c: float
        Compressor efficiency
    eta_comb: float
        Combustion efficiency
    cpa: float
        Air's heat capacity
    cpe: float
        Exhaust gas' heat capacity
    cpc: float
        heat capacity of gas in compressor
    cpt: float
       heat capazcity of gas in turbine
    rev_cy: float
        Number of revolution per cycle
    # different inertia
    ITC: float
        Turbo charger inertia
    I_e: float
        Engine inertia
    I_sh: float
        Shaft inertia
    I_p: float
        Propeller inertia
    rc: float
        Engine's compression ratio
    V_D: float
        Engine's displacement volume
    cd: float
        Discharge coefficient
    A_eq: float
        Equivalent cylinders flow area (Aeq)
    lhv: float
        Lower heating value
    cvir: float
        Inlet receiver Cv
    cver: float
        Exhaust receiver Cv
    V_ir: float
        Inlet receiver volume
    V_er: float
        Exhaust receiver volume
    mdot_c: float
        Compressor's mass flow
    mdot_t: float
        Turbine's mass flow
    """

    def setup(self):
        """Cosapp setup of the model

        Populates cosapp's system of engine's MVEM model. It uses the PortDerivOut and the PortInputMVeng
        classes as ports, and the other variables as cosapp's inwards.
        """
        self.add_input(PortInputMVeng, 'var_in')    # define a new input port
        self.add_output(PortDerivOut, 'var_out')  # define a new output port
        # Solitary variables added,
        # parameters that are not produced by other submodels:
            #air pressure and temperature
        self.add_inward('p_a')
        self.add_inward('T_a')
        self.add_inward('R_a')
        #temperature water
        self.add_inward('T_w')
        
        #exhaust system geometry characteristics
        self.add_inward('k_ep')
        #air cooler effectiveness coefficient epsilon = kAC0+  kAC1*mdot_a+ kAC2*mdot_a^2
        self.add_inward('k_ac0')
        self.add_inward('k_ac1')
        self.add_inward('k_ac2')
        
        self.add_inward('kf0')
        self.add_inward('kf1')
        self.add_inward('kf2')
        
        self.add_inward('kpac')
        
        self.add_inward('kpaf')
     
        
            
        #ratio of specific heat of ambiant and exhaust
        self.add_inward('gamma_a')
        self.add_inward('gamma_e')
        
    
        
          #turbine efficiency NOTE it is a function of turbine pressure ratio but for the moment taken as a constant  input
        self.add_inward('eta_t')
        #shaft efficiency
        self.add_inward('eta_sh')
        
        #compressor efficiency constant
        self.add_inward('eta_c')
        #combustion efficiency
        self.add_inward('eta_comb')
        
        
           #heat capacity of air and exhaust
        self.add_inward('cpa')
        self.add_inward('cpe')
        self.add_inward('cpc')
        self.add_inward('cpt')
           #number of revolutions per cycle
        self.add_inward('rev_cy')
        
   
            #different inertia
        self.add_inward('ITC')
        self.add_inward('I_e')
        self.add_inward('I_sh')
        self.add_inward('I_p')
      
        

   
     
 
        
    
        
     
        # engine compression ratio
        self.add_inward('rc')
  
        #the engine displacement volume
        self.add_inward('V_D')
        
        #discharge coefficient
        self.add_inward('cd')
        #equivalent cylinders flow area (Aeq)
        self.add_inward('A_eq')

        
        #lower heating value
        self.add_inward('lhv')
        
      

        #inlet and exhaust receivers characteristics: cv and volume
        self.add_inward('cvir')
        self.add_inward('cver')
        
        self.add_inward('V_ir')
        self.add_inward('V_er')
        
        
        #to read in tables
        self.add_inward('mdot_c')
        self.add_inward('mdot_t')



 
          
      
    
        
    
        
############################ Computation ####################################################     
#####################What this submodel do is explicited in the System_05b_F_EngineMVEM file################################          
        

    def compute(self): # `compute` method defines what the system does
       
            """The compute method to run the simulation.

            It sets up the Enginederivatives class to run the simulation (we could probably optimize it with a bit of OOP)

            """
            #What this submodel do is explicited in the System_05b_F_EngineMVEM file  
            ED=Enginederivatives(self.p_a, self.k_ep, self.gamma_e,self.eta_sh, self.eta_t, self.eta_c, self.I_e, self.I_sh, 
                                 self.I_p,
                                 self.rev_cy, 
                                 self.V_D, self.ITC, self.cpa, self.cpe, self.var_in.Q_p, 
                                 self.var_in.mf_dot, self.mdot_c, self.mdot_t,
                                 self.T_a,
                                 self.A_eq, self.gamma_a, self.var_in.p_IR, self.var_in.T_IR, self.var_in.W_dot,
                                 self.var_in.Qw_dot, self.eta_comb, self.lhv,
                                self.var_in.p_exh,self.kpac, self.kpaf, self.cpc, self.cpt, self.cver, self.cvir, self.V_ir, self.V_er, self.cd, 
                                self.R_a,
                                 self.var_in.NE_Current,self.kf0, self.kf1, self.kf2, self.k_ac0, self.k_ac1, self.k_ac2, self.T_w,
                                self.var_in.T_exh, self.var_in.deltat, self.var_in.NTC, self.var_in.prc, self.var_in.prt)
            

            
            ED.Mass_flow()
            
            ED.Combining_cycles()
            print('mapping')
            ED.Mapping()
            
            ED.pressure_ratios()
            
            
            ED.Temp_and_pressure()
    
            ED.Rotation_derivatives(self.var_in.NTC)
          
            ED.Inlet_receiver()
            
            ED.Exhaust_receiver()
            
            ED.checkresult()

         
            #setting the outputs: both speed derivatives
            self.var_out.dNE_dt = ED.dNE_dt
            self.var_out.dNTC_dt = ED.dNTC_dt
            
            self.var_out.p_IR = ED.p_IR
            self.var_out.T_IR = ED.T_IR
            self.var_out.p_exh = ED.p_exh
            self.var_out.T_exh = ED.T_exh
            
            self.var_out.prc = ED.prc
            self.var_out.prt = ED.prt
            
            self.var_out.pstor =  np.append(self.var_in.pstor,ED.p_IR)
            self.var_out.Tstor = np.append(self.var_in.Tstor,ED.T_IR)
            self.var_out.mestor = np.append(self.var_in.mestor,ED.mir)
            
        
        
        
   