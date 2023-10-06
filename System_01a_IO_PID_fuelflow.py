# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 08:54:34 2022

@author: hbusson
"""

from cosapp.base import System, Port
import numpy as np
from System_01b_F_PID_fuelflow import *

#File a: determine the input and output necessary for the submodel to work
     #input NE, NE_current, NTC_Current, Nord, p_ER
     #output p_inl, T_inl, dma, dmf, T_c, xr


############################ INPUTS ####################################################
#Input port to determine what are the necessary input produced previously from other submodels

class PortInputP1WS(Port):
    """Handles the input cosapp Port for the fuel flow model.

    The attributes are initiated with cosapp’s setup method.

    Attributes
    ----------

        NE: np.ndarray
            Contains the engine's speed

        NE_Current:
            Current engine's speed
        Nord: np.ndarray
            Contains the speed order
        time_t: float
            time t
        count: int
            iteration number
    """
    def setup(self):
        #vector of  engine speed (past and present value)
        self.add_variable('NE' , dtype= np.ndarray, value= np.zeros(60))
        #current rotation speed of engine and turbo charge
        
        self.add_variable('NE_Current')
        #vector of order of boat speed
        self.add_variable('Nord' , dtype= np.ndarray, value=np.zeros(60))
        #time t
        self.add_variable('time_t')
        #iteration number 
        self.add_variable('count')
     

############################ OUTPUTS ####################################################
#Output port to determine what are the output produced by this model           
class PortInlet(Port):
    """Handles the output cosapp Port for the fuel flow model.

    The attributes are initiated with cosapp’s setup method.

    Attributes
    ----------

    dmf: float
        Fuel mass flow
    t_inj: float
        timing of injeciton
    mf_cy: float
        Quantity of fuel per cycle
    """
    
    def setup(self):
    
        #fuel mass flow
   
        self.add_variable('dmf')
       
        #rack position
        self.add_variable('t_inj')
        
          
        #quantity of fuel per cycle
        self.add_variable('mf_cy')
  
#Compressor and air cooler class model, include both input and output port, additional parameters variables            

class PID_sys(System):
    """The engine model implementation as a cosapp System.

    The attributes are initiated with cosapp’s setup method.

    Attributes
    ----------
    var_in: PortInputP1WS
        The input Port of the system
    var_out: PortInlet
        The output Port of the system
    rev_cy : float
        The number of revolutions per cycle
    zc : int
        The engine's number of cylinder

    kp : float
        Linear coefficient for PID regulation
    kd : float
        Derivative's coefficient for PID regulation
    ki : float
        Integral's coefficient for PID regulation
    xr_o :float
        Original time of injection
    """

    def setup(self):
        """Cosapp setup of the model

        This function populates cosapp's system of the PID's model. It uses PortInputP1WS as an input port and
        PortInlet as an output port. The other variables are cosapp's inward.
        """
        self.add_input(PortInputP1WS, 'var_in')    # define a new input port
        self.add_output(PortInlet, 'var_out')  # define a new output port
        # Solitary variables added,
        # parameters that are not produced by other submodels:
            #number of revolution per cycle
            
        self.add_inward('rev_cy')
            #number of engine cylinders
        self.add_inward('zc')
 
        

   
            #constants for PID controller
        self.add_inward('kp')
        self.add_inward('kd')
        self.add_inward('ki')
      
            #original time of injection (at the beginning of the simulation)
        self.add_inward('tinj_o')
        
        
        

############################ Computation ####################################################     
#####################What this submodel do is explicited in the System_01b_F_Compressor file################################        
      
     #compute function determine what this sub model do
    def compute(self): # `compute` method defines what the system does
       
       """The compute method to run the simulation.

       It sets up the Controller_fuel class to run the simulation (we could probably optimize it with a bit of OOP)
       """
        
     
       PID_ff =Controller_fuel(self.var_in.NE_Current, self.kp,self.ki,self.kd, self.var_in.NE,self.var_in.Nord ,
                          self.tinj_o,self.var_in.time_t,  self.rev_cy, self.zc, self.var_in.count)
        

        
       PID_ff.running()
        
        #setting the outputs: inlet pressure and temperature, mass flows, rack position and cooler temperature

        
      
       self.var_out.dmf = PID_ff.mdot_f
        
       self.var_out.t_inj = PID_ff.t_inj
        
       self.var_out.mf_cy = PID_ff.mf_cy
       
       PID_ff.checkresult()
        
            
            
            

       
        
        
        