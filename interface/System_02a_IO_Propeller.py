# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 14:18:35 2022

@author: hbusson
"""





from cosapp.base import System, Port
import numpy as np
from System_02b_F_Propeller import *
from scipy.interpolate import CubicSpline



#File a: determine the input and output necessary for the submodel to work
     #input Vs, NE, time_t
     #output p_inl, T_inl, dma, dmf, T_c, xr


############################ INPUTS ####################################################
#Input port to determine what are the necessary input produced previously from other submodels

class PortInputProp(Port):
    """Handles the cosapp input Port for the propeller's model

    The attributes are initiated with cosapp’s setup method.

    Attributes
    ----------

    Vs: float
        This is the current speed of the ship model
    NE: float
        This is the current rotaiton speed of the engine
    R_w : float
        Wind resistance
    R_awl: float
        Wave resistance
    """
   
    def setup(self):
        """Cosapp setup of the Port

        Populates the input Port of the propeller model.
        """
        self.add_variable('Vs')
        self.add_variable('NE')
        self.add_variable('time_t')
      


############################ OUTPUTS ####################################################
#Output port to determine what are the output produced by this model          
class PortOutputProp(Port):
    """Handles the cosapp output Port for the propeller's model

    The attributes are initiated with cosapp’s setup method.

    Attributes
    ----------
    dVs_dt: float
        Ship acceleration
    Tm: float
        Torque of the motor
    Power_Motor: float
        This is the power coming from the motor
    """
    def setup(self):
        """Cosapp setup of the Port

        Populates the output Port of the propeller model.
        """
        self.add_variable('dVs_dt')
        self.add_variable('Tm')
        self.add_variable('Power_Motor')
        

     
        
        
    
        
#Propeller class model, include both input and output port, additional parameters variables 
class Propeller_sys(System):
    """The propeller model implementation as a cosapp System.

    The attributes are initiated with cosapp's setup method.

    Attributes
    ----------

    var_in: PortInputProp
        The input Port of the system
    var_out: PortOutputProp
        The output Port of the system
    z_p: float
        Propeller's blade number
    A_e_A_o: float
        Disk area coefficient
    p_D_p: float
        The pitch to diameter ratio
    D_p: float
        The propeller diameter,
    omega: float
        The ship wake fraction
    rho_sw: float
        Seawater density
    sh_displ: float
        Ship displacement
    m_hydro: float
        An added virtual mass, which is used in order to take into account the hydrodynamic force arising due to the acceleration of a body in a fluid
    c_1: float
        This is a value
 
    cq: float
        Constant for the wageningen propeller series
    ct: float
        Constant for the wageningen propeller series
    s1: float
        Constant for the wageningen propeller series
    t1: float
        Constant for the wageningen propeller series
    u1: float
        Constant for the wageningen propeller series
    v1: float
        Constant for the wageningen propeller series
    s2: float
        Constant for the wageningen propeller series
    t2: float
        Constant for the wageningen propeller series
    u2: float
        Constant for the wageningen propeller series
    v2: float
        Constant for the wageningen propeller series
    n1: float
        Constant for the wageningen propeller series
    n2: float
        Constant for the wageningen propeller series
    thrust: float
        Constant for the wageningen propeller series

    """

    def setup(self):
        self.add_input(PortInputProp, 'var_in')    # define a new input port
        self.add_output(PortOutputProp, 'var_out')  # define a new output port
        # Solitary variables can also be added,
        # as either `inward` or `outward` variables:
             #the number of propeller blades
        self.add_inward('z_p')
        #disk area coefficient
        self.add_inward('A_e_A_o')
        #the pitch to diameter ratio
        self.add_inward('p_D_p')

        # the propeller diameter,
        self.add_inward('D_p')
        # the ship wake fraction
        self.add_inward('omega')
        #The mass, ms, is the mass of the ship, which is calculated by multiplying the ship displacement sh_displ with the sea water density rhp_sw
        self.add_inward('rho_sw')
        self.add_inward('sh_displ')
        # is an added virtual mass, which is used in order to take into account the hydrodynamic force arising due to the acceleration of a body in a fluid
        self.add_inward('m_hydro')
        #
        self.add_inward('c_1')
        
    
        
        #constants for wageningen propeller series
        self.add_inward('cq')
        self.add_inward('ct')
        self.add_inward('s1')
        self.add_inward('t1')
        self.add_inward('u1')
        self.add_inward('v1')
        self.add_inward('s2')
        self.add_inward('t2')
        self.add_inward('u2')
        self.add_inward('v2')
        self.add_inward('n1')
        self.add_inward('n2')
        self.add_inward('thrust')
        
        self.add_inward('kt',  dtype= np.ndarray,value=np.array([]))
        self.add_inward('kq', dtype= np.ndarray,value=np.array([]))
     
        
        

############################ Computation ####################################################     
#####################What this submodel do is explicited in the System_02b_F_Propeller file################################

    def compute(self): # `compute` method defines what the system does

          """Cosapp setup of the model

          This function populates cosapp’s system of the compressor’s model. It uses PortInputProp as an input port and
           PortOutputProp as an output port. The other variables are cosapp’s inward.
          """
          Prop=Propeller(self.var_in.Vs, self.var_in.NE, self.z_p,self.A_e_A_o, self.p_D_p,  self.D_p, self.omega,  
                           self.rho_sw, self.sh_displ, self.m_hydro, self.c_1, self.thrust, 
                           self.cq, self.ct,self.s2,self.t2,self.u2,self.v2,
                           self.s1,self.t1,self.u1,self.v1,
                           self.n1,self.n2, self.var_in.time_t)
          
          
          
          Prop.prop_comput()
          Prop.Velocity_derivative()
          
          Prop.Motor()
          
          Prop.checkresult()
          
          
          
          #setting the outputs: boat speed, motor power
          
          self.var_out.Tm = Prop.Tm
          self.var_out.dVs_dt = Prop.dVs_dt
          self.var_out.Power_Motor = Prop.Power_Motor
          
   
       
        
        
    
#input = Vs, NE
#output = dvs/dt, Tm   