# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 10:05:48 2022

@author: hbusson
"""

from cosapp.base import System, Port
import numpy as np
from System_04b_F_Engine0d import *
from main import *


import math

class PortInputengine(Port):
    """Handles the input cosapp Port for the engine's 0d model.

    The attributes are initiated with cosapp's setup method.

    Attributes
    ----------

    p_IR
        Pressure at the inlet receiver
    T_IR
        Temperature at the inlet receiver
    NE_current
        Current engine speed
    mf_cy
        Mass of fuel injected per cycle
    t_inj
        injection duration
    """
    def setup(self):
        """Cosapp setup of the Port

        Populates the input Port of the engine's 0d model.
        """
        self.add_variable('p_IR')
        self.add_variable('T_IR')
        self.add_variable('NE_Current')
        self.add_variable('mf_cy')
        self.add_variable('t_inj')
        
class PortExhaust(Port):
    """Handles the output cosapp Port for the engine's 0d model.

    The attributes are initiated with cosapp's setup method.

    Attributes
    ----------

    p_exh
        This is a vector regrouping all the successive pressures during the simulation of the combustion, last one is pressure at exhaust receiver
    Qf
        This is a vector regrouping all the successive pressures during the simulation of the combustion,, successive heat created by combustion
    T_exh
        This is a vector regrouping all the successive temperatures during the simulation of the combustion, last one is temperature at exhaust receiver
    V_out
        This is a vector regrouping all the successive volumes during the simulation of the combustion
    W_dot
        This is a vector regrouping all the successive works (pressure * volume) during the simulation of the combustion
    Qw_dot
        This is a vector regrouping all the successive heat exchanges during the simulation of the combustion
  
    """
    def setup(self):
        """Cosapp setup of the Port

        Populates the output Port of the engine's 0d model.
        """
        self.add_variable('p_exh',  dtype= np.ndarray,value=np.array([]))
  
        self.add_variable('Qf', dtype= np.ndarray,value=np.array([]))
        self.add_variable('T_exh',  dtype= np.ndarray,value=np.array([]))
        self.add_variable('V_out',  dtype= np.ndarray,value=np.array([]))

        self.add_variable('W_dot')
        self.add_variable('Qw_dot')
      
        
     
        

class InCylinder(System):
    """The engine's 0d model implementation as a cosapp System.

    The attributes are initiated with cosapp's setup method.

    Attributes
    ----------

    var_in: PortInputengine
        The input Port of the system
    ex: PortExhaust
        The output Port of the system
    deltat
        The timestep of the simulation
    N_ref
        Rotation speed reference
    lamb_ref
        Stochiometric reference ratio
    delta_phi_comb_ref
        The reference duration of the combustion (expressed in rotation angle of the crankshaft)
    m_ivc_ref
        Reference mass at inlet valve closing
    phi_id_ref
        The ignition delay reference (expressed in rotation angle of the crankshaft)
    Vc
        The clearance volume
    lhv
        The lower heating value
    Tw
        The cylinder's wall's temperature
    ma_ini
        The initial mass of air
    mWieberef
        One of the Wiebe curves used for combustion simulation
    acomb
        One of the Wiebe curves used for combustion simulation
    bcomb
        One of the Wiebe curves used for combustion simulation
    am
        One of the Wiebe curves used for combustion simulation
    bm
        One of the Wiebe curves used for combustion simulation
    cm
        One of the Wiebe curves used for combustion simulation
    a_id
        One of the Wiebe curves used for combustion simulation
    b_id
        One of the Wiebe curves used for combustion simulation
    c_id
        One of the Wiebe curves used for combustion simulation
    delta_m
        One of the Wiebe curves used for combustion simulation
    mWiebe
        One of the Wiebe curves used for combustion simulation
    stroke
        The stroke length
    rc
        The compression ratio
    Tref
        The reference temperature used or energy calculation
    Sw
        The wall surface of exchange
    re
        The gaz state constant divided by the exhaust gas molar mass
    ra
        The gaz state constant divided by the air molar mass
    rf
        The gaz state constant divided by the fuel molar mass
    l
        The length of rod lenght

    ar
        The crank radius

    Cstoich
        The stochiometric coefficient for the modeled combustion reaction
    R
        The gaz constant (maybe we can remove it)
    test
        This is a variable
    phi_soi
        The crank angle at which the injection is started
    bore
        The bore of the engine's cylinder
    kf0
        The friction mean effective pressure coefficient
    kf1
        The friction mean effective pressure coefficient
    kf2
        The friction mean effective pressure coefficient
    rev_cy
        The number of revolutions per cycle
        
    topen
        timing for injection model: delay between soi and maximum mass fuel flow
    tclose
        timing for injection model: delay between  maximum mass fuel flow and going back to 0g/s
    SMFR
        maximum fuel flow during injection


    """
    def setup(self):
        """Cosapp setup of the model

        Populates cosapp's system of engine's 0d model. It uses the PortExhaust and the PortInputengine
        classes as ports, and the other variables as cosapp's inwards.

        """
        self.add_input(PortInputengine, 'var_in')    # define a new input port
        self.add_output(PortExhaust, 'ex')  # define a new output port
        # Solitary variables can also be added,
        # as either `inward` or `outward` variables:
        self.add_inward('deltat')
        self.add_inward('N_ref')
        self.add_inward('lamb_ref')
        
        self.add_inward('delta_phi_comb_ref')
        self.add_inward('m_ivc_ref')
   
        self.add_inward('phi_id_ref')
        self.add_inward('Vc')
        self.add_inward('mWieberef')
        self.add_inward('lhv')
 
        
        self.add_inward('Tw')
  
        self.add_inward('ma_ini')
        
        self.add_inward('a_id')
        self.add_inward('b_id')
        self.add_inward('c_id')
        self.add_inward('stroke')
        
        self.add_inward('rc')
        self.add_inward('acomb')
        self.add_inward('bcomb')
        self.add_inward('am')
        
        self.add_inward('bm')
        self.add_inward('cm')
        self.add_inward('Tref')
        self.add_inward('delta_m')
        
        self.add_inward('Sw')
        self.add_inward('re')
        self.add_inward('ra')
        self.add_inward('rf')
        
        self.add_inward('l')
        self.add_inward('ar')
        self.add_inward('Cstoich')
        self.add_inward('R')
        self.add_inward('test')


        self.add_inward('phi_soi')
        self.add_inward('bore')
        
        
        self.add_inward('kf0')
        self.add_inward('kf1')
        self.add_inward('kf2')
        self.add_inward('rev_cy')
     
        self.add_inward('topen')
        self.add_inward('tclose')
        self.add_inward('SMFR')
     

        
        self.add_inward('T_Test', dtype= np.ndarray,value=np.array([]))
    
    
     
      
        
        
     
    
    def compute(self): # `compute` method defines what the system does
         
         """The compute method to run the simulation.

          It sets up the Close_cycle class to run the simulation (we could probably optimize it with a bit of OOP)

          """
  
        

      #  print('yuki pression ini',self.var_in.p_IR, self.var_in.T_IR)
  
         load = "25%"
         combustion_model(load)
         file_paths = MainProgram.pages[0].return_file_paths()
         trial_path = file_paths[10]
         df = pd.read_csv(rf"{trial_path}")
         
         
         print(len(df))
         
         flag=-1
         flag2=-1
         mark=-1
         
         for i in range(len(df)):
          if df['CA (deg)'][i]<0 and flag == -1: 
              flag = 1 
              mark = i
          if    df['CA (deg)'][i]>=179 and flag == 1 and flag2==-1:  
              flag2=1
              mark2=i
              
         df = df[mark:mark2+2]  
         
       
      

         self.ex.p_exh =df['P(bar)'].to_numpy()
         self.ex.T_exh = df['T(K)'].to_numpy()
         self.ex.V_out= df['V (m3)'].to_numpy()
  
         self.ex.Qf = df['Heat Release (MW)'].to_numpy()
         self.ex.Qw_dot = df['heat_loss [MW]'].sum()
         
         time = df['time (s)'].to_numpy()
         
         diff_time = time[-1] - time[0]
         
         W = trapz(df['dWv_dt (Vmean)'], time)
         
         power= W*12 / diff_time / 1000

         
        
    
         print('work done',power)
    
         self.ex.W_dot = power
      
       
        
         #InCylinder_cycle.checkresults()    
        
        #self.ex.W_dot =4900000
        
     
        

        #○check if revolution per cycle
   
