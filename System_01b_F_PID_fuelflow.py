# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 16:36:31 2022

@author: hbusson
"""

import math

# Model computed with the help of the paper Development of a combined mean value-zero dimensional model and
# application for a large marine four-stroke Diesel engine simulation #
# by Baldi et al #

#Problem for determining mass air flow #

class Current:
    """Stores a current state.

    Attributes
    ----------

    NE_Current: float
        Engine speed
    NTC_Current: float
        Turbo charger speed
    Vs_Current: float
        This is a value
    """
    def __init__(self, NE_Current, NTC_Current, Vs_Current):

        """Builds a Current object"""
        #initial conditions
        self.NE_Current = NE_Current
        self.NTC_Current = NTC_Current
        self.Vs_Current = Vs_Current
        
    def update(self,NE, NTC, Vs):
        """Updates a Current object"""
        self.NE_Current = NE
        self.NTC_Current = NTC
        self.Vs_Current = Vs


class Controller_fuel:
    """Holds the computation for the fuel flow control

    Attributes
    ----------

    xr_inter: float
        Intermediary time injection (computed at object build)
    t: float
        Timestep for the computation
    xr_o: float
        Original time of injection
    kp: float
        PID coefficient
    ki: float
        PID coefficient
    kd: float
        PID coefficient
    NE: np.ndarray
        Engine's rotation speed
    Nord: np.ndarray
        Engine's rotation speed order
    NE_Current: float
        Current engine's rotation speed
    zc: int
        Engine's number of cylinder
    rev_cy: float
        Revolutions per cycle
    count: int
        Current number of iteration of the simulation
    """
    def __init__(self,NE_Current, kp,ki,kd, NE,Nord , tinj_o,t, rev_cy, zc,count):

        """Builds the object"""
    

        
        self.tinj_inter = tinj_o+ kp*60/(2*math.pi)*(NE[count]-Nord[count])
        self.t = t
        self.tinj_o = tinj_o
     
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.NE = 60/(2*math.pi)*NE
        self.Nord =60/(2*math.pi)* Nord
        self.NE_Current = NE_Current
       

        self.zc = zc
        self.rev_cy = rev_cy
        self.count = count

        
        
        
    def Integ(self):
        
        """Computes the integrative part for the control"""
        inte=0
      
       
        for i in range(1,self.count+1):
           
             b= self.Nord[i-1]-self.NE[i-1]
             B =self.Nord[i]-self.NE[i]
             inte= inte+( b+B)/2
           
        return inte
    
    def Deriv(self):
        """Computes the derivative part for the control"""
        count= self.count
        der=0
        if count>1:
            der = -(self.NE[count]-self.Nord[count]) + (self.NE[count-1]-self.Nord[count-1])
      
        return der
            
    def final_position(self):
        """Computes, sets and return the final time of injection"""
        self.t_inj = self.tinj_inter +self.kd*self.Deriv()+self.ki*self.Integ()
        #print('final position', 'inter', self.xr_inter , self.kd, 'deriv',self.Deriv(),self.ki, 'integ',self.Integ())
        #print('NE',self.NE[self.count], 'Nord', self.Nord[self.count])
        return self.t_inj

 
    
    def Mass_fuel_cycle(self):
    
        """Computes the mass of fuel for a cycle"""
        mf_cy = 14500*self.t_inj
        
        return mf_cy
    
    
    def running(self):

          """Run a step of the control loop."""
          self.t_inj = self.final_position()

          self.mf_cy = self.Mass_fuel_cycle()
       

          self.mdot_f = self.zc*self.mf_cy*self.NE_Current/(60*self.rev_cy)
          
         
         
    def checkresult(self):
         
          print('CHECK RESULT PPID')
          print(' t_inj',self.t_inj, 'self.mf_cy', self.mf_cy,'self.mdot_f ', self.mdot_f )
          print('tinj inter', self.tinj_inter,self.tinj_o)
    
          