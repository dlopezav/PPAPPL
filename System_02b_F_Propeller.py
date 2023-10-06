# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 17:25:23 2022

@author: hbusson
"""
#TODO
#gearbox, motor and wageningen
# formula for ms and mhydro to be checked

# Model computed with the help of the paper Ship Propulsion Plant Transient Response Investigation using a Mean Value Engine Model #
# by G. P. Theotokatos #

from math import *
import math
import pandas as pd
from scipy.interpolate import CubicSpline

class Propeller:
    """This is the propeller model.

    Attributes
    ----------

    rho_sw: float
        Seawater density
    D_p: float
        Propeller diameter
    omega: float
        The ship wake fraction
    V_s: float
    ms: float * sh_displ
        rho_sw * sh_displ
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
    p_D_p: float
        The pitch to diameter ratio
    A_e_A_o: float
        Disk area coefficient
    z_p: float
        Propeller's blade number
    R_w: float
        Wind resistance
    R_awl
        Water resistance
    
    """
  
    def __init__(self,V_s, NE, z_p,A_e_A_o, p_D_p,  D_p, omega,  
                 rho_sw, sh_displ, m_hydro, c_1,thrust,
                 cq, ct,s2,t2,u2,v2,
                  s1,t1,u1,v1,
                  n1,n2,  time_t):

          """Builds a Propeller object"""
   
          self.rho_sw= rho_sw
          self.D_p = D_p
          self.omega = omega
          self.V_s = V_s
   
        
          self.V_s = V_s
        #formula for ms and mhydro to be checked
          self.ms = rho_sw* sh_displ
          self.m_hydro = m_hydro
        
          self.c_1 = c_1
          self.thrust = thrust

        
          self.NE = NE
          self.ct = ct
          self.cq = cq
          self.s1 = s1
          self.t1 = t1
          self.u1 = u1
          self.v1 = v1
          self.s2 = s2
          self.t2 = t2
          self.u2 = u2
          self.v2 = v2
          self.n1 =n1.astype(int)
          self.n2 = n2.astype(int)
          self.p_D_p = p_D_p
          self.A_e_A_o = A_e_A_o
          self.z_p = z_p
          #wave and wind resistance initialised at 0

          self.time_t = time_t
          #
        
          
          
        
    def Gearbox(self):
        """Gearbox reduction. To be implemented."""
     
        self.Np = self.NE

         
    def Wageningen(self):
        """Compute the k_T and k_Q values for the propeller."""
        self.k_T = 0
        self.k_Q = 0
 
        for i in range(max(self.n1)):
            
            self.k_T = self.k_T + self.ct[i]*pow(self.J,self.s1[i])* pow(self.p_D_p, self.t1[i])*  pow(self.A_e_A_o, self.u1[i])*  pow(self.z_p, self.v1[i])
            
        for j in range(max(self.n2)): 
            self.k_Q = self.k_Q + self.cq[j]*pow(self.J,self.s2[j])* pow(self.p_D_p, self.t2[j])*  pow(self.A_e_A_o, self.u2[j])*  pow(self.z_p, self.v2[j])
        self.Reynolds_number()
        if self.Rey > 2*pow(10,6): 
            self.k_T=self.k_T +     self.delta_kt
            self.k_Q=self.k_Q +     self.delta_kq
            
 
     
    def Reynolds_number(self):
        """Compute the reynolds number"""
        C075R = 2.073 *self.D_p * self.A_e_A_o/self.z_p
        nu = 1.188 * pow(10,-6)
        
        
        n = self.Np /(2*math.pi)
        
        self.Rey= C075R * math.sqrt( pow(self.V_A,2)+ pow(0.75*math.pi*n*self.D_p,2))/nu
        
        self.L = math.log10( self.Rey ) -0.301
        
        self.delta_kt =0
        self.delta_kq =0
        
        df2 = pd.read_excel (r'C:/Users/hbusson/Documents/TNTM/PythonScripts/specific ship code/champs_elysees/wageningen_reynolds.xlsx')

     
        n1 = df2['n1'].to_numpy()
        n1 = n1.astype(int)
        CTn = df2['CTn'].to_numpy()
        sn = df2['sn'].to_numpy()
        tn = df2['tn'].to_numpy()
        un = df2['un'].to_numpy()
        vn = df2['vn'].to_numpy()
        wn = df2['wn'].to_numpy()


        n2 = df2['n2'].to_numpy()
        n2 = n2.astype(int)
        CQn = df2['CQn'].to_numpy()
        sn2 = df2['sn2'].to_numpy()
        tn2 = df2['tn2'].to_numpy()
        un2 = df2['un2'].to_numpy()
        vn2 = df2['vn2'].to_numpy()
        wn2 = df2['wn2'].to_numpy()
        
        for i in range(9):
            self.delta_kt = self.delta_kt+ CTn[i]*pow(self.L,sn[i])*pow(self.A_e_A_o,tn[i])*pow(self.p_D_p,un[i])*pow(self.J,vn[i])*pow(self.z_p,wn[i])
            
        for j in range(13):
            self.delta_kq = self.delta_kq+CQn[j]*pow(self.L,sn2[j])*pow(self.A_e_A_o,tn2[j])*pow(self.p_D_p,un2[j])*pow(self.J,vn2[j])*pow(self.z_p,wn2[j])
        
    
    def reading_file(self):
        
        self.kt = pd.read_excel(r'C:/Users/hbusson/Documents/TNTM/ship data/champs_elysees/BV/first delivery/KT.xlsx')
        self.kq = pd.read_excel(r'C:/Users/hbusson/Documents/TNTM/ship data/champs_elysees/BV/first delivery/KQ.xlsx')



        self.spl_kt = CubicSpline(self.kt['J'], self.kt['KT'])
        self.spl_kq = CubicSpline(self.kq['J'], self.kq['KQ'])
        steadystate= pd.read_excel(r'C:/Users/hbusson/Documents/TNTM/PythonScripts/specific ship code/champs_elysees/new version/steady_states.xlsx')

        #order of speed 

    



        #wind resistance and wave resistance
        self.rwi= steadystate['Rwind']
        self.rwa =steadystate['Rawl']

        splrwi = CubicSpline(steadystate['x'], steadystate['Rwind'])  
        splrwa = CubicSpline(steadystate['x'], steadystate['Rawl'])
        
        self.R_w = splrwi(self.time_t)
        self.R_awl = splrwa(self.time_t)


    
    def prop_comput(self):
        """Compute propulsion values."""
        
        self.Gearbox()
        
   
        
        self.Lwl = 406.08
        self.Bwl =61.3
        #ship afward draft
        self.Draught=16
        
        displ_vol=212404
        self.Cb = displ_vol /(self.Lwl * self.Bwl * self.Draught)
        
        # torque and thrust computed thanks to k_T and k_Q propelleries series
        self.omega = 0.5*self.Cb-0.05
       
        self.V_A = (1-self.omega)* self.V_s
        n = self.Np/(60) 
        self.n =n
        self.J = abs(self.V_A/(n*self.D_p))
        self.J=0.2
        
        #print('J',self.J,'V_A', self.V_A)
        
        self.Wageningen()
        
        self.reading_file()
        
        self.k_T = self.spl_kt(self.J)
        self.k_Q = self.spl_kq(self.J)
        
        
        self.eta_p = self.k_T * self.J/(self.k_Q*2*math.pi)
        
        
        self.Qp = self.k_Q* self.rho_sw * pow(n,2)* pow(self.D_p,5)
        self.Tp = self.k_T* self.rho_sw* pow(n,2)* pow(self.D_p,4)
        

        

        
    def Motor(self):
        #to be modified
        self.Tm = self.Qp
        self.Power_Motor = self.Tm*self.Np


    def Velocity_derivative(self):
        """Compute the ship acceleration"""
        #computation of the ship speed derivative
        
        self.rt =  420.33*pow(self.V_s,2) + 712.74*self.V_s + 12158
        
        #ship breadth 
        self.B =61.3
        
        #ship afward draft
        self.TA=16
        
        #propeller diameter
        self.Dp =10
        #Longitudinal center of buoyancy
        self.lcb =0.5
        #ship lenght waterline
        self.L = 406.08
        #midway between the foremost and the aftmost perpendiculars
        self.AM = self.B*self.TA
        #ship lenght waterline
        self.Lwl = 406.08
        displ_vol=212404
        self.Cp = displ_vol / (self.AM * self.Lwl)
        #  the single screw after body form
        self.Cstern = -10
        
        self.thrust = 0.25014*pow(self.B/self.L,0.28956)*pow((math.sqrt(self.B/self.TA)/self.Dp),0.2624)
        self.thrust= self.thrust / pow((1-self.Cp+0.225*self.lcb),0.1762) +0.0015*self.Cstern
    
     
   
        self.dVs_dt = (self.Tp*(1-self.thrust) - self.rt )/(self.ms+ self.m_hydro) 
        return self.dVs_dt
    
 
    
    def checkresult(self):
        """Check the values obtained by printing them"""
        print('CHECK RESULT PROPELLER')
        
        print('power motor here', self.Power_Motor ,'Tm', self.Tm,'Np',self.Np )
        print('checking velocity  dVs',self.dVs_dt,'Vs', self.V_s, (self.Tp*(1-self.thrust) - self.rt ),(self.ms+ self.m_hydro) )
        print('checking qp', self.Qp, 'k_Q', self.k_Q, 'np', self.Np, 'J', self.J,'va' ,self.V_A,'n',self.n, 'dp',self.D_p, self.V_s)
        
        print('checking tp', self.Tp,'k_T',  self.k_T, self.J)
        
        print('thrust', self.thrust, 'wake fraction', self.omega)
        
        print('J', self.J, self.n, self.Dp, self.V_A, self.ms)
        
        
        print('thrust', self.thrust)
        
        print('qp', self.Qp , 'k_Q', self.k_Q,'rho sw',self.rho_sw , 'n', self.n, 'D_p',self.D_p)
    
    
 #input = Vs, NE
#output = dvs/dt, Tm    
        