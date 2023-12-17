# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 17:33:56 2022

@author: hbusson
"""
from math import *
import math
import pandas as pd
from SAC_PLS import *
from main import *

#TODO
#calcul heat exchange in exhaust receiver

class Enginederivatives:
    """This class is used to compute the full cycle of the engine.


    Attributes
    ----------

    p_a : float
        The ambiant pressure
    k_ep : float
        coefficient for computing The pressure after turbine
    R_a :
         R divided by air mass
    gamma_e :
        Ratio of specific heat for the exhaust gas
    gamma_a :
        Ratio of specific heat for the air
    V_D :
        Displacement Volume
    eta_t :
        Turbine efficiency
    eta_comb :
        Combustion efficiency
    eta_sh :
        Shaft efficiency
    I_e :
        Engine inertia
    I_sh :
        Shaft inertia
    I_p :
        Propeller inertia
    rev_cy :
        Revolution per cycle
    mdot_f :
        Mass flow for something
    ITC :
        Turbo charger inertia
    cpa :
        Air's heat capacity
    cpe :
        Exhaust gas' heat capacity
    deltat :
        Timestep for the simulation
    mdot_t :
        Turbine's mass flow
    T_a :
        Ambiant air's temperature
    Q_p :
        Torque of the propeller transmitted to the engine
    A_eq :
        Equivalent cylinders flow area (Aeq)
    W_dot: float
        The derivative of the work
    Qw_dot: float
        The derivation of the heat transfer
    cd :
        Discharge value
    lhv :
        Low heating value
    T_exh :
        Exhaust gas temperature
    p_exh :
        Exhaust gas pressure
    T_IR :
        Inlet receiver temperature
    p_IR :
        Inlet receiver pressure
    mdot_c :
        Compressor's mass flow
    kpac: float
        coefficient to determine pressure drop in air cooler as a function of mass flow
    kpaf: float
        coefficient to determine pressure drop in air filter as a function of mass flow
    cpc :
        heat capacity for working medium in the compressor
    cpt :
        heat capacity for working medium in the turbine
    eta_c :
        Charger efficiency
    cvir: float
        Inlet receiver Cv
    cver: float
        Exhaust receiver Cv
    V_ir: float
        Inlet receiver volume
    V_er: float
        Exhaust receiver volume
    NE_Current : float
        Engine speed
   kf0: float
       friction pressure as a function of NE polynome's coefficient 
   kf1: float
      friction pressure as a function of NE polynome's coefficient 
   kf2: float
       friction pressure as a function of NE polynome's coefficient 
    k_ac0: float
        Air cooler efficiency's polynome's coefficient
    k_ac1: float
        Air cooler efficiency's polynome's coefficient
    k_ac2: float
        Air cooler efficiency's polynome's coefficient
    T_w :
        Sea water temperature
    prc :
        pressure ratio compressor
    prt :
        pressure ratio turbine
    """
    def __init__(self,p_a, k_ep,  gamma_e, eta_sh, eta_t, eta_c, I_e, I_sh, I_p, rev_cy, V_D,
                        ITC, cpa, cpe,  Q_p,  mdot_f,mdot_c, mdot_t, T_a,
                        A_eq, gamma_a, p_IR, T_IR, W_dot, Qw_dot, eta_comb, lhv, p_exh,
                        kpac, kpaf, cpc, cpt, cver, cvir, V_ir, V_er, cd, R_a, NE_Current, 
                        kf0, kf1 ,kf2,
                        k_ac0, k_ac1, k_ac2, T_w, T_exh,
                        deltat, NTC, prc, prt):

       self.p_a = p_a
       self.k_ep = k_ep
       self.R_a = R_a
       self.gamma_e = gamma_e
       self.gamma_a = gamma_a
       self.V_D = V_D
       self.eta_t = eta_t
       #function of equivalence ratio
       self.eta_comb = eta_comb
       self.eta_sh = eta_sh
       self.I_e = I_e
       self.I_sh = I_sh
       self.I_p = I_p
       self.rev_cy = rev_cy
       self.mdot_f = mdot_f
       self.ITC = ITC
       self.cpa = cpa
       self.cpe = cpe
       self.deltat = deltat     
       self.mdot_t = mdot_t
       self.T_a = T_a
       self.Q_p = Q_p      
       self.A_eq = A_eq
       self.W_dot = W_dot
       self.Qw_dot = Qw_dot
       self.cd = cd
       self.lhv = lhv
       self.T_exh = T_exh
       self.p_exh = p_exh
       self.T_IR = T_IR
       self.p_IR = p_IR
       self.mdot_c = mdot_c
       self.k_ep = k_ep
       self.kpac = kpac
       self.kpaf = kpaf
     
       self.cpc = cpc
       self.cpt =cpt
       self.eta_c = eta_c
       self.cvir = cvir
       self.cver = cver
       self.V_ir = V_ir
       self.V_er = V_er
       self.NE_Current = NE_Current
       self.kf0 =kf0
       self.kf1=kf1
       self.kf2 = kf2
       self.k_ac0 =k_ac0
       self.k_ac1=k_ac1
       self.k_ac2 = k_ac2
       self.T_w = T_w
       self.NTC = NTC
       self.Tref = 0
       self.prt = prt
       self.prc = prc

 #rajouter calcul ma
 
    def Training_PLS(self, rpm_res, pr_c_res):
       
       
       file_paths = main_page.return_file_paths()
       r_path = file_paths[2]
       TC_path = file_paths[3]
       df_w =  pd.read_excel(rf"{r_path}")
       df_w = df_w['T_cooling_water']
       df_map = pd.read_excel(rf"{TC_path}")

       

       if isinstance(rpm_res, int): 
           n=1
       else:    
           n = rpm_res.shape[0]
           
       n2 = len(df_map)
       h1 = 4
       h2 = 3
       f1 = 0
       f2 = 0
       i = 1

       t_w_ac = np.array([ ])
       p_atm_res = np.array([ ])
       t_atm_res = np.array([ ])
       p_sc_res = np.array([ ])
       t_sc_res = np.array([ ])
       t_sc_1_res = np.array([ ])
       p_c_res = np.array([ ])
       t_c_res = np.array([ ])
       m_a_res = np.array([ ])
       m_ac_res = np.array([ ])
       m_ac_result = np.array([ ])
       epsi_ac_res = np.array([ ])
       epsi_ac_real_res = np.array([ ])
  

       t_w_ac_res = np.array(df_w)
       error_ep_res = np.array([ ])
       error_tsc_res = np.array([ ])
       
       arr1 = []
       arr2 = []
       arr1b = []
       arr2b = []
       
       while len(arr1) < n:
           arr1 = np.append(arr1, len(arr1))
           
       while len(arr2) < n:
           arr2 = np.append(arr2, len(arr2))

       while len(arr1b) < len(df_map):
           arr1b = np.append(arr1b, len(arr1b))
           
       while len(arr2b) < len(df_map):
           arr2b = np.append(arr2b, len(arr2b))
           
           
     
       
       while i <= h1:
           k = 0
           while k <= i:
               f1 += 1
               k += 1
           i += 1
       i = 1
       while i <= h2:
           k = 0
           while k <= i:
               f2 += 1
               k += 1
           i += 1
           
       X1 = np.zeros((n,f1))
      
       X2 = np.zeros((n,f2))
       X1b = np.zeros((n2,f1))
       X2b = np.zeros((n2,f2))
       
       rpm_map = np.array(df_map.iloc[:,0])
       pr_map = np.array(df_map.iloc[:,1])
       
       for s1 in arr1:
           s1 = int(s1)
           r = 0
           i = 1
           
         
           while i <= h1:
               k = 0
               while k <= i:
                 
                   X1[s1][r] = rpm_res**(i-k)*pr_c_res**k
                   k += 1
                   r += 1
               i += 1
         
       for s2 in arr2:
           s2 = int(s2)
           r = 0
           i = 1
           while i <= h2:
               k = 0
               while k <= i:
                   X2[s2][r] = rpm_res**(i-k)*pr_c_res**k
                   k += 1
                   r += 1
               i += 1       
                   
       for s1b in arr1b:
          s1b = int(s1b)
          r = 0
          i = 1
          while i <= h1:
              k = 0
              while k <= i:
                  X1b[s1b][r] = rpm_map[s1b]**(i-k)*pr_map[s1b]**k
                  k += 1
                  r += 1
              i += 1

       for s2b in arr2b:
           s2b = int(s2b)
           r = 0
           i = 1
           while i <= h2:
               k = 0
               while k <= i:
                   X2b[s2b][r] = rpm_map[s2b]**(i-k)*pr_map[s2b]**k
                   k += 1
                   r += 1
               i += 1

      




  
           
       

       y1 = np.array(df_map.iloc[:,3])
       y2 = np.array(df_map.iloc[:,4])
       Yb = np.array([y1, y2])
       Yb = Yb.T
       
      
       pls1 = PLSRegression(n_components = 14) #Last change at 8-9 components (size of X1)
       pls1.fit(X1b,Yb)
       pls2 = PLSRegression(n_components = 9) #Limit at 14 components (size of X2)
       pls2.fit(X2b,Yb)
    
       Y1 = pls1.predict(X1)
     
       m_dot_c = Y1[:,1]
       Y2 = pls2.predict(X2)
       eta_c = Y2[:,0]    
      
       return m_dot_c, eta_c
           
 
    def Mapping(self):
        
        

      self.mdot_c, self.eta_c = self.Training_PLS(self.NTC, self.prc)
 

      print('self.NTC, self.prc',self.NTC, self.prc)
      print('HEEEEEERE',self.mdot_c, self.eta_c)
      
     

      
      if self.prt <3:
          self.mdot_t = -3.7*pow(self.prt,4)+49.775*pow(self.prt,3)-246.253*pow(self.prt,2)+539.385*self.prt-160.629
          
      else:
          self.mdot_t = 280.965
    
      self.eta_t =   5.5693*pow(10,-9)*pow(self.NTC,2) + 0.026216271*pow(self.prt,2)
      


    
      print('ntc', self.NTC, 'mdot_c', self.mdot_c, 'etac', self.eta_c, self.NTC, self.prc)
      
      print('ntc', self.NTC, 'mdot_t', self.mdot_t, 'etat', self.eta_t)
      
    
      
          
    def Mass_flow(self):
         
      
          
    
       #determinator of air mass flow
       
       #self.mdot_a=  self.cd*self.A_eq*self.p_IR/ (math.sqrt(self.R_a* self.T_IR)) *math.sqrt( 2*self.gamma_a/(self.gamma_a-1)* pow((self.p_IR)/self.p_exh , 2/self.gamma_a) -pow( (self.p_IR)/self.p_exh, (self.gamma_a+1)/self.gamma_a)) 
     
       
       #conservation of mass to determine exhaust mass flow
       cdaeq=15000000
       prcy = self.p_exh/self.p_IR
      # f=self.p_IR/math.sqrt(self.R_a*self.T_IR) * math.sqrt( 2*self.gamma_a/(self.gamma_a-1)*( pow(prcy,2/self.gamma_a)- pow(prcy,(self.gamma_a+1)/self.gamma_a)))
      #self.mdot_a = cdaeq*f
       self.mdot_a = 270

       self.mdot_e = 1000*self.mdot_f + self.mdot_a
       
       
      
      
          
     
    def Combining_cycles(self):
        
        
           print('wdot' ,self.W_dot )
           print('qwdot',self.Qw_dot)
        
           self.He_dot = self.mdot_a*self.lhv*self.eta_comb + self.mdot_a*self.cpa*self.T_IR - self.W_dot - self.Qw_dot
           
           
           
           pumping_work =  self.V_D* (self.p_IR  -self.p_exh)*pow(10,5)
   
          
           self.Pw_indicated = (self.W_dot )
       
           NErpm = self.NE_Current*60/(2*math.pi)
           self.p_f = self.kf0+ self.kf1*NErpm + self.kf2* pow(NErpm,2)
           self.Pw_f =     0.3695 + 4.7853e-05*NErpm +       self.Pw_indicated*2.9877e-02 
           self.Pwb_barre = self.Pw_indicated - self.Pw_f
            
           pi = math.pi
     
            
         
           self.Qe= self.Pwb_barre*30/(self.NE_Current*math.pi)
           
           
           print('p_IR', self.p_IR, 'p_exh', self.p_exh)
           print('Pwb', self.Pwb_barre, 'Pw_indicated', self.Pw_indicated, 'wdot', self.W_dot , 'friction', self.Pw_f )
           print('NE', self.NE_Current, self.rev_cy, 'NErpm', NErpm)
           print('pressure indicated', self.Pw_indicated/(pow(10,5)*self.V_D*NErpm/60) ,'pf', self.p_f)
        
    
           
    def pressure_ratios(self):
        
      
        self.delta_pac = self.kpac*pow(self.mdot_c,2)
        self.delta_paf = self.kpaf*pow(self.mdot_c,2)
        
        
        self.prc = (self.p_IR + self.delta_pac) /(self.p_a - self.delta_paf)
        #print('delta_pac', self.delta_pac, 'p_ir',self.p_IR)
        
        #print( 'prc', self.prc,(self.p_IR + self.delta_pac) )
                
        self.delta_ep = self.k_ep*pow(self.mdot_e,2)
        self.prt = (self.p_a + self.delta_ep)/self.p_exh
        
       # print('self.delta_ep', self.delta_ep )
        
        
    def Temp_and_pressure(self):    
     
        
        self.T_c = self.T_a*(1+ (pow(self.prc,(self.gamma_a-1)/self.gamma_a)-1)/self.eta_c )
        
      
        
        
        self.epsi = self.k_ac0 + self.k_ac1*self.mdot_c + self.k_ac2*pow(self.mdot_c,2)
       
        self.T_AC = self.T_c - self.epsi * (self.T_c- self.T_w)
        
      #  print('epsi', self.epsi)


        self.Ttd= self.T_exh* (1- self.eta_t*(1-pow(self.prt,(self.gamma_e-1)/self.gamma_e)))
        
      #  print('T_c', self.T_c, 'prc', self.prc,'eta_c',self.eta_c, 'pir', self.p_IR, 'p_a', self.p_a)
        
        self.ptd = (self.p_a + self.delta_ep)
        
      #  print('Ttd', self.Ttd, 'T_exh',self.T_exh, 'ptd', self.ptd,'eta_t',self.eta_t, 'pexh', self.p_exh, 'delta_ep', delta_ep)
     
   
         
    def Rotation_derivatives(self,NTC):
         pi=math.pi
       
         #computation or turbine and compressor torque
         self.Q_t = (self.mdot_t)*self.cpt*(self.T_exh-self.Ttd)/(pi*NTC/30)
         self.Q_c = self.mdot_c * self.cpc*(self.T_c-self.T_a)/(pi*NTC/30) 
    
         #computation of NTC and NE derivatives
         self.dNTC_dt = 30/math.pi*(self.Q_t - self.Q_c)/self.ITC 
         
    
      
         self.dNE_dt = 30/math.pi*(self.eta_sh*self.Qe - self.Q_p)/(self.I_e + self.I_sh + self.I_p)
   
      #   print('TORQUE qe', self.Qe, 'qp', self.Q_p)
      #   print( 'Qt', self.Q_t, 'Q_c', self.Q_c, 'te', self.T_exh,'ttd',self.Ttd,'tc',self.T_c,'ta',self.T_a)
      #   print('dnc',   self.dNTC_dt )
      #   print( 'dNTC',   self.dNTC_dt)
      #   print('dNE', self.dNE_dt, 'Q_e',self.Qe,'Q_p', self.Q_p,'inertia',self.I_e + self.I_sh + self.I_p )
         return self.dNTC_dt, self.dNE_dt
     

     
    def Inlet_receiver(self):
         
        #variation in mass flow = mass flow in (compressor) - mass flow out (inlet receiver mdot_a)
         dmt = self.mdot_c - self.mdot_a
         
         #internal energy of inlet receiver
         u = self.cvir*(self.T_IR-self.Tref)
         #molar mass of inlet receiver
         M_a = 28.96
         #mass in inlet receiver using perfect gas law
         mir = self.p_IR*pow(10,5)*self.V_ir/(0.287*self.T_IR)
         #density in inlet receiver
         rho= mir/self.V_ir
        #variation in temperature = (variation in enthalpy - udmt)/ cv*mass
         self.dTir = (self.mdot_c*self.cpc*self.T_AC - self.mdot_a*self.cpa*self.T_IR -u*dmt)/(self.cvir*mir)
         
        # print('dtir', dTir,'tac', self.T_AC, 'epsi', self.epsi, 'mccpctac',self.mdot_c*self.cpc*self.T_AC, 'macpatir', self.mdot_a*self.cpa*self.T_IR, 'udm',u*dmt )
         #update of temperature
         self.T_IR = self.T_IR + self.deltat*self.dTir
         
       
         
         self.p_IR = self.T_IR* rho*0.287/pow(10,5)
         
       #  print('pir',self.p_IR)
         
       #  print('dTir', self.deltat*dTir, 'T_AC', self.T_AC, 'T_IR', self.T_IR, 'mdot_c', self.mdot_c, 'mdot_a', self.mdot_a  )
         
         self.mir = mir
         
         
         
    def Exhaust_receiver(self):
        #variation in mass flow = mass flow in (compressor) - mass flow out (inlet receiver mdot_a)
        dmt = self.mdot_e - self.mdot_t
        
        #molar mass of exhaust receiver
        M_e = 62.00089485
        #mass in exhaust receiver using perfect gas law
        mer = self.p_exh*pow(10,5)*self.V_er/(0.13*self.T_exh)
        #internal energy of exhaust receiver
        u=self.cver*(self.T_exh-self.Tref)
        #density in inlet receiver
        rho= mer/self.V_er
        
     
        
        #variation in temperature = (variation in enthalpy - udmt)/ cv*mass
        #TO ADD: heat exchange
        self.dTer = (self.mdot_e*self.cpe*self.T_exh - self.mdot_t*self.cpt*self.Ttd -u*dmt)/(self.cver*mer)
        #update of temperature
        
       
        
        self.T_exh = self.T_exh + self.deltat*self.dTer 
        
    
        
        self.p_exh = self.T_exh* rho*0.13/pow(10,5)
        
        self.mer = mer
        
       # print('dTer', dTer, 'me and Te',self.mdot_e, self.T_exh , 'mt and Ttd',self.mdot_t,self.Ttd ,u*dmt)
        

        
   
    def checkresult(self):
        
        print('TORQUE qt',self.Q_t,'qc',self.Q_c, 'qe',self.Qe,'Qp', self.Q_p)
        
 
        

        
        print('POWER pb barre', self.Pwb_barre, 'pw indicated',  self.Pw_indicated , 'power friction', self.Pw_f)
        print('brake pressure', self.Pwb_barre/(pow(10,5)*self.V_D*self.NE_Current/(2*math.pi)),       'pressure indicated', self.Pw_indicated/(pow(10,5)*self.V_D*self.NE_Current/(2*math.pi)) ,'pf', self.p_f)
        
        print( 'Temperature Tc', self.T_c, 'Tac', self.T_AC, 'Ttd', self.Ttd, 'Tir', self.T_IR, 'T_exh', self.T_exh)
        
        print('pressure',  self.prc  , 'ptd',  self.ptd, 'pir',  self.p_IR, 'p_exh',self.p_exh)
        
        print('mass flows ma', self.mdot_a, 'mc', self.mdot_c, 'mt', self.mdot_t)
        
        print('ntc', self.NTC, 'prc',self.prc)


        print('efficiencies epsi ac', self.epsi, 'etat', self.eta_t, 'eta_c', self.eta_c)
        
        print('derivatives inlet and exhaust', self.deltat*self.dTer ,self.deltat*self.dTir)
        
        
#input Qp, Qe, Tc, mdot_f, p_IR, T_IR,p_exh, Ne,NTC xr, W_dot, Qw_dot, Wtot
#output dnedt, dncdt, p_exh, p_IR, T_IR

   