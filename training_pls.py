# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 11:50:35 2023

@author: hbusson
"""
import numpy as np
from SAC_PLS import *


def Training_PLS( rpm_res, pr_c_res):
       
       
    
       df_w =  pd.read_excel(r'C:\Users\hbusson\Desktop\jokin\results.xlsx')
       df_w = df_w['T_cooling_water']
       df_map = pd.read_excel(r'C:\Users\hbusson\Desktop\jokin\TC_maps2.xlsx')

       

       #print('rpm_res',rpm_res)

       if isinstance(rpm_res, int) or isinstance(rpm_res, float): 
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
       
       print('pls',m_dot_c, eta_c)
      
       return m_dot_c[0], eta_c[0]
           