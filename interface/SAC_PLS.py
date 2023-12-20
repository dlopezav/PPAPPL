# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 10:55:59 2023

@author: jirciogaray
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math as m
import scipy as sc
from sklearn.cross_decomposition import PLSRegression
from main import main_program_instance as MainProgram

file_paths = MainProgram.pages[0].return_file_paths()
folder_path = file_paths[0]
# Caminho da pasta
folder_path = file_paths[0]

# Lista de valores numéricos
engine_values = [25, 30, 35, 40, 45, 50, 55, 60, 65, 70]

# Dicionário para armazenar os DataFrames
dfs = {}

# Lendo os arquivos Excel para cada valor numérico
for value in engine_values:
    file_path = f"{folder_path}\\steady_states_engine_{value}_2.xlsx"
    df_name = f"df_{value}"
    dfs[df_name] = pd.read_excel(file_path)

# Lendo os demais DataFrames
df_w = pd.read_excel(f"{folder_path}\\results.xlsx")
df_w = df_w['T_cooling_water']
df_map = pd.read_excel(f"{folder_path}\\TC_maps2.xlsx")


# Criando a lista df
df = [dfs[f"df_{value}"] for value in engine_values]

def objective(x, k_c_1, k_c_2, k_c_0):
 return k_c_1 * x + k_c_2 * x**2 + k_c_0

n = len(df)
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
pr_c_res = np.array([ ])
rpm_res = np.array([ ])
t_w_ac_res = np.array(df_w)
error_ep_res = np.array([ ])
error_tsc_res = np.array([ ])

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

arr1 = []
arr2 = []
arr1b = []
arr2b = []

while len(arr1) < len(df):
    arr1 = np.append(arr1, len(arr1))
    
while len(arr2) < len(df):
    arr2 = np.append(arr2, len(arr2))

while len(arr1b) < len(df_map):
    arr1b = np.append(arr1b, len(arr1b))
    
while len(arr2b) < len(df_map):
    arr2b = np.append(arr2b, len(arr2b))

rpm_map = np.array(df_map.iloc[:,0])
pr_map = np.array(df_map.iloc[:,1])


for df in df:
    p_atm =  np.average(df.iloc[:,43])
    p_sc =  np.average(df.iloc[:,38]*0.001)+p_atm
    p_c = p_sc + 0.01
    pr_c = p_c/p_atm
    t_atm = np.average(df.iloc[:,45]+273.15)                  #K
    rpm = np.average((df.iloc[:,34]+ df.iloc[:,35] + df.iloc[:,36])/3)
    t_sc = np.average(df.iloc[:,47]+273.15)
    
    p_atm_res = np.append(p_atm_res, p_atm)
    p_sc_res = np.append(p_sc_res, p_sc)
    p_c_res = np.append(p_c_res, p_c)
    pr_c_res = np.append(pr_c_res, pr_c)
    rpm_res = np.append(rpm_res, rpm)
    t_atm_res = np.append(t_atm_res, t_atm)
    t_sc_res = np.append(t_sc_res, t_sc)
    

for s1 in arr1:
    s1 = int(s1)
    r = 0
    i = 1
    while i <= h1:
        k = 0
        while k <= i:
            X1[s1][r] = rpm_res[s1]**(i-k)*pr_c_res[s1]**k
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
            X2[s2][r] = rpm_res[s2]**(i-k)*pr_c_res[s2]**k
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
m_a_res = Y1[:,1]
Y2 = pls2.predict(X2)
eta_c = Y2[:,0]
t_st = 298.15                                         # [K]
p_st = 100.0                                            # [kPa]
gamma_air = 1.4
i =0

for e in eta_c:
    t_c_res = np.append(t_c_res,t_atm_res[i] * (1 + (1/e)*((pr_c_res[i])**((gamma_air - 1)/gamma_air)-1)))
    m_ac_res = np.append(m_ac_res,m_a_res[i] * ((t_atm_res[i])/t_st)**(1/2)/(p_atm_res[i]/(p_st/100)))
    epsi_ac_res = np.append(epsi_ac_res, (t_c_res[i]-t_sc_res[i])/(t_c_res[i]-(t_w_ac_res[i]+273.15)))
    i += 1
  
m_ac_fit = m_ac_res
epsi_fit = epsi_ac_res
t_sc_fit = t_sc_res

k_count = [1,2,3,4,5,6,7,8,9]

for k in k_count:
    
    if m_ac_fit[k] < m_ac_fit[k-1]:
        m_ac_fit_change = m_ac_fit[k-1]
        m_ac_fit[k-1] = m_ac_fit[k]
        m_ac_fit[k] = m_ac_fit_change
       
        epsi_fit_change = epsi_fit[k-1]
        epsi_fit[k-1] = epsi_fit[k]
        epsi_fit[k] = epsi_fit_change
        
        t_sc_fit_change = t_sc_fit[k-1]
        t_sc_fit[k-1] = t_sc_fit[k]
        t_sc_fit[k] = t_sc_fit_change
        
   

x, y = m_ac_fit, epsi_fit
popt, _ = sc.optimize.curve_fit(objective, x, y)
k_1, k_2, k_0 = popt 
j = 0

for m_ac in m_ac_res:
    epsi_ac_real = k_0 + k_1 * m_ac + k_2 * m_ac ** 2
    t_sc_1 = t_c_res[j] - epsi_ac_real*(t_c_res[j] - (t_w_ac_res[j]+273.15)) 
    epsi_ac_real_res = np.append(epsi_ac_real_res,epsi_ac_real)   
    t_sc_1_res = np.append(t_sc_1_res,t_sc_1)
    j += 1
   
m_ac_fit2 = m_ac_res
epsi_real_fit = epsi_ac_real_res
t_sc1_fit = t_sc_1_res
error_ep_res = (abs(epsi_ac_real_res - epsi_ac_res)/epsi_ac_res)
error_tsc_res = (abs(t_sc_1_res-t_sc_res)/t_sc_res)     
for k in k_count:
    
    if m_ac_fit2[k] < m_ac_fit2[k-1]:
        m_ac_fit2_change = m_ac_fit2[k-1]
        m_ac_fit2[k-1] = m_ac_fit2[k]
        m_ac_fit2[k] = m_ac_fit2_change
       
        epsi_real_fit_change = epsi_real_fit[k-1]
        epsi_real_fit[k-1] = epsi_real_fit[k]
        epsi_real_fit[k] = epsi_real_fit_change
        
        t_sc1_fit_change = t_sc1_fit[k-1]
        t_sc1_fit[k-1] = t_sc1_fit[k]
        t_sc1_fit[k] = t_sc1_fit_change    


plt.plot(m_ac_fit, epsi_real_fit, '-r',label = 'Predicted')
plt.plot(m_ac_fit, epsi_fit, 'xk', label = 'real')
plt.title('Effectiveness')
plt.xlabel('Mass folw rate (kg/s)')
plt.legend()
plt.show()

plt.plot(m_ac_fit, t_sc_fit,'xk',label = 'Data') 
plt.plot(m_ac_fit, t_sc1_fit, 'or',label = 'Predicted')
plt.title('Scav Temperature')
plt.xlabel('Mass folw rate (kg/s)')
plt.legend()
plt.show()

case = ['25%', '30%', '35%', '40%','45%', '50%', '55%', '60%', '65%', '70%']
data_sacpls = {'CASE':case,
        'rpm': rpm_res,
        'P atm [bar]': p_atm_res,
        'T atm[K]': t_atm_res,
        'P compressor [bar]': p_c_res,
        'T compressor [K]': t_c_res,
        'Pressure ratio [-]': pr_c_res,
        'P scav [bar]': p_sc_res,
        'T scav [K]': t_sc_res,
        'T scav calc [K]': t_sc_1_res,
        'T scav error':error_tsc_res,
        'T cooling water [K]': t_w_ac_res,
        'Effectiveness': epsi_ac_res,
        'Efectiveness clalc': epsi_ac_real_res,
        'Effectiveness error': error_ep_res,
        'Mass flow [kg/s]': m_ac_res,
               }


df_results_sacpls = pd.DataFrame(data_sacpls)
df_results_sacpls['vol_flow'] = df_results_sacpls['Mass flow [kg/s]']/1.204
df_results_sacpls['u_c'] = (df_results_sacpls['rpm']/60)*m.pi*0.924
df_results_sacpls['Mach'] = df_results_sacpls['u_c']/(1.4*df_results_sacpls['T atm[K]']*287)**(1/2)
df_results_sacpls['phi'] = df_results_sacpls['vol_flow']/(df_results_sacpls['u_c']*0.670554)
df_results_sacpls['psi'] = (df_results_sacpls['Pressure ratio [-]']**(2/7)-1)/(0.4*df_results_sacpls['Mach']**2)
df_plot = pd.DataFrame()
df_plot['CASE'] = df_results_sacpls['CASE']
df_plot['phi'] = df_results_sacpls['phi']
df_plot['psi'] = df_results_sacpls['psi']
df_plot = df_plot.sort_values('phi', ascending=True)


plt.plot(df_plot['phi'],df_plot['psi'])


