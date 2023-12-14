# -*- coding: utf-8 -*-
"""
Created on Wed May 11 09:27:23 2022

@author: hbusson
"""


from wave_functions import *
from scipy.interpolate import CubicSpline
from training_pls import *
from SAC_PLS import *
from scipy.optimize import minimize, rosen, rosen_der
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
import pandas as pd
import numpy as np
import sympy as sp
import math
phi = sp.Symbol('phi')


def Area_i(phi):

    A_i = 1
    return A_i


def Area_e(phi):
    A_e = 1
    return A_e


def Inte_calc(phi):

    A_i = Area_i(phi)
    A_e = Area_e(phi)

    y = A_i*A_e/(math.sqrt(pow(A_i, 2)*pow(A_e, 2)))

    return y


class Initialisation:

    def __init__(self, rev_cy, zc, NE_ini, NTC_ini, Vs_ini):
        # engine geometry
        self.rev_cy = rev_cy  # revolution per cycle
        self.zc = zc  # number of engine cylinders
        # initial conditions
        self.NE_ini = NE_ini
        self.NTC_ini = NTC_ini
        self.Vs_ini = Vs_ini

    def update(self, NE, NTC, Vs):

        self.NE_ini = NE
        self.NTC_ini = NTC
        self.Vs_ini = Vs


class Fuel_flow_rate:

    def __init__(self, mf_cy, rev_cy, I_init):

        self.mdot_f = I_init.zc*mf_cy*I_init.NE_ini/(60*I_init.rev_cy)


class Propeller:

    def __init__(self, V_s, NE, z_p, A_e_A_o, p_D_p, k_Q, k_T, D_p, omega,
                 rho_sw, sh_displ, m_hydro,  I_p, pi, Water_resist):
        self.NE = NE
        self.D_p = D_p
        self.rho_sw = rho_sw
        self.Res_vec = np.array([])
        self.Res_wind = np.array([])

        self.Water_resist = Water_resist

        self.V_A = (1-omega) * V_s/1.94
        if NE > 0:
            self.J = self.V_A/(NE*D_p/60)
        else:
            self.J = 0

        self.K_T_real = -0.0881798*pow(self.J, 2) - 0.385233*self.J + 0.498566
        self.K_Q_real = -0.0182731170 * \
            pow(self.J, 2) - 0.0352535981*self.J + 0.0638759
        self.eta_p = self.K_T_real * self.J/(self.K_Q_real*2*pi)
        print('J', self.J, 'vs', V_s, 'Ne', NE)
        # k_Q  and k_T to be actually computed
        self.Qp = self.K_Q_real * rho_sw * pow(NE/60, 2) * pow(D_p, 5)
        self.Tp = self.K_T_real * rho_sw * pow(NE/60, 2) * pow(D_p, 4)

        # do not know how to have Ipa

        self.V_s = V_s
        # ms and mhydro to be determined
        self.ms = rho_sw * sh_displ
        self.m_hydro = m_hydro
      

        print('Qp', self.Qp, 'kq', self.K_Q_real, 'j', self.J)
        print('rho', rho_sw, ' ne', NE, 'dp', D_p)
        print('Tp', self.Tp, self.K_T_real, self.J)

    def Wave_Resistance(self, Water_resist, alpha, Hs, Tp, Vs, method_str):

        A = 1

        if method_str == 'Stawave1':
            self.Rwave = self.Water_resist.Stawave1(Hs)
        else:
            self.Rwave = self.Water_resist.irregular_wave(
                A, Hs, Tp, Vs, alpha, method_str)

        return self.Rwave

    def wind_resistance(self,  Wind_Dir, Wind_Speed):
        draft = 16
        height = 65
        Loa = 420
        Af = (height - draft)*61.3
        Al = (height - draft)*406
        Cl = 180

        CDl = 0.922-0.507*Al/(Loa*self.B)-1.162*Cl/Loa

        if CDl < 0:
            CDl = 0.3

        V1 = self.V_s + Wind_Speed*np.cos(Wind_Dir)
        V2 = Wind_Speed*np.sin(Wind_Dir)

        Baw = np.arctan(V1/V2)
        Vwr = np.sqrt(np.power(V1, 2)+np.power(V2, 2))

        rho_air = 1.293
        Rwind = 0.5*rho_air*Af*CDl*np.power(Vwr, 2)*np.cos(Baw)
        print('wind res eval', CDl, Vwr, Baw)

        return Rwind

    def Velocity_derivative(self, V_s, alpha, Hs, Tp, dummy, method_str,  Wind_Dir, Wind_Speed):

        self.rt = 1000*(6.983189029*pow(self.V_s, 2) -
                        101.63*self.V_s + 784.31)
        # estimated with the sea trials
        self.rt = 1000*(10.883*pow(self.V_s, 2) - 208.45*self.V_s + 1697.58)

        # ship breadth
        self.B = 61.3

        # ship afward draft
        self.TA = 16

        # propeller diameter
        self.Dp = 10
        # Longitudinal center of buoyancy
        self.lcb = 0.5
        # ship lenght waterline
        self.L = 406.08
        # midway between the foremost and the aftmost perpendiculars
        self.AM = self.B*self.TA
        # ship lenght waterline
        self.Lwl = 406.08
        displ_vol = 212404
        self.Cp = displ_vol / (self.AM * self.Lwl)
        #  the single screw after body form
        self.Cstern = -10
        # amesim formula
        self.thrust = 0.25014*pow(self.B/self.L, 0.28956) * \
            pow((math.sqrt(self.B/self.TA)/self.Dp), 0.2624)
        self.thrust = self.thrust / \
            pow((1-self.Cp+0.225*self.lcb), 0.1762) + 0.0015*self.Cstern

        self.ms = 200000

        self.m_hydro = 100000
        self.thrust = 0.2

        self.Res_vec = np.append(self.Res_vec, self.Wave_Resistance(
            self.Water_resist, alpha, Hs, Tp, V_s, method_str))
        self.Res_wind = np.append(
            self.Res_vec,  self.wind_resistance(Wind_Dir, Wind_Speed))
        self.wind = self.wind_resistance(Wind_Dir, Wind_Speed)

        dVs_dt = 1/1.94*(self.Tp*(1-self.thrust) - self.rt-dummy *
                  self.Res_vec[-1] - self.wind)/(self.ms + self.m_hydro)

        print('dvs_dt', dVs_dt, 'Tp', self.Tp, 'Tp(1-t)', self.Tp *
              (1-self.thrust), 'rt', self.rt, 'masses', (self.ms + self.m_hydro))
        print('self.V_s', self.V_s)
        print('thrust', self.thrust)
        print('res_vect',   self.Res_vec[-1])
        return dVs_dt, self.Res_vec[-1],  self.Res_wind[-1]


class governor_modelling:

    def __init__(self, zc, NE_Current, kp, ki, kd, NE, Nord, xr_o, t):

        self.tinj_o = xr_o
        self.t_inj_p = 10*kp*60/(2*math.pi)*(Nord[t]-NE[t])
        self.count = t
        self.xr_o = xr_o
        self.kp = kp
        self.ki = ki
        self.kd = kd
        self.NE = NE
        self.NE_Current = NE_Current
        self.Nord = Nord
        self.zc = zc

    def Integ(self):
        """Computes the integrative part for the control"""
        inte = 0

        for i in range(1, self.count+1):

            b = self.Nord[i-1]-self.NE[i-1]
            B = self.Nord[i]-self.NE[i]
            inte = inte+(b+B)/2

        return inte

    def Deriv(self):
        """Computes the derivative part for the control"""
        count = self.count
        der = 0
        if count > 1:
            der = -(self.NE[count]-self.Nord[count]) + \
                (self.NE[count-1]-self.Nord[count-1])

        return der

    def final_position(self):
        """Computes, sets and return the final time of injection"""
        self.t_inj = self.tinj_o+self.t_inj_p + self.kd*self.Deriv()+self.ki * \
            self.Integ()

        if self.t_inj < 0:
            self.t_inj = 0
        print('final position', 'inj o', self.tinj_o, 'inj p', self.t_inj_p,
              self.kd, 'deriv', self.Deriv(), self.ki, 'integ', self.Integ())
        print(self.t_inj)
        print('NE', self.NE[self.count], 'Nord', self.Nord[self.count])
        return self.t_inj, self.tinj_o, self.t_inj_p, self.Deriv()

    def Mass_fuel_cycle(self):
        """Computes the mass of fuel for a cycle"""
        mf_cy = 14500*self.t_inj

        return mf_cy

    def running(self):
        """Run a step of the control loop."""
        self.t_inj, tinj_o, t_inj_p, t_inj_d = self.final_position()

        self.mf_cy = self.Mass_fuel_cycle()

        self.mdot_f = self.zc*self.mf_cy*self.NE_Current/(60*2)


class Engine_parameters:

    def __init__(self, k_ac0, k_ac1, k_ac2, T_w, f_ac, R_a, p_a, A_ac, c_d, gamma_a, NTC, kc1,
                 kc2, kc3,
                 T_a, eta_c, zc, NE, Q_p,  ITC, I_e, I_sh, I_p, eta_sh, cpe, cpa, kf0, kf1, kf2,
                 pimax, rev_cy, V_D, H_l, mdot_f,  R_e,
                 x_r):

        self.k_ac0 = k_ac0
        self.k_ac1 = k_ac1
        self.k_ac2 = k_ac2
        self.T_w = T_w
        self.f_ac = f_ac
        self.R_a = R_a
        self.p_a = p_a
        self.A_ac = A_ac
        self.c_d = c_d
        self.gamma_a = gamma_a
        self.NTC = NTC
        self.kc1 = kc1
        self.kc2 = kc2
        self.kc3 = kc3
        self.T_a = T_a
        self.eta_c = eta_c
        self.zc = zc
        self.NE = NE
        self.Q_p = Q_p
        self.ITC = ITC
        self.I_e = I_e
        self.I_sh = I_sh
        self.I_p = I_p
        self.eta_sh = eta_sh
        self.cpe = cpe
        self.cpa = cpa
        self.eta_comb = eta_comb
        self.x_r = x_r
        self.pimax = pimax
        self.kf0 = kf0
        self.kf1 = kf1
        self.kf2 = kf2
        self.rev_cy = rev_cy
        self.kz0 = kz0
        self.kz1 = kz1
        self.V_D = V_D
        self.H_l = H_l
        self.mdot_f = mdot_f
        self.gamma_e = gamma_e
    
        self.R_e = R_e
       

    def compressor(self):

        self.prc = self.kc1*pow(self.NTC, 2)+self.kc2*self.NTC + self.kc3
        self.T_c = self.T_a * \
            (1+(pow(self.prc, ((self.gamma_a-1)/self.gamma_a))-1)/self.eta_c)
        print('compressor')
        print('prc', self.prc, 'T_c', self.T_c, 'etac',
              self.eta_c, 'ta', self.T_a, self.gamma_a)

    def Tinlet(self, mdot_a):
        # do not know how to have mdot_a at this stage
        self.epsi = self.k_ac0 + self.k_ac1*mdot_a + self.k_ac2*pow(mdot_a, 2)

        self.T_inlet = self.T_c - self.epsi * (self.T_c - self.T_w)

        return self.T_inlet

    def Pinlet(self, mdot_a):
        self.delta_pac = self.f_ac*self.R_a*self.T_c * \
            pow(mdot_a, 2)/(2*self.prc*self.p_a*pow(self.A_ac, 2))
        self.p_inl = self.prc*self.p_a - self.delta_pac
        print('Pinlet')
        print(self.p_inl, 'delta pac', self.delta_pac, 'prc',
              self.prc, 'pa', self.p_a, 'tc', self.T_c)
        return self.p_inl

    def Chemical_energy(self):
        pi_barre = self.x_r/0.01*self.pimax*self.eta_comb
        pf_barre = self.kf0+self.kf1*self.NE + self.kf2*pi_barre
        self.pb_barre = pi_barre - pf_barre
        self.dzeta = self.kz0 + self.kz1*self.pb_barre

        print('chemical energy')
        print('self.eta_comb', self.eta_comb,
              'pimax', self.pimax, 'xr', self.x_r)
        print('pb barre', self.pb_barre, 'pi barre',
              pi_barre, 'pf barre', pf_barre)
        print('dzeta', self.dzeta)

    def Main_Equation_1(self, x):

        pi = math.pi
        mdot_a = x[0]
        p_exh = x[1]

        self.epsi = self.k_ac0 + self.k_ac1*mdot_a + self.k_ac2*pow(mdot_a, 2)

        print('self epsi', self.epsi)

        self.epsi = self.k_ac0 + self.k_ac1*5 + self.k_ac2*pow(5, 2)

        print('self epsi', self.epsi)

        if self.epsi < 0.5 and self.epsi > 1:
            y = 100000000000

        else:

            Integrals_phi = sp.integrate(Inte_calc(phi), (phi, 0, 2*pi))

            self.A_eq = self.zc / (2*pi)*Integrals_phi

        #self.delta_pac = self.f_ac*self.R_a*self.T_c*pow(mdot_a,2)/(2*self.prc*self.p_a*pow(self.A_ac,2))
            self.delta_pac = 0.1
            self.p_inl = self.prc*self.p_a - self.delta_pac

            p_exh = self.p_inl - 0.2

            self.T_inlet = self.T_c - self.epsi * (self.T_c - self.T_w)

            #print('pinl', self.p_inl,self.delta_pac, 'tinlet', self.T_inlet, self.T_w, self.epsi, self.T_c)

            pr = p_exh/self.p_inl

            self.c_d = 150

        # print('plouf',self.T_inlet)

        #print('plaf',pr,2*self.gamma_a/(self.gamma_a-1), 2*self.gamma_a/(self.gamma_a-1) *  pow(pr, 2/self.gamma_a) , pow(pr, (self.gamma_a+1)/self.gamma_a) )
            f = 2*self.gamma_a / \
                (self.gamma_a-1) * (pow(pr, 2/self.gamma_a) -
                                    pow(pr, (self.gamma_a+1)/self.gamma_a))

            f2 = pow(self.c_d*self.A_eq * self.p_inl /
                     (math.sqrt(self.R_a*self.T_inlet)), 2)*f
            print('f2', f2, 'f', f, 'pr', pr,
                  'p_exh', p_exh, 'pinl', self.p_inl)
            print('f2 mdot_a', mdot_a)
            y1 = pow(mdot_a, 2) - f2
            # y2= (mdot_a*self.cpa*self.T_inlet + self.eta_comb*self.mdot_f*self.H_l*self.dzeta)*self.eta_ex - (mdot_a+self.mdot_f)* self.cpe*T_exh

            #y3 = mdot_a + self.mdot_f - self.At_eff *p_exh/ (math.sqrt(self.R_e* T_exh)) * math.sqrt( 2*self.gamma_e/(self.gamma_e-1)*pow(2*self.gamma_e/(self.gamma_e-1),((self.gamma_e+1)/(self.gamma_e-1))))

            print('y1', y1)

            y = abs(y1)
            # +abs(y2)+abs(y3)

        return y

    def Main_Equation_3(self, x):

        # review how calibrated

        p_exh = x[0]
        T_exh = x[1]

        self.delta_pac = 0.1
        self.p_inl = self.prc*self.p_a - self.delta_pac

        p_exh = self.p_inl - 0.15
        self.At_eff = 400

        y3 = 1000*self.mdot_a + self.mdot_f - self.At_eff * p_exh / (math.sqrt(self.R_e * T_exh)) * math.sqrt(
            2*self.gamma_e/(self.gamma_e-1)*pow(2*self.gamma_e/(self.gamma_e-1), ((self.gamma_e+1)/(self.gamma_e-1))))
        y3 = abs(y3)

        #print('y3',y3, 1000*self.mdot_a + self.mdot_f, self.At_eff *p_exh/ (math.sqrt(self.R_e* T_exh)) * math.sqrt( 2*self.gamma_e/(self.gamma_e-1)*pow(2*self.gamma_e/(self.gamma_e-1),((self.gamma_e+1)/(self.gamma_e-1)))))

        return y3

    def Update(self, mdot_a, p_exh, T_exh):

        self.mdot_a = mdot_a
        self.p_exh = p_exh
        self.T_exh = T_exh

    def Rotation_derivatives(self, t, NTC, NE):

        # include ttd and ptd here
        self.f_ep = 4.20*10e-12
        ptd = self.p_a + self.f_ep*pow(1000*self.mdot_a+self.mdot_f, 2)
        print('pexh', p_exh, 'ptd', ptd, self.p_inl)
        self.prt = self.p_exh/ptd

        if self.prt < 1:
            self.prt = self.prc - 0.2

        print('self prt', self.prt, 'self.prc', self.prc)
        self.eta_t = 5.5693*pow(10, -9)*pow(self.NTC, 2) + \
            0.026216271*pow(self.prt, 2)
        print('self.eta_t', self.eta_t, 'self.T_exh',
              self.T_exh, 'ptd', ptd, 'self.p_exh', self.p_exh)
        Ttd = self.T_exh * \
            (1 - self.eta_t*(1-pow(ptd/self.p_exh, (self.gamma_e-1)/self.gamma_e)))
        Q_t = (self.mdot_a + self.mdot_f/1000) * \
            self.cpe*(self.T_exh-Ttd)/(math.pi*NTC/30)
        Q_c = self.mdot_a * self.cpa*(self.T_c-self.T_a)/(pi*NTC/30)
        pi_barre = self.x_r/0.01*self.pimax*self.eta_comb
        pf_barre = self.kf0+self.kf1*self.NE + self.kf2*pi_barre
        self.pb_barre = pi_barre - pf_barre

        print('pb barre', self.pb_barre, 'vd', self.V_D, 'pi*rev_cy',
              (2*math.pi*self.rev_cy), 'revcy', self.rev_cy)

        Q_eng = 100000*self.pb_barre * self.V_D / (2*math.pi*self.rev_cy)

        dNTC_dt = (Q_t - Q_c)/self.ITC
        self.eta_sh = 0.98
        print('prosper yoplaboum', 'qeng', self.eta_sh*Q_eng, 'self.Q_p',
              self.Q_p, 'inertia', self.I_e + self.I_sh + self.I_p)
        print(self.eta_sh*Q_eng, self.Q_p, self.I_e, self.I_sh, self.I_p)
        dNE_dt = (self.eta_sh*Q_eng - self.Q_p) / \
            (self.I_e + self.I_sh + self.I_p)

        print('Qt', Q_t, 'mdota + mdot f', (1000*self.mdot_a +
              self.mdot_f), 'cpe', self.cpe, 'texh-ttd', self.T_exh-Ttd)
        print('qC', Q_c, 'mdot_a', 1000*self.mdot_a, 'tc-ta',
              (self.T_c-self.T_a), 'tc', self.T_c, 'ta', self.T_a)

        return dNTC_dt, dNE_dt, Q_t, Q_c, Q_eng, self.prt, self.prc, self.T_exh, Ttd, self.T_c


class Diff_equation_shaft:

    def __init__(self, current_time, delta_t, I_ini):

        self.dNE = []
        self.dNTC = []
        self.dVs = []
        self.current_time = 0
        self.delta_t = delta_t
        self.NE = I_ini.NE_ini
        self.NTC = I_ini.NTC_ini
        self.V_s = I_ini.Vs_ini
        self.NE0 = I_ini.NE_ini
        self.NTC0 = I_ini.NTC_ini
        self.V_s0 = I_ini.Vs_ini

    def rungekutta4(f, yn, tn1, tn, args=()):

        h = tn1 - tn
        k1 = f(yn, tn, *args)
        k2 = f(yn + k1 * h / 2., tn + h / 2., *args)
        k3 = f(yn + k2 * h / 2., tn + h / 2., *args)
        k4 = f(yn + k3 * h, tn + h, *args)
        yn1 = yn + (h / 6.) * (k1 + 2*k2 + 2*k3 + k4)
        return yn1

    def update(self, dNE, dNTC, dVs, current_time):

        self.dNE.append(dNE)
        self.dNTC.append(dNTC)
        self.dVs.append(dVs)
        self.current_time = current_time + self.delta_t

    def Next_timestep(self):

        self.NE = self.NE + self.dNE[-1]

        if self.NE <= 2:
            self.NE = 2

        self.NTC = self.NTC + self.dNTC[-1]
        if self.NTC < 5000:
            self.NTC = 5000
        self.V_s = self.V_s + self.dVs[-1]

        if self.V_s < 0:
            self.V_s = 0


#------------------------------------------------- Algorithme -------------------------- #
# data
df_me = pd.read_excel(
    'C:/Users/hbusson/Documents/TNTM/ship data/champs_elysees/BV/second delivery/recapData_BV_comparison.xlsx', sheet_name='ME')
df_stw = pd.read_excel(
    'C:/Users/hbusson/Documents/TNTM/ship data/champs_elysees/BV/second delivery/recapData_BV_comparison.xlsx', sheet_name='STW')


df_wave = pd.read_excel(
    'C:/Users/hbusson/Documents/TNTM/ship data/champs_elysees/BV/second delivery/recapData_BV_comparison.xlsx', sheet_name='wave')
# interpolations


df_wind = pd.read_excel(
    'C:/Users/hbusson/Documents/TNTM/ship data/champs_elysees/BV/second delivery/recapData_BV_comparison.xlsx', sheet_name='wind')
# interpolations

for i in range(len(df_me)):

    if np.isnan(df_me['me_ne'][i]):
        df_me['me_ne'][i] = df_me['me_ne'][i-1]
        df_me['me_load'][i] = df_me['me_load'][i-1]
        df_me['me_shaftp'][i] = df_me['me_shaftp'][i-1]

    if np.isnan(df_me['me_load'][i]):
        df_me['me_load'][i] = df_me['me_load'][i-1]

    if df_me['me_ne'][i] <= 2:
        df_me['me_ne'][i] = 2

lmin = min(df_me['me_load'])
lmax = max(df_me['me_load'])

nemin = min(df_me['me_ne'])
nemax = max(df_me['me_ne'])

wtmin = min(df_wave['Tp'])
wtmax = max(df_wave['Tp'])

whmin = min(df_wave['Height'])
whmax = max(df_wave['Height'])

almin = min(df_wave['Dir'])
almax = max(df_wave['Dir'])

wdmin = min(df_wind['Wind_Dir'])
wdmax = max(df_wind['Wind_Dir'])

wsmin = min(df_wind['Wind_Speed'])
wsmax = max(df_wind['Wind_Speed'])


df_me['interp'] = df_me['min_cumul']/100
df_stw['interp'] = df_stw['min_cumul']/100
df_wave['interp'] = df_wave['min_cumul']/100
df_wind['interp'] = df_wind['min_cumul']/100


spl_ne = CubicSpline(df_me['interp'], df_me['me_ne'])
spl_load = CubicSpline(df_me['interp'], df_me['me_load'])

spl_tp = CubicSpline(df_wave['interp'], df_wave['Tp'])
spl_hei = CubicSpline(df_wave['interp'], df_wave['Height'])
spl_alph = CubicSpline(df_wave['interp'], df_wave['Dir'])

spl_wsp = CubicSpline(df_wind['interp'], df_wind['Wind_Speed'])
spl_wdir = CubicSpline(df_wind['interp'], df_wind['Wind_Dir'])


df_ne = np.clip(spl_ne(df_stw['interp']), nemin, nemax)
df_ne = df_ne
df_load = np.clip(spl_load(df_stw['interp']), lmin, lmax)

df_tp = np.clip(spl_tp(df_stw['interp']), wtmin, wtmax)
df_hei = np.clip(spl_hei(df_stw['interp']), whmin, whmax)

df_wdir = np.clip(spl_wdir(df_stw['interp']), wdmin, wdmax)
df_wsp = np.clip(spl_wsp(df_stw['interp']), wsmin, wsmax)

df_alpha = np.clip(spl_alph(df_stw['interp']), almin, almax)
df_alpha = 0.0174533*df_alpha
df_wdir = 0.0174533*df_wdir
print(df_alpha)

for i in range(len(df_alpha)):

    if df_alpha[i] > math.pi:
        df_alpha[i] = 2*math.pi - df_alpha[i]
    if df_wdir[i] > math.pi:
        df_wdir[i] = 2*math.pi - df_wdir[i]

df_power = df_load/100*64000000


for i in range(1028):

    if df_ne[14000+i] < 65:
        df_ne[14000+i] = 65


df_ne = np.array(df_ne[2000:14000])
df_power = np.array(df_power[2000:14000])
df_load = np.array(df_load[2000:14000])
stw = df_stw['STW']
df_stw = np.array(stw[2000:14000])
df_wdir = np.array(df_wdir[2000:14000])
df_wsp = np.array(df_wsp[2000:14000])
print(df_wave['Tp'])



Nord = np.zeros((len(df_ne)))

Nord = df_ne

Vs_real = np.zeros((len(df_ne)))


Vs_real = df_stw


Power_real_tot = np.zeros((len(df_ne)))


Power_real_tot = df_power


#--------------------------- engine subsystem -------------------------------#


# number of engine cylinders
zc = 12
# initial rotation speed engine
NE_ini = df_ne[0]
# turbo compressor speed
NTC_ini = 7000
# original engine rack position
xr_o = 66/14500
# engine rotation 


# the temperature of the air cooler coolant medium
T_w = 30+273
# friction factor air cooler
f_ac = 0.01
# gas constant, resistance air
R_a = 287.058
# pressure air
p_a = 1
# Area air cooler
A_ac = 1000
# discharge coefficient
c_d = 1
# ratio of specific heats air
gamma_a = 1004.68506 / (1004.68506-287)

# kc = constant for compressor pressure ratio estimated from the steady states
kc1 = 0.0000000834
kc2 = -0.0007303296
kc3 = 3.20287563
# air Temperature
T_a = 300
# compressor efficiency
eta_c = 0.81
# Different Intertia: turbocharger, engine, shaft, propeller

ITC = 6
I_e = 500000
I_sh = 1
I_p = 1
# shaft efficiency
eta_sh = 0.99
# heat capacities, air and exhaust
cpe = 1100
cpa = 1004.68506
# combustion efficiency
eta_comb = 0.7
# maximum indicated mean effective pressure of the engine
pimax = 20
# friction mean effective pressure is considered function of the indicated mean effective pressure and the engine crankshaft speed, nthese are the coefficients

kf0 = 0.3695
kf1 = 4.7853e-05
kf2 = 2.9877e-02
# The proportion of the chemical energy of the fuel contained in the exhaust gas is considered linear function of the engine mean effective pressure, these are the coefficients
kz0 = 1
kz1 = 1
# diametre cylindre [m]
bore = 0.92
# stroke
stroke = 3.47
# volume displacement
V_D = zc*math.pi * pow(bore, 2)*stroke/4
# lower heating value
H_l = 50
# ratio heat capacity exhaust
gamma_e = 1.375
# thrust coefficient
thrust = 0.22


# exhuast gas constant
R_e = R_a


# air cooler effectiveness is assumed to be a function of the air mass flow rate, these are the coefficients
k_ac0 = 0.6929756
k_ac1 = 0.01474564
k_ac2 = -0.00019314


#-----------------------------------PID*-------------------------- #

# parameters for PID controller
# proportionnal
kp = 5.5*pow(10, -6)
# integral
ki = 6*pow(10, -5)
# derivatives
kd = 4*pow(10, -7)


###################### Propeller system ########################################

# parametres de base
# the number of propeller blades, zp
z_p = 5
# the disk area coefficient, AE/Ao
A_e_A_o = 0.57
# the pitch to diameter ratio, p/Dp
p_D_p = 0.665
# ship displacement
sh_displ = 280607
# add virtual mass
m_hydro = 1
# initial ship speed
Vs_ini = df_stw[0]

# water resistance included or not
dummy = 1


# ship wake fraction, w, which is considered constant taking values in the range from 0.20 to 0.45 for ships with a single propeller
omega = 0.4
#




# density of sea water
rho_sw = 1026
# non-dimensional coefficients  calculated from polynomial equations for Wageningen B propeller series
# integrate this
k_Q = pd.read_excel(
    'C:/Users/hbusson/Documents/TNTM/ship data/champs_elysees/BV/first delivery/KQ.xlsx')
k_T = pd.read_excel(
    'C:/Users/hbusson/Documents/TNTM/ship data/champs_elysees/BV/first delivery/KT.xlsx')

########## wave resistance ##############

#lenght to the waterline
Lwl = 406.8
#diameter propeller
D_p = 10
#breadth of ship
B = 61.3
#lenght between perpendicular
Lpp = 393.90
#coefficients for water resistance
LR = 49.8565
LE = 100.9574
#draft aft
Ta = 16
#draft 
Tf = 16
#middle draft
T = 16
#
kyy = 0.25*Lpp
#block coefficient
Cb = 0.7

WR = water_resistance(Lwl, D_p, B, sh_displ, rho_sw,
                      Cb, Lpp, Ta, Tf, T, kyy, LR, LE)


# revolutions per cycle
rev_cy = 1


####################### additonal ################


#timestep of the simulation
delta_t = 1





pi = math.pi




I_init = Initialisation(rev_cy, zc, NE_ini, NTC_ini, Vs_ini)


DE_S = Diff_equation_shaft(0, delta_t, I_init)


NE_tot = [NE_ini]
NTC_tot = [NTC_ini]
Ttd_tot = []
Vs_tot = [Vs_ini]
Qp_tot = []
Tp_tot = []
Qe_tot = []
Qt_tot = []
Qc_tot = []
Qt_tot = []
prt_tot = []
prc_tot = []
Texh_tot = []
Tc_tot = []
Ttd_tot = []
Jstore = []
t_tot = []
t_tot2 = [0]
Power_tot = []
Power_tot3 = []
Res_vec = []
Res_windv = []

t_inj_tot = []
t_inj_p_tot = []


method_str = 'SNNM'


niter = 5000
for t in range(niter):
    current_time = t
    delta_t = 1
    print('time', t)
    print(t)

    # comment est calculé Nord ????
    Governor = governor_modelling(
        zc, NE_ini, kp, ki, kd, NE_tot, Nord, xr_o, t)

    print('Beginning')

    Governor.running()

    print('plaf')
    t_inj, tinj_o, t_inj_p, t_inj_d = Governor.final_position()
    print('plif', current_time)

    print('tinj', t_inj)

    t_inj_tot.append(t_inj)

    t_inj_p_tot.append(t_inj_p)

    mf_cy = Governor.Mass_fuel_cycle()

    print('t_inj', t_inj, 'mf_cy', mf_cy)

    mdot_f = Governor.mdot_f

    print('mdot_f', mdot_f)
    print('plif', I_init.NE_ini)
    Fuel_flow = Fuel_flow_rate(mf_cy, rev_cy, I_init)

    print('fuel flow', Fuel_flow)
    print('ne ini', I_init.NE_ini)

    Propeller_1 = Propeller(I_init.Vs_ini, I_init.NE_ini, z_p, A_e_A_o, p_D_p,  k_Q, k_T, D_p, omega,
                            rho_sw, sh_displ, m_hydro,  I_p, pi, WR)

    print('propeller1',    Propeller_1.Qp)
    print('Vs_ini',   I_init.Vs_ini)

    Jstore.append(Propeller_1.J)

    Engine_param = Engine_parameters(k_ac0, k_ac1, k_ac2, T_w, f_ac, R_a, p_a, A_ac, c_d, gamma_a,
                                     I_init.NTC_ini, kc1, kc2, kc3, T_a, eta_c, zc, I_init.NE_ini,
                                     Propeller_1.Qp, ITC, I_e, I_sh, I_p, eta_sh, cpe, cpa, kf0, kf1, kf2, pimax, rev_cy, V_D, H_l, Fuel_flow.mdot_f,  R_e,
                                      t_inj )

    Engine_param.compressor()
    Engine_param.Chemical_energy()

    print('check here', I_init.NTC_ini,   Engine_param.prc)

    mdot_a, eta_c = Training_PLS(I_init.NTC_ini,   Engine_param.prc)

    print('check', mdot_a, eta_c)

    #x0 = [5,1.5]
    #bnds = ((1, 100), (1.3,1.6))
    # print(bnds)
    #Sol= minimize(Engine_param.Main_Equation_1, x0,  method='Nelder-Mead', bounds=bnds)
    #print('Sol x')
    # print(Sol)
   # mdot_a = Sol.x[0]

    Engine_param.Update(mdot_a, 2, 500)
    x0 = [2, 500]
    bnds2 = ((1, 20), (300, 2500))

    Sol2 = minimize(Engine_param.Main_Equation_3, x0,
                    method='Nelder-Mead', bounds=bnds2)
    print('Sol2 x')
    print(Sol2)

    p_exh = Engine_param.p_inl - 0.15
    T_exh = Sol2.x[1]

    Engine_param.Update(mdot_a, p_exh, T_exh)

    print('mdot a ', mdot_a, 'pexh ', p_exh, 'T_exh', T_exh)

    dNTC, dNE, Qt, Qc, Qe, prt, prc, Texh, Ttd, Tc = Engine_param.Rotation_derivatives(
        t, I_init.NTC_ini, I_init.NE_ini)
    print('dNTC', dNTC, 'dNe', dNE)

    print(df_tp[0])
    dVs, Res_water, Res_wind = Propeller_1.Velocity_derivative(
        I_init.Vs_ini, df_alpha[t], df_hei[t], df_tp[t], dummy, method_str,  df_wdir[t], df_wsp[t])

    DE_S.update(dNE, dNTC, dVs, current_time)

    DE_S.Next_timestep()

    print('derivatives')
    print('dne', dNE, 'dntc', dNTC, 'dvs', dVs, 'time', current_time)

    I_init.update(DE_S.NE, DE_S.NTC, DE_S.V_s)

    NE_tot.append(I_init.NE_ini)
    NTC_tot.append(I_init.NTC_ini)

    Power_tot.append(I_init.NE_ini*Qe / 9.5488)
    Power_tot3.append(I_init.NE_ini * Propeller_1.Qp / 9.5488)
    Vs_tot.append(I_init.Vs_ini)

    print(Propeller_1.Qp)

    Qp_tot.append(Propeller_1.Qp)
    Tp_tot.append(Propeller_1.Tp)
    Qe_tot.append(Qe)
    Qt_tot.append(Qt)
    Qc_tot.append(Qc)
    prt_tot.append(prt)
    prc_tot.append(prc)
    Texh_tot.append(Texh)
    Tc_tot = np.append(Tc_tot, Tc)
    Ttd_tot = np.append(Ttd_tot, Ttd)
    t_tot.append(t)
    Res_vec.append(Res_water)
    Res_windv.append(Res_wind)


niter = 5000
real_time_tot = [i for i in range(niter)]
print('real_time_tot', real_time_tot)


print(NTC_tot)


plt.plot(t_tot, Jstore)
plt.title('J')
plt.show()


plt.plot(t_tot, Qp_tot, 'r')
plt.plot(t_tot, Qe_tot, 'b')
plt.title('qe and qp')
plt.show()
plt.plot(t_tot, Tp_tot)
plt.title('Tp')
plt.show()

plt.show()
plt.plot(t_tot, Qc_tot, 'b')
plt.title('qc')


plt.plot(t_tot, Qt_tot, 'r')
plt.title('qt')
plt.show()
plt.plot(t_tot, prc_tot, 'r')

plt.plot(t_tot, prt_tot, 'k')
plt.title('prt and prc')
plt.show()


del NE_tot[-1]

print(NE_tot)
print(Nord)


plt.plot(t_tot, NE_tot, 'r')
plt.plot(real_time_tot[0:real_time_tot[-1]+1],
         Nord[0:real_time_tot[-1]+1], 'b')
plt.title('NE')
plt.show()
plt.plot(t_tot[2:len(t_tot)], NTC_tot[2:t_tot[-1]+1], 'r')

plt.title('NTC')
plt.show()



plt.plot(t_tot, t_inj_tot)
plt.title('t inj')
plt.show()

del Vs_tot[-1]

plt.plot(t_tot, Vs_tot, 'r')
plt.plot(t_tot,df_stw[0:niter], 'b')
plt.title('speed')
plt.show()


Power_tot2 = [1.5*x for x in Power_tot]

plt.plot(t_tot[0:niter], Power_tot[0:niter], 'r', label='Simulated power')
plt.plot(t_tot[0:niter], Power_tot3[0:niter],
         'k', label='Simulated power shaft')
plt.plot(real_time_tot[0:niter], Power_real_tot[0:niter], 'b', label='Measure')

plt.title('Power')
plt.figlegend()
plt.xlabel('Minutes')
plt.ylim(0, 6.4e7)
plt.show()



plt.plot(t_tot, t_inj_p_tot)
plt.title('t inj p')
plt.show()

g = 9.81

Fr = np.zeros((len(t_tot)))
for i in range(len(t_tot)):
    Fr[i] = df_stw[i]/math.sqrt(g*406)



plt.plot(t_tot, Res_windv)
plt.title('wind resistance')
plt.show()
plt.plot(t_tot, Res_vec)
plt.title('water resistance')
plt.show()
