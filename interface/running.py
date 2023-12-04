# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 15:37:58 2022

@author: hbusson
"""
#initial mass in cylinder : ideal gas law
#initial mass receivers : ideal gas law as we have tir pir for both


#libraries needed


import numpy as np
from scipy.interpolate import CubicSpline
from System_01a_IO_PID_fuelflow import *
from System_02a_IO_Propeller import *
from System_03a_IO_Powermanagement import *
from System_04a_IO_Engine0d import *
from System_05a_IO_EngineMVEM import *
from Init_system import *
from math import *
from main import MainProgram as MP

class Current:
    
    def __init__(self, NE_Current, NTC_Current, Vs_Current):


        #initial conditions
        self.NE_Current = NE_Current
        self.NTC_Current = NTC_Current
        self.Vs_Current = Vs_Current
        
    def update(self,NE, NTC, Vs):
        
        self.NE_Current = NE
        self.NTC_Current = NTC
        self.Vs_Current = Vs
        


class Run_program:

    PID=PID_sys(name='PID')
    baldi0d=InCylinder(name='baldi0d')
    Prop = Propeller_sys('Prop')
    EngineRot = Rotationspeed('EngineRot')
    PowerMag = Powermanagement('PowerMag')

    ############################### MODEL PARAMETERS #######################################################
    ############################### MODEL PARAMETERS #######################################################
    ############################### MODEL PARAMETERS #######################################################
    ############################### MODEL PARAMETERS #######################################################
    ############################### MODEL PARAMETERS #######################################################

    (time_t,total_time,deltat2,NTC_ini,prc,prt,Vs_ini,total_time,
    N_cycle,sh_displ,rho_sw,mass_of_ship,D_p,thrust,
    c_1,z_p,A_e_A_o,p_D_p,omega,n1_nump,ct_nump,s1_nump,
    t1_nump,u1_nump,v1_nump,n2_nump,cq_nump,s2_nump,t2_nump,
    u2_nump,v2_nump,m_hydro,zc,kp,kd,ki,xr_o,Nord,
    deltat,rev_cy,Vc,lhv,Sw,stroke,l,ar,R,Tw,p_ivc,
    T_ivc,ma_ini,rc,Tref,Cstoich,re,ra,rf,acomb,bcomb,
    am,bm,cm,delta_m,a_id,b_id,c_id,mWiebe,lamb_ref,
    phi_id_ref,m_ivc_ref,N_ref,mWieberef,delta_phi_comb_ref,SMFR,
    topen,tclose,p_a,T_a,gamma_a,gamma_e,eta_sh,
    eta_comb,eta_t,eta_c,p_ER_ini,T_ER_ini,T_IR_ini,
    p_IR_ini,cpe,cpa,kpac,kpaf,cvir,cver,V_ir,V_er,
    kf0,kf1,kf2,k_ac0,k_ac1,k_ac2,bore,zc,ITC,I_e,
    I_sh,I_p,cd,A_eq,R_a,T_w,Q_ab_max,Q_ab_min,Max_power,
    LHV1,LHV2,LHV3,LHV4,LHV5,MCR1,MCR2,MCR3,MCR4,MCR5,
    an,an2,an3,an4,an5,bn,bn2,bn3,bn4,bn5,
    cn,cn2,cn3,cn4,cn5,PME1_max,PME1_min,nu1ME1,nu2ME1,nuelME1,
    PME2_max,PME2_min,nu1ME2,nu2ME2,nuelME2,
    PAE1_max,PAE1_min,nu1AE1,nu2AE1,nuelAE1,
    PAE2_max,PAE2_min,nu1AE2,nu2AE2,nuelAE2) = (MP.values)

    V_D = zc*math.pi* pow(bore,2)*stroke/4
    k_ep=4.20*pow(10,-12)
    NE_ini =49*2*math.pi/60
    NE_Current = NE_ini
    #initial pressure and temperature in exhaust receiver
    T_exh = T_ER_ini
    p_exh = p_ER_ini
    #initial pressure and temperature in exhaust receiver
    T_IR = T_IR_ini
    p_IR = p_IR_ini





    ##############################Initialisation#######################################






    C_Current= Current(NE_ini, NTC_ini, Vs_ini)


    NE_tot = np.array([NE_ini])
    NTC_tot = np.array([NTC_ini])

    #Putting the different parameters in every single sub models
    Initialisations(  
                    #PID controller model #
                    PID, rev_cy, kp, kd, ki, 
                    xr_o, time_t, T_a, p_a, k_ac0,
                    k_ac1, k_ac2, T_w,  V_D,  cd,
                    A_eq, R_a, gamma_a, eta_c,
                    
                    #Engine 0d Baldi #
                    baldi0d,deltat, zc,
                    delta_phi_comb_ref, lamb_ref, N_ref, m_ivc_ref,phi_id_ref,
                    Vc,
                    mWieberef, lhv,  Tw, p_ivc, T_ivc, ma_ini, mWiebe,
                    a_id,b_id,c_id,stroke,rc,acomb, bcomb, am,bm,cm, 
                    Tref, delta_m, Sw, ra,re,rf, l,ar, Cstoich,R, kf0, kf1, kf2,bore,
                    SMFR, topen, tclose,
                    #Propeller model #
                    Prop,  z_p,A_e_A_o, p_D_p, D_p, omega,  
                    rho_sw, sh_displ, m_hydro, c_1, I_p,pi,cq_nump,ct_nump,s2_nump,t2_nump,u2_nump,v2_nump,s1_nump, 
                    t1_nump,u1_nump,v1_nump,n1_nump,n2_nump,thrust,
                    #Optimisation energy model #
                    PowerMag,  PME1_max,PME2_max, PAE1_max,  PAE2_max,
                    PME1_min,PME2_min, PAE1_min,  PAE2_min,
                    LHV1, LHV2, LHV3, LHV4, LHV5, 
                    MCR1, MCR2, MCR3, MCR4, MCR5, 
                    an,  an2,  an3,  an4,  an5,
                    bn,bn2,bn3,bn4,bn5,
                    cn, cn2, cn3, cn4, cn5, 
                    nu1ME1,  nu1ME2, nu1AE1, nu1AE2, 
                    nu2ME1,nu2ME2,nu2AE1,nu2AE2,
                    nuelME1,   nuelME2,   nuelAE1,   nuelAE2, 
                    Q_ab_max, Q_ab_min,
                    #MVEM model #
                    EngineRot,  k_ep, gamma_e, eta_sh, eta_t, I_e, I_sh,
                    ITC, cpa, cpe, eta_comb, kpac, kpaf, cvir, cver, V_ir, V_er )







    count = 0

    #storage of pressure, temperature, time and mass flow
    pstor = np.array(p_IR_ini)
    Tstor = np.array(T_IR_ini)
    mestor = np.array(90919)
    tstor = np.array([0])




    while time_t <total_time:
    
    ########################################### Sub model PID Controller ###########################################
    

    #Inputs for PiD fuel flow #
        #current time
        PID.var_in.time_t=time_t
        #Vector of engine speeds (present and past values)
        PID.var_in.NE= NE_tot
        #Current engine speed
        PID.var_in.NE_Current= C_Current.NE_Current

        #vector of desired speed
        PID.var_in.Nord = Nord
        #counting
        PID.var_in.count = count
        
        #running PID fuel flow model found in System_01a_IO_Compressor: "compute" function is called
        PID.run_once()



    ########################################### Sub model Propeller ###########################################

    

    #Inputs for Propeller #

        #Current time of simulation
        Prop.var_in.time_t = time_t
        # Current speed of ship 
        Prop.var_in.Vs = C_Current.Vs_Current
        # Current engine speed
        Prop.var_in.NE=C_Current.NE_Current

    #running Propeller model found in System_02a_IO_Propeller: "compute" function is called
        Prop.run_once()


    ########################################### Sub model Engine 0d cycle ###########################################


    #Inputs for Engine 0d cycle #
    
        #Temperature inlet 
        baldi0d.var_in.T_IR=T_IR
        #Pressure inlet 
        baldi0d.var_in.p_IR=p_IR
        #Current engine speed 
        baldi0d.var_in.NE_Current=C_Current.NE_Current
        
        baldi0d.var_in.mf_cy=PID.var_out.mf_cy
        
        baldi0d.var_in.t_inj = PID.var_out.t_inj

    #running 0d model found in System_04a_IO_Engine0d: "compute" function is called    
        baldi0d.run_once()
        
        
        



    ########################################## Sub model Engine MVEM ############################################


    #Inputs for Engine MVEM cycle #

        #Power requested requested Torque 
        EngineRot.var_in.Q_p =  Prop.var_out.Power_Motor/C_Current.NE_Current

        EngineRot.var_in.deltat= deltat
    
    
        
        #Temperature inlet receiver
        EngineRot.var_in.T_IR= T_IR
        #pressure inlet receiver
        EngineRot.var_in.p_IR= p_IR
        #pressure exhaust receiver
        EngineRot.var_in.p_exh= p_exh
        #temperature exhaust receiver
        EngineRot.var_in.T_exh= T_exh
        
        
        #pressure ratio for turbine
        EngineRot.var_in.prt= prt
        #pressure ratio for compressor
        EngineRot.var_in.prc= prc
        # fuel flow 
        EngineRot.var_in.mf_dot = PID.var_out.dmf
        # current engine speed
        EngineRot.var_in.NE_Current= C_Current.NE_Current
        # current turbo charger speed
        EngineRot.var_in.NTC=C_Current.NTC_Current
    



        

        EngineRot.var_in.W_dot= baldi0d.ex.W_dot
        EngineRot.var_in.Qw_dot= baldi0d.ex.Qw_dot

    #running Engine Rotation found in System_05a_IO_EngineMVEM: "compute" function is called     
        EngineRot.run_once()
        
        Prop.run_once()

    ########################################## Sub model Next time step ############################################
        
        #taking the different derivatives computed: output of previous sub models
        #derivative of TC speed output of engine MVEM sub model
        dNTC = EngineRot.var_out.dNTC_dt
        #derivative of engine speed output of engine MVEM sub model
        dNE = EngineRot.var_out.dNE_dt

        # updating current engine and turbo charger speed and the ship speed
        NE_Current = C_Current.NE_Current*60/(2*math.pi) + dNE*deltat2
        NTC_Current = C_Current.NTC_Current + dNTC*deltat2
        
        NE_Current = 2*math.pi/60*NE_Current
    
        C_Current.update(NE_Current, NTC_Current, Vs_ini)
        
        p_IR = EngineRot.var_out.p_IR
        p_exh = EngineRot.var_out.p_exh
        T_IR = EngineRot.var_out.T_IR
        T_exh = EngineRot.var_out.T_exh
        
        prc = EngineRot.var_out.prc
        prt = EngineRot.var_out.prt
        
        # time t
        
        time_t = time_t +deltat2
        tstor = np.append(tstor, time_t)
        pstor = np.append(pstor, p_IR)
        Tstor = np.append(Tstor, T_IR)

        
    

        
        # putting actual engine speed at the end of the vector of speeds
        NE_tot=np.append(NE_tot,NE_Current)
        NTC_tot=np.append(NTC_tot,NTC_Current)
        
        
        count = count+1
        
        print(9.5493 *NE_tot)
        
        print(NTC_tot)

    
  





