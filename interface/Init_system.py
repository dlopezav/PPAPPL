# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 11:19:48 2022

@author: hbusson
"""

import numpy as np


#Memo
#1 check if VD = V_D



def Initialisations(PID,  rev_cy, kp, kd, ki,
                  
                  tinj_o, time_t,  T_a, p_a, k_ac0,
                  k_ac1, k_ac2, T_w,   V_D,  cd,
                  A_eq, R_a, gamma_a, eta_c,
                  baldi0d,deltat, zc,
                  delta_phi_comb_ref, lamb_ref, N_ref, m_ivc_ref,phi_id_ref,
                   Vc,
                  mWieberef, lhv,   Tw, p_ivc, T_ivc, ma_ini, mWiebe,
                  a_id,b_id,c_id,stroke,rc,acomb, bcomb, am,bm,cm, 
                  Tref, delta_m, Sw, ra,re,rf, l,ar, Cstoich,R, kf0, kf1, kf2,bore,
                  SMFR, topen, tclose,
                  Prop,  z_p,A_e_A_o, p_D_p, D_p, omega,  
                                   rho_sw, sh_displ, m_hydro, c_1, I_p,pi,cq,ct,
                                   s2_nump,t2_nump,u2_nump,v2_nump,s1_nump,
                                   t1_nump,u1_nump,v1_nump,n1,n2,thrust,
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
                                   EngineRot,  k_ep, gamma_e, eta_sh, eta_t, I_e, I_sh,
                                                      
                                                            ITC, cpa, cpe, eta_comb ,
                                                            kpac, kpaf, cvir, cver, V_ir, V_er
                                                           
                                                           ):
    
    


########################### SYSTEM COMPRESSOR AIR COOLER #################################
        



    PID.rev_cy=rev_cy
    PID.zc=zc
    
    PID.kp=kp
    PID.kd=kd
    PID.ki=ki
     

  
    
    PID.tinj_o=tinj_o




########################### SYSTEM ENGINE 0d #################################



    baldi0d.deltat=deltat
    baldi0d.N_ref=N_ref
    baldi0d.delta_phi_comb_ref= delta_phi_comb_ref
    baldi0d.m_ivc_ref=m_ivc_ref

    baldi0d.lamb_ref= lamb_ref

    baldi0d.phi_id_ref= phi_id_ref
    baldi0d.Vc=Vc
    baldi0d.mWieberef=mWieberef
    baldi0d.lhv=lhv
   
    
    baldi0d.Tw= Tw
    #tobe considered 
    baldi0d.ma_ini= ma_ini

    baldi0d.a_id= a_id
    baldi0d.b_id= b_id
    baldi0d.c_id= c_id
    baldi0d.stroke= stroke
    
    baldi0d.rc= rc
    baldi0d.acomb= acomb
    baldi0d.bcomb= bcomb
    baldi0d.am= am
 
    baldi0d.bm= bm
    baldi0d.cm= cm
    baldi0d.Tref= Tref
    baldi0d.delta_m= delta_m
    
    baldi0d.Sw= Sw
    baldi0d.re= re
    baldi0d.ra= ra
    baldi0d.rf= rf
 
    baldi0d.l= l
    baldi0d.ar= ar
    baldi0d.Cstoich= Cstoich
    baldi0d.R= R

      
  

    baldi0d.bore= bore
    
    baldi0d.kf0=kf0
    baldi0d.kf1=kf1
    baldi0d.kf2=kf2
   
    baldi0d.rev_cy=rev_cy
    
    baldi0d.SMFR=SMFR
    baldi0d.topen=topen
    baldi0d.tclose=tclose

    
    


    Prop.z_p=z_p
    Prop.A_e_A_o=A_e_A_o
    Prop.p_D_p=p_D_p

    Prop.D_p=D_p
    Prop.omega=omega
    Prop.rho_sw=rho_sw
    Prop.sh_displ= sh_displ
    Prop.m_hydro=m_hydro
    Prop.c_1=c_1

 
    
    Prop.cq=cq
    Prop.ct=ct
    Prop.s1=s1_nump
    Prop.t1=t1_nump
    Prop.u1=u1_nump
    Prop.v1=v1_nump
    Prop.s2=s2_nump
    Prop.t2=t2_nump
    Prop.u2=u2_nump
    Prop.v2=v2_nump
    Prop.n1=n1
    Prop.n2=n2
    Prop.thrust=thrust
   

    




    EngineRot.p_a=p_a
    EngineRot.k_ep=k_ep

    EngineRot.gamma_e = gamma_e

    EngineRot.eta_t = eta_t
    EngineRot.cpa = cpa
    EngineRot.cpe = cpe
    EngineRot.cpc = cpa
    EngineRot.cpt = cpe
 

    EngineRot.rev_cy= rev_cy

    EngineRot.ITC = ITC
    EngineRot.I_e = I_e
    EngineRot.I_sh = I_sh
    EngineRot.I_p = I_p
    EngineRot.eta_sh = eta_sh
    EngineRot.T_a= T_a
 
    EngineRot.k_ac0=k_ac0
    EngineRot.k_ac1=k_ac1
    EngineRot.k_ac2=k_ac2
    EngineRot.T_w=T_w

    EngineRot.R_a=R_a
   
     
    EngineRot.rc=rc

    EngineRot.V_D=V_D
 
    EngineRot.cd=cd
    EngineRot.A_eq=A_eq


        
    EngineRot.gamma_a=gamma_a
    EngineRot.eta_c=eta_c
    EngineRot.R_a=R_a
    
    EngineRot.eta_comb=eta_comb
    EngineRot.lhv=lhv
    
            
    EngineRot.kpac=kpac
    EngineRot.kpaf=kpaf
    EngineRot.cvir=cvir
    
    EngineRot.cver=cver
    EngineRot.V_ir=V_ir
    EngineRot.V_er=V_er
    
    EngineRot.kf0=kf0
    EngineRot.kf1=kf1
    EngineRot.kf2=kf2
    

 



 



    PowerMag.PME1_max=PME1_max
    PowerMag.PME2_max=PME2_max
    PowerMag.PAE1_max=PAE1_max
    PowerMag.PAE2_max=PAE2_max


    PowerMag.PME1_min=PME1_min
    PowerMag.PME2_min=PME2_min
    PowerMag.PAE1_min=PAE1_min
    PowerMag.PAE2_min=PAE2_min
    

    PowerMag.LHV1=LHV1
    PowerMag.LHV2=LHV2
    PowerMag.LHV3=LHV3
    PowerMag.LHV4=LHV4
    PowerMag.LHV5=LHV5

    PowerMag.MCR1=MCR1
    PowerMag.MCR2=MCR2
    PowerMag.MCR3=MCR3
    PowerMag.MCR4=MCR4
    PowerMag.MCR5=MCR5
            
    PowerMag.an=an
    PowerMag.an2=an2
    PowerMag.an3=an3
    PowerMag.an4=an4
    PowerMag.an5=an5
    
    PowerMag.bn=bn
    PowerMag.bn2=bn2
    PowerMag.bn3=bn3
    PowerMag.bn4=bn4
    PowerMag.bn5=bn5

    PowerMag.cn=cn
    PowerMag.cn2=cn2
    PowerMag.cn3=cn3
    PowerMag.cn4=cn4
    PowerMag.cn5=cn5

    PowerMag.nu1ME1=nu1ME1
    PowerMag.nu1ME2=nu1ME2
    PowerMag.nu1AE1=nu1AE1
    PowerMag.nu1AE2=nu1AE2
         
       
    PowerMag.nu2ME1=nu2ME1
    PowerMag.nu2ME2=nu2ME2
    PowerMag.nu2AE1=nu2AE1
    PowerMag.nu2AE2=nu2AE2
     
      
        
    PowerMag.nuelME1=nuelME1
    PowerMag.nuelME2=nuelME2
    PowerMag.nuelAE1=nuelAE1
    PowerMag.nuelAE2=nuelAE2
     
        
    PowerMag.Q_ab_max=Q_ab_max
    PowerMag.Q_ab_min=Q_ab_min