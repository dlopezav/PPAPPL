# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 17:20:37 2022

@author: hbusson
"""
from __future__ import division
import numpy as np


import matplotlib.pyplot as plt

    
    


from pyomo.environ import (Binary, ConcreteModel, Constraint,
                           NonNegativeReals, Objective,
                           RangeSet, Var, minimize)
from pyomo.common.collections import ComponentMap



def Engines(nbmoteur1, nbmoteur2, mw1, mw2,an,bn,cn,an2,bn2,cn2,load_ini, load_ini2, storage, 
            timeperiod,  nch, ndch,PchMax, PdchMax, Emax,Eini,deltat, MCR, MCR2):
    list1 = np.array([]) 
    i=0
    #building an object from class engine at each iteration and adding to the list list1
    while i < nbmoteur1:
    #current load is chosen at 0.8 by default
        P=mw1*np.zeros((timeperiod))
        load =np.zeros((timeperiod)) 
        load[0]=load_ini
        P[0] = mw1*load_ini
        list1= np.append( list1,Engine( mw1,load_ini,P[0],an,bn,cn, MCR,nbmoteur1) )
        i=i+1
    list2 = np.array([])
    i=0
    #building an object from class engine at each iteration and adding to the list list2
    while i < nbmoteur2:
        P2=mw2*np.zeros((timeperiod))
        load2 =np.zeros((timeperiod)) 
        load2[0]=load_ini2
        P2[0] = mw2*load_ini2
        list2= np.append( list2,Engine( mw2,load_ini2,P2[0],an2,bn2,cn2, MCR2,nbmoteur2) )
        i=i+1

        
        
    stor=np.empty(
            (timeperiod),
            dtype=object
           )
    totlist=np.empty(
            (nbmoteur1+nbmoteur2),
            dtype=object
           )
 
        
    for t in range(timeperiod):
        if storage == 1:
            if t==0:
                stor[t]=Storage(0, nch, ndch,PchMax, PdchMax, Emax, Eini )
            if t>0:
                stor[t]=Storage(stor[t-1].Energystorage(deltat) , nch, ndch,PchMax, PdchMax, Emax,Eini )
        
    for i in range(nbmoteur1):
        totlist[i]=list1[i]
    for i in range(nbmoteur2):
        totlist[nbmoteur1+i]=list2[i]
        

       

    
   
    return totlist,list1,list2, stor

class Engine:
     def __init__(self, Pmax,load,P,an,bn,cn, MCR,n ):
         #max power of the engine
            self.Pmax = Pmax
         #current load of the engine
            self.load=load
         # current power
            self.P =P

        #maximum continuous rating
            self.MRC = MCR
        #number of components for this type of engines
            self.n = n
        #efficiency
            self.nu= cn+bn*self.load+an*pow(self.load,2)
            self.an=an
            self.bn=bn
            self.cn=cn
     def mass_flow(self,LHV):
         
         self.dm_i = self.load * self.MCR * self.n / (self.nu*LHV)
         
     def maxrate(self,Rimax):
        
        self.Rimax = Rimax
        
            
     # function returning the CSP with current load  
     def CSP(self,t,a,b,c):
         #a,b and c are parameters used to define the CSP curve as a function
         # of load: CSP=  a*load^2 + b*load + c 
        self.a = a
        self.b = b
        self.c = c
        return self.c+self.b*self.P[t]+self.a*pow(self.P[t],2)
     def fuelconsumption(self,t):
         return (self.c+self.b*self.load+self.a*pow(self.load,2))*self.P
     #function plotting the current CSP with current load
     def CSPplot(self,t):
        # affiche le graphe
        yp = self.c+self.b*self.load[t]+self.a*pow(self.load[t],2)
        plt.plot(self.load[t],yp, color='r', marker='+')
      #function plotting the CSP curve for all possible load
     def CSPcurve(self):
        xs = np.arange(0, 1,0.1)
        ys = self.c+self.b*xs+self.a*pow(xs,2)
        # dessine
        plt.plot(xs, ys, color='b')

        

class WHR_system:
    def __init__(self, Q_exh, T_exh, m_exh , cp):
        
        self.Q_exh= Q_exh
        self.T_exh = T_exh
        self.m_exh = m_exh
        self.cp= cp
        
    def compute_Qexh(self, Tmin):
        
        self.q_exh = self.m_exh*self.cp*(self.T_exh-Tmin)
        
    

class Storage:
    def __init__(self, EESSm1, nch, ndch, PchMax, PdchMax, Emax, Eini,uess ):
        self.Pch = 0
        self.Pdch =0
        self.PchMax = PchMax
        self.PdchMax = PdchMax
        self.uess = uess
        self.EESSm1 = EESSm1
        self.nch = nch
        self.ndch = ndch
        self.Emax = Emax
        self.Eini = Eini
    def Energystorage(self,deltat):          
        EESS =self.EESSm1 +self.uess*self.Pch*deltat*self.nch -(1-self.uess)*(self.Pdch/self.ndch)*deltat
        return EESS



class MINLP(ConcreteModel):
    """Convex MINLP problem Assignment 6 APSE."""

    def __init__(self,Input_power, PME1_max,PME2_max, PAE1_max,  PAE2_max,
                 PME1_min,PME2_min, PAE1_min,  PAE2_min,
                  LHV1, LHV2, LHV3, LHV4, LHV5, 
                  MCR1, MCR2, MCR3, MCR4, MCR5, 
                  an,  an2,  an3,  an4,  an5,
                  bn,bn2,bn3,bn4,bn5,
                  cn, cn2, cn3, cn4, cn5, 
                  nu1ME1,  nu1ME2, nu1AE1, nu1AE2, 
                  nu2ME1,nu2ME2,nu2AE1,nu2AE2,
                  nuelME1,   nuelME2, nuelAE1, nuelAE2, 
                  Q_ab_max, Q_ab_min,
                 *args, **kwargs):
        """Create the problem."""
        kwargs.setdefault('name', 'MINLP')
        super(MINLP, self).__init__(*args, **kwargs)
        m = self

        """Set declarations"""
        I = m.I = RangeSet(1,13 , doc='continuous variables')
        J = m.J = RangeSet(1, 5, doc='discrete variables')

        # initial point information for discrete variables
        initY = {
            'sub1': {1: 1, 2:1, 3: 1, 4:1, 5:1},
            'sub2': {1: 1, 2:1, 3: 1, 4:1, 5:1},
            'sub3': {1: 1, 2:0, 3: 1, 4:0, 5:0},
        }
        # initial point information for continuous variables
        initX = {1: 5.0, 2: 5.0, 3:5.0, 4:5.0, 5:5.0, 6:5.0, 7:5.0, 8:5.0, 9:5.0, 10:5.0, 11:5.0, 12:5.0,
                 13:5.0}

        """Variable declarations"""
        # DISCRETE VARIABLES
        Y = m.Y = Var(J, domain=Binary, initialize=initY['sub1'])
        # CONTINUOUS VARIABLES
        X = m.X = Var(I, domain=NonNegativeReals, initialize=initX)

        """Constraint definitions"""
        # CONSTRAINTS
        m.const1 = Constraint(expr=-PME1_max*(m.X[1]+ m.X[2] + m.X[3]) + PME1_min*m.Y[1] <= 0)
        m.const2 = Constraint(expr=PME1_max*(m.X[1] + m.X[2] + m.X[3])  - PME1_max*m.Y[1] <= 0)
        m.const3 = Constraint(expr=PME2_max*(m.X[4] + m.X[5] + m.X[6])  - PME2_max*m.Y[2] <= 0)
        m.const4 = Constraint(expr=-PME2_max*(m.X[4]+ m.X[5] + m.X[6]) + PME2_min*m.Y[2] <= 0)
        
        m.const5 = Constraint(expr=-PAE1_max*(m.X[7]  + m.X[8] + m.X[9])    + PAE1_min*m.Y[3] <= 0)
        m.const6 = Constraint(expr= PAE1_max*(m.X[7]  + m.X[8] + m.X[9])    - PAE1_max*m.Y[3] <= 0)
        m.const7 = Constraint(expr= PAE2_max*(m.X[10] + m.X[11] + m.X[12])  - PAE2_max*m.Y[4]  <= 0)
        m.const8 = Constraint(expr=-PAE2_max*(m.X[10] + m.X[11] + m.X[12])  + PAE2_min*m.Y[4]  <= 0)
        
        m.const9 = Constraint(expr=    Input_power[0]*0.99 -  PME1_max*m.X[1]*nu1ME1 -  PME2_max*m.X[4]*nu1ME2     -  PAE1_max*m.X[7]*nu1AE1 -  PAE2_max*m.X[10]*nu1AE2 <= 0)
        m.const12 = Constraint(expr=   Input_power[0]*1.01 -  PME1_max*m.X[1]*nu1ME1 -  PME2_max*m.X[4]*nu1ME2     -  PAE1_max*m.X[7]*nu1AE1 -  PAE2_max*m.X[10]*nu1AE2 >= 0)
        m.const10 = Constraint(expr=   Input_power[1]*0.99 -  PME1_max*m.X[2]*nu2ME1 -  PME2_max*m.X[5]*nu2ME2     -  PAE1_max*m.X[8]*nu2AE1 -  PAE2_max*m.X[11]*nu2AE2 <= 0)
        m.const13 = Constraint(expr=   Input_power[1]*1.01 -  PME1_max*m.X[2]*nu2ME1 -  PME2_max*m.X[5]*nu2ME2     -  PAE1_max*m.X[8]*nu2AE1 -  PAE2_max*m.X[11]*nu2AE2 >= 0)
        m.const11 = Constraint(expr=   Input_power[2]*0.99    -  PME1_max*m.X[3]*nuelME1 - PME2_max*m.X[6]*nuelME2 - PAE1_max*m.X[9]*nuelAE1 - PAE2_max*m.X[12]*nuelAE2 <= 0)
        m.const14 = Constraint(expr=   Input_power[2]*1.01    -  PME1_max*m.X[3]*nuelME1 - PME2_max*m.X[6]*nuelME2 - PAE1_max*m.X[9]*nuelAE1 - PAE2_max*m.X[12]*nuelAE2 >= 0)
        #check that
        m.const15 = Constraint(expr= Input_power[3]*0.99 - Q_ab_max*m.X[13] <= 0)
        m.const16 = Constraint(expr= Input_power[3]*1.01 - Q_ab_max*m.X[13] >= 0)
        m.const17 = Constraint(expr= Q_ab_max*m.Y[5] - Q_ab_max*m.X[13] >= 0)
        m.const18 = Constraint(expr= Q_ab_min*m.Y[5] - Q_ab_max*m.X[13] <= 0)

        """Cost (objective) function definition"""
        m.objective = Objective(expr=(m.X[1]+m.X[2]+m.X[3])*MCR1*Y[1]/((cn + bn*(m.X[1]+m.X[2]+m.X[3])+ an*pow(m.X[1]+m.X[2]+m.X[3],2))*LHV1)
                                +(m.X[4]+m.X[5]+m.X[6])*MCR2*Y[2]/(cn2 + bn2*(m.X[4]+m.X[5]+m.X[6])+ an2*pow(m.X[4]+m.X[5]+m.X[6],2)*LHV2)
                                +(m.X[7]+m.X[8]+m.X[9])*MCR3*Y[3]/(cn3 + bn3*(m.X[7]+m.X[8]+m.X[9])+ an3*pow(m.X[7]+m.X[8]+m.X[9],2)*LHV3)
                                +(m.X[10]+m.X[11]+m.X[12])*MCR4*Y[4]/(cn4 + bn4*(m.X[10]+m.X[11]+m.X[12])+ an4*pow(m.X[10]+m.X[11]+m.X[12],2)*LHV4)
                                +(m.X[13])*MCR5*Y[5]/(cn5 + bn5*(m.X[13])+ an5*pow(m.X[13],2)*LHV5)
                                ,sense=minimize)
        """Bound definitions"""


def func_test(X,Y, PME1_max, PME1_min, PME2_max, PME2_min, PAE1_max, PAE1_min, PAE2_min, PAE2_max ,MCR1,MCR2,MCR3,MCR4,MCR5,LHV1,LHV2,LHV3,LHV4,LHV5,an, an2, an3, an4, an5,bn, bn2, bn3, bn4, bn5,cn, cn2, cn3, cn4, cn5,ME1nu1, ME2nu1, AE1nu1, AE2nu1,ME1nu2, ME2nu2, AE1nu2, AE2nu2,ME1nuel, ME2nuel, AE1nuel, AE2nuel,Pprop1, Pprop2, Pel):

    y=(X[1]+X[2]+X[3])*MCR1*Y[1]/((cn + bn*(X[1]+X[2]+X[3])+ an*pow(X[1]+X[2]+X[3],2))*LHV1)
    +(X[4]+X[5]+X[6])*MCR2*Y[2]/(cn2 + bn2*(X[4]+X[5]+X[6])+ an2*pow(X[4]+X[5]+X[6],2)*LHV2)
    +(X[7]+X[8]+X[9])*MCR3*Y[3]/(cn3 + bn3*(X[7]+X[8]+X[9])+ an3*pow(X[7]+X[8]+X[9],2)*LHV3)
    +(X[10]+X[11]+X[12])*MCR4*Y[4]/(cn4 + bn4*(X[10]+X[11]+X[12])+ an4*pow(X[10]+X[11]+X[12],2)*LHV4)
    +(X[13])*MCR5*Y[5]/(cn5 + bn5*(X[13])+ an5*pow(X[13],2)*LHV5)
    
    return y




   
        
