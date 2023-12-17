# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 08:48:04 2022

@author: hbusson
"""

from cosapp.base import System, Port
from pyomo.core.expr.calculus.diff_with_sympy import differentiate_available

import pyomo.common.unittest as unittest

from System_03b_F_Powermanagement import *


from pyomo.environ import SolverFactory, value, maximize
from pyomo.solvers.tests.models.LP_unbounded import LP_unbounded
from pyomo.solvers.tests.models.QCP_simple import QCP_simple
from pyomo.opt import TerminationCondition

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


class PortPower(Port):
    """Handles the cosapp input Port for the power management model.

    The attributes are initiated with cosapp's setup method.

    Attributes
    ----------
    Input_Power: float
        The power incoming to the management system to be distributed among the different engines and/or boilers
    """
    def setup(self):
        self.add_variable('Input_Power')
        
class PortOptim(Port):
    """Handles the cosapp output Port for the power management model.

    The attributes are initiated with cosapp's setup method.

    Attributes
    ----------
    fin_val: float
        This is the final value of the optimization of power management: ie the minimal value reached by the algorithm
    X: np.ndarray
        Continue variables to be optimised: load for each engine/boiler that is minimizing fuel consumption
    Y: np.ndarray
        Discrete variables to be optimised: 0 for engine turned off and 1 for engine on
    term_cond: str
        Terminal condition of the optimisation
    """
    def setup(self):
        self.add_variable('fin_val')
        self.add_variable('X', dtype= np.ndarray,value=np.array([]))
        self.add_variable('Y', dtype= np.ndarray,value=np.array([]))
        self.add_variable('term_cond', dtype=str, value= 'not run yet')

class Powermanagement(System):

    def setup(self):
        """Defines system structure: max and min power that can be reached, lower heating values, efficiencies
        different coefficients that describes the fuel consumption as a function of load"""
        self.add_input(PortPower, 'power_in')    # define a new input port
        self.add_output(PortOptim, 'power_out')  # define a new output port
        # Solitary variables can also be added,
        # as either `inward` or `outward` variables:
        self.add_inward('PME1_max', 10.0)
        self.add_inward('PME2_max', 10.0)
        self.add_inward('PAE1_max', 15.0)
        self.add_inward('PAE2_max', 15.0)
    
        
        self.add_inward('PME1_min', 0.0)
        self.add_inward('PME2_min', 0.0)
        self.add_inward('PAE1_min', 0.0)
        self.add_inward('PAE2_min', 0.0)
   
        
        self.add_inward('LHV1',100.0)
        self.add_inward('LHV2',100.0)
        self.add_inward('LHV3',100.0)
        self.add_inward('LHV4',100.0)
        self.add_inward('LHV5',100.0)
        
        self.add_inward('MCR1', 10.0)
        self.add_inward('MCR2', 10.0)
        self.add_inward('MCR3', 15.0,)
        self.add_inward('MCR4', 15.0)
        self.add_inward('MCR5', 10.0)
        
        self.add_inward('an', 1.0)
        self.add_inward('an2', 1.0)
        self.add_inward('an3', 1.0)
        self.add_inward('an4', 1.0)
        self.add_inward('an5', 1.0)
        
        self.add_inward('bn', 1.0)
        self.add_inward('bn2', 1.0)
        self.add_inward('bn3', 1.0)
        self.add_inward('bn4', 1.0)
        self.add_inward('bn5', 1.0)
        
        self.add_inward('cn', 1.0)
        self.add_inward('cn2', 1.0)
        self.add_inward('cn3', 1.0)
        self.add_inward('cn4',1.0)
        self.add_inward('cn5', 1.0)
        
        self.add_inward('nu1ME1' ,0.8)
        self.add_inward('nu1ME2' ,0.8)
        self.add_inward('nu1AE1' ,0.8)
        self.add_inward('nu1AE2' ,0.8)
      
        
        self.add_inward('nu2ME1', 0.8)
        self.add_inward('nu2ME2', 0.8)
        self.add_inward('nu2AE1', 0.8)
        self.add_inward('nu2AE2', 0.8)
     
      
        
        self.add_inward('nuelME1', 0.8)
        self.add_inward('nuelME2', 0.8)
        self.add_inward('nuelAE1', 0.8)
        self.add_inward('nuelAE2', 0.8)
     
        
        self.add_inward('Q_ab_max', 30.0)
        self.add_inward('Q_ab_min', 3)
        self.add_outward('hop', 0.0)
        

    
    def compute(self): # `compute` method defines what the system does
    
        """ performs the optimisation """
        model = MINLP(self.power_in.Input_Power, self.PME1_max,self.PME2_max, self.PAE1_max,  self.PAE2_max,
                      self.PME1_min,self.PME2_min, self.PAE1_min,  self.PAE2_min,
                      self.LHV1, self.LHV2, self.LHV3, self.LHV4, self.LHV5, 
                      self.MCR1, self.MCR2, self.MCR3, self.MCR4, self.MCR5, 
                      self.an,  self.an2,  self.an3,  self.an4,  self.an5,
                      self.bn,self.bn2,self.bn3,self.bn4,self.bn5,
                      self.cn, self.cn2, self.cn3, self.cn4, self.cn5, 
                      self.nu1ME1,  self.nu1ME2, self.nu1AE1, self.nu1AE2, 
                      self.nu2ME1,self.nu2ME2,self.nu2AE1,self.nu2AE2,
                      self.nuelME1,   self.nuelME2,   self.nuelAE1,   self.nuelAE2, 
                      self.Q_ab_max, self.Q_ab_min)


 

        required_solvers = ('ipopt', 'glpk')
        if all(SolverFactory(s).available() for s in required_solvers):
            subsolvers_available = True
        else:
            subsolvers_available = False


        #@unittest.skipIf(not subsolvers_available,
        #                 'Required subsolvers %s are not available'
        #                 % (required_solvers,))
        #@unittest.skipIf(not differentiate_available,
        #                 'Symbolic differentiation is not available')

        results = SolverFactory('mindtpy').solve(model, strategy='OA',
                                            init_strategy='rNLP',
                                            mip_solver=required_solvers[1],
                                            nlp_solver=required_solvers[0],
                                            iteration_limit=50000,
                                            heuristic_nonconvex=True
                                            )

        self.power_out.fin_val = value(model.objective.expr)
        X= np.array([])
        Y= np.array([])
      
        for i in range(len(model.X) ):
            X = np.append(X,model.X[i+1].value)
        for i in range(len(model.Y) ):  
            Y = np.append(Y,model.Y[i+1].value)
      
        self.power_out.X = X  
        self.power_out.Y = Y
        self.power_out.term_cond =results.solver.termination_condition
        
    def checkresults(self)   : 
        print('CHECK POWER MANAGEMENT')
        
     
        
        print(' loads ', self.power_out.X)
        print('engine on/off dummies', self.power_out.Y)
        print('terminal conditions',   self.power_out.term_cond)
        
        print('POWER to propeller',  self.power_out.X[0] *self.PME1_max*self.nu1ME1+
              self.power_out.X[3] *self.PME2_max*self.nu1ME2+
              self.power_out.X[6] *self.PAE1_max*self.nu1AE1+
              self.power_out.X[9] *self.PAE2_max*self.nu1AE2)
        
        print('POWER input', self.power_in.Input_Power)
        
      
        

        

        
        

        
        
        
        

