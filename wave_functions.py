# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 10:48:25 2023

@author: hbusson
"""

import scipy.integrate as integrate

import scipy.special as special

import math

import numpy as np

class water_resistance:
    
      
    def __init__(self,Lwl, Dp, B, displ, rho_sw, Cb, Lpp, Ta, Tf, T, kyy, LR, LE):
        
        
   
    #### ship characteristics ######

        self.Lwl = Lwl   #lenght at waterline
        self.Dp = Dp #propeller diameter
        self.B = B #breadth
        self.displ = displ #displacement
        self.rho_sw = rho_sw
        self.Lpp = Lpp #lenght between perpendicular
        self.Ta = Ta #draught at AP
        self.Tf = Tf #draught at AP
        self.T = T #draught at midship
        

        self.Cb = Cb #block coefficient
        self.kyy = kyy
        #relative froude number
        self.E1 = math.atan(0.495*B/LE)
        self.E2 = math.atan(0.495*B/LR)
        
    def Stawave1(self,Hs):
        
        g=9.80665
        Rstawave1=1/16*self.rho_sw*g*pow(Hs,2)*self.B*math.sqrt(self.B/self.Lwl)
        
        return Rstawave1
  
    def Stawave2(self,Vs,omega):
        g=9.80665
        f = omega/(2*math.pi)
        period = 1/f
        self.Fr = Vs/math.sqrt(9.80665 *self.Lpp) #froude number
        omega_barre = math.sqrt(self.Lpp/g)* pow(self.kyy,1/3)*omega/(1.17*pow(self.Fr,-0.143))
        a1=  60.3*pow(self.Cb,1.34)
        if omega_barre<1:
            b1=11
            d1=14
        else: 
            b1=-8.5
            d1=-566*pow(self.Lpp/self.B,-2.66)
        
        raw = pow(omega_barre,b1)*math.exp(b1/d1*(1-pow(omega_barre,d1)))*a1*pow(self.Fr,1.5)*math.exp(-3.5*self.Fr)
        
        Rstawave2a = 4 *self.rho_sw* g* pow(self.B,2)/(self.Lpp) * raw
        
        
        
        lamb = g*pow(period,2)/(2*math.pi)

        k = 2*math.pi/lamb
        
        I1 = special.iv(1,1.5*k*self.T)
        K1 = special.kn(1,1.5*k*self.T)
        
        f1=0.692*pow(Vs/(math.sqrt(self.T*g)),0.769)+ 1.81*pow(self.Cb,6.95)
        alpha1 = pow(math.pi,2)*pow(I1,2)/( pow(math.pi,2)*pow(I1,2) +pow(math.pi,2)*pow(K1,2)*f1)
        
        Rstawave2b= 0.5*self.rho_sw*g*self.B*alpha1
     
        Rstawave2 = Rstawave2a+Rstawave2b
        return Rstawave2


    def SNNM(self, alpha, Vs,f):
        
        
        
        E2 = self.E2
        E1 = self.E1
        
        omega = 2*math.pi*f
        period = 1/f
        #amplitude
        Hw =1 
        kyy= 0.25*self.Lpp
        
        #to have kyy, omega, c , E1, E2
        
        g=9.80665
        Vg=g/(2*omega)
        #The wavelength (L) is calculated using: L = gT2/2π, where g=9.8 m/s2
        lamb = g*pow(period,2)/(2*math.pi)
        self.Fr = Vs/math.sqrt(9.80665 *self.Lpp) #froude number
        self.Frrel = (Vs-Vg)/math.sqrt(g*self.Lpp)
        
     
        Tdeep = max(self.Ta, self.Tf)
    
        
        omega_barre=2.142*pow(self.kyy,1/3)*math.sqrt(self.Lpp/lamb)*pow((self.Cb/0.65),0.17)
        

     
        
        
        omega_barre = omega_barre*(1-0.111/self.Cb*(math.log(self.B/Tdeep)-math.log(2.75)))*     \
        ((-1.377*pow(self.Fr,2)+1.157*self.Fr)*abs(math.cos(alpha))+0.618*(13+math.cos(2*alpha))/14)
        
      
    
        
        
        if alpha <= math.pi/2 and alpha >=0:
            a1= pow( (0.87/self.Cb), (1+self.Fr)*math.cos(alpha))/(math.log(self.B/Tdeep))*(1+2*math.cos(alpha))/3
            if self.Fr <0.12:
                a2=0.0072+0.1676*self.Fr
            else:
                a2=pow(self.Fr,1.5)*math.exp(-3.5*self.Fr)
      
            
          
            
            
        if alpha == math.pi  :  
          if Vs >Vg and self.Frrel >=0.12  :
              a1= pow( (0.87/self.Cb), (1+self.Frrel))/(math.log(self.B/Tdeep))
          else : 
              a1=(0.87/self.Cb)/(math.log(self.B/Tdeep))
           
              
          if Vs <= Vg/2: 
              a2=0.0072*(2*Vs/Vg -1)
          elif Vs>Vg/2 and self.Frrel < 0.12:
              a2=0.0072+0.1676*self.Frrel 
          else :
              a2=pow(self.Frrel,1.5)*math.exp(-3.5*self.Frrel)
              
              
        if alpha > math.pi/2 and alpha <math.pi:
            
            if Vs >Vg and self.Frrel >=0.12  :
            
                a1 = alpha*(2-2/math.pi)* pow( (0.87/self.Cb), (1+self.Fr)*math.cos(alpha))/(math.log(self.B/Tdeep))*(1+2*math.cos(alpha))/3 + (alpha*2/math.pi-1)*pow( (0.87/self.Cb), (1+self.Frrel))/(math.log(self.B/Tdeep))
            else :
                a1 = alpha*(2-2/math.pi)* pow( (0.87/self.Cb), (1+self.Fr)*math.cos(alpha))/(math.log(self.B/Tdeep))*(1+2*math.cos(alpha))/3 +  (alpha*2/math.pi-1)*(0.87/self.Cb)/(math.log(self.B/Tdeep))
            
            if self.Fr <0.12:
                
                a2 = alpha*(2-2/math.pi)*0.0072+0.1676*self.Fr
            else:
                
                a2 = alpha*(2-2/math.pi)*pow(self.Fr,1.5)*math.exp(-3.5*self.Fr)
                
            if  Vs <= Vg/2: 
                a2=a2+(alpha*2/math.pi-1)*0.0072*(2*Vs/Vg -1)
                
            elif Vs>Vg/2 and self.Frrel < 0.12:
                a2=a2+(alpha*2/math.pi-1)*(0.0072+0.1676*self.Frrel )
                
            else :
                a2=a2+(alpha*2/math.pi-1)*pow(self.Frrel,1.5)*math.exp(-3.5*self.Frrel)
            
            
            
        a3=1.0+28.7*math.atan( abs(self.Ta-self.Tf)/self.Lpp)    

         
        if omega_barre <1:
            b1=11
            d1=566*pow(self.Lpp*self.Cb/self.B,-2.66)
        else:
            b1=-8.5
            d1=-566*pow(self.Lpp*self.Cb/self.B,-2.66)*(4-125*math.atan(abs(self.Ta-self.Tf)/self.Lpp))
            
            
        
        Rawm = 3859.2*self.rho_sw*g*pow(self.B,2)/self.Lpp*pow(self.Cb,1.34)*pow(self.kyy,2)*a1*a2*a3*pow(omega_barre,b1)
        
    
        Rawm = Rawm*math.exp(b1/d1*(1-pow(omega_barre,d1)) )
        
        Rawr1=0
        Rawr2=0
        Rawr3=0
        Rawr4=0
        
        if alpha >=0 and alpha <= self.E1:
        
            falpha = math.cos(alpha)
        elif alpha > self.E1:
            falpha=0
            
        
        
        
        if alpha >=0 and alpha <= math.pi - self.E1:
           T= Tdeep
           if lamb/self.Lpp <=2.5:
               alpha_T=1-math.exp(-4*math.pi*(T/lamb-T/(2.5*self.Lpp)))    
           else:
               alpha_T = 0
           Rawr1= 2.25/4*self.rho_sw*g*self.B*pow(Hw,2)*alpha_T*\
              (pow( math.sin(E1+alpha),2)+ 2*omega*Vs/g*(math.cos(alpha)-math.cos(E1)*math.cos(E1+alpha)) \
                                                 *pow((0.87/self.Cb),(1+4*math.sqrt(self.Fr))*falpha) )
    
        if alpha >=0 and alpha <= self.E1:   
            T= Tdeep
            if lamb/self.Lpp <=2.5:
                alpha_T=1-math.exp(-4*math.pi*(T/lamb-T/(2.5*self.Lpp)))    
            else:
                alpha_T = 0
            Rawr2= 2.25/4*self.rho_sw*g*self.B*pow(Hw,2)*alpha_T*\
                (pow( math.sin(E1-alpha),2)+ 2*omega*Vs/g*(math.cos(alpha)-math.cos(E1)*math.cos(E1-alpha))\
                 *pow((0.87/self.Cb),(1+4*math.sqrt(self.Fr))*falpha) )
        if alpha >=self.E2 and alpha <= math.pi: 
            if self.Cb <=0.75:
                T=Tdeep*(4+math.sqrt(abs(math.cos(alpha))))/5
            else:
                T=Tdeep*(2+math.sqrt(abs(math.cos(alpha))))/3
            if lamb/self.Lpp <=2.5:
                alpha_T=1-math.exp(-4*math.pi*(T/lamb-T/(2.5*self.Lpp)))    
            else:
                alpha_T = 0    
            Rawr3 =- 2.25/4*self.rho_sw*g*self.B*pow(Hw,2)*alpha_T*\
                (pow( math.sin(E2-alpha),2)+ 2*omega*Vs/g*(math.cos(alpha)-math.cos(E2)*math.cos(E2-alpha)))
               
        if alpha >=self.E2 and alpha <= math.pi:    
            if self.Cb <=0.75:
                T=Tdeep*(4+math.sqrt(abs(math.cos(alpha))))/5
            else:
                T=Tdeep*(2+math.sqrt(abs(math.cos(alpha))))/3
            
            if lamb/self.Lpp <=2.5:
                alpha_T=1-math.exp(-4*math.pi*(T/lamb-T/(2.5*self.Lpp)))    
            else:
                alpha_T = 0
            Rawr4=- 2.25/4*self.rho_sw*g*self.B*pow(Hw,2)*alpha_T*\
                (pow(math.sin(E2+alpha),2)+ 2*omega*Vs/g*(math.cos(alpha)-math.cos(E2)*math.cos(E2+alpha)))
            
        
        Rawr = Rawr1+Rawr2+Rawr3+Rawr4
        return Rawr
        
        
        
    def spectre(self,f,Hp,Tp):

        S_f=5/16*pow(Hp,2) *pow(Tp,4)*pow( f,-5)* math.exp(-5/4*pow((f/Tp),-4))
        
   
        
        
        
        return S_f

    def integrand(self,A,Hp,Tp, Vs,f,alpha,method_str):
        
        
        A= 1
        if method_str == 'SNNM':
            I=self.spectre(f,Hp,Tp)*self.SNNM(alpha, Vs,f)/pow(A,2)
        if method_str == 'Stawave2':
            omega = 2*math.pi*f
            I=self.spectre(f,Hp,Tp)*self.Stawave2(Vs,omega)/pow(A,2)
        
        return I

    def irregular_wave(self, A,Hs,Tp,Vs, alpha, method_str):
        
       
        f = 1/Tp
     
       
        
        
        Riw= 2*integrate.quad(lambda x: self.integrand(A,Hs,Tp,Vs,x, alpha,method_str), 0, np.inf)
        
     
        
     
        
        return Riw[0]
    