"""
P6
July 2020
A. Takken
"""

import numpy as np
from scipy.optimize import minimize_scalar

import SUB_muPandkP
import SUB_cp
import SUB_HtoT
import SUB_TtoH

def fun1(Dh,DmeanM,Lch,Aheat,mdot,mdotch,Tpi,TRAC,propellant,channellayout,NISTP,limitsP,MMP,Acs):

    mu = lambda Tb: SUB_muPandkP.fun1(propellant,Tb)[0]                             #[Pa s], dynamic viscosity
    k = lambda Tb: SUB_muPandkP.fun1(propellant,Tb)[1]                              #[W/m/K], thermal conductivity
    cp = lambda Tb: SUB_cp.fun1(Tb,NISTP,limitsP,MMP)                               #[J/kg/K], heat capacity at constant pressure
    Tpo1 = lambda Tb: SUB_HtoT.fun1(NISTP,limitsP,MMP,P6(Tb)/mdot,Tpi,TRAC)         #[K], output temperature
  
    def Tpo2(Tb):
        return 2*Tb-Tpi                                                         #[K], output temperature
    
    def PrP(Tb):
        return mu(Tb)*cp(Tb)/k(Tb)                                              #[-], channel Prandtl number
    
    def ReD(Tb):
        return mdotch*Dh/(Acs*mu(Tb))                                           #[-], channel Reynolds number
    
    def hP(Tb):
        return NuP(Tb)*k(Tb)/Dh                                                 #[W/m2/K], convective heat transfer coefficient

    if Tpi == TRAC:
        def P6(Tb):
            return 0.0
    else:
        def P6(Tb):
            return hP(Tb)*Aheat*(Tpo2(Tb)-Tpi)/np.log((TRAC-Tpi)/(TRAC-Tpo2(Tb))) #[W], channel propellant convection

    Rei = ReD(Tpi)                                                              #[-], inflow Reynolds number
    
    if channellayout == 0: #straight channels
        if Rei < 2300: #laminar flow                 
            def NuP(Tb):
                return (3.657+0.0677*(ReD(Tb)*PrP(Tb)*Dh/Lch)**1.33/            #[-], Nusselt number (Stephan)
                        (1+0.1*PrP(Tb)*(ReD(Tb)*Dh/Lch)**0.3))  
        elif Rei < 5.0e6: #turbulent flow
            def fDB(Tb):    
                return (0.790*np.log(ReD(Tb))-1.64)**-2                         #[-], friction factor (Bergman)
            def NuP(Tb):
                return (fDB(Tb)/8*(ReD(Tb)-1000)*PrP(Tb)/(1+12.7*(fDB(Tb)/8)    #[-], Nusselt number (Gnielinski (enhanced))
                           **(1/2)*(PrP(Tb)**(2/3)-1))*(1+(Dh/Lch)**(2/3)))    
     
    elif channellayout == 1: #spiral channels
        if Rei < 1.0e4: #laminar flow           
            def NuP(Tb):
                return 0.913*(ReD(Tb)*(Dh/DmeanM)**0.5)**0.476*PrP(Tb)**0.2     #[-], Nusselt number (Kalb & Seader)
        elif Rei < 1.0e5: #turbulent flow      
            def NuP(Tb):
                return 0.023*ReD(Tb)**0.85*PrP(Tb)**0.4*(Dh/DmeanM)**0.1        #[-], Nusselt number (Seban & McLaughlin)
            
    def RESULTANT(Tb):
        return abs(Tpo1(Tb)-Tpo2(Tb))                                                #resulting formula, should approach 0
  
    if Tpi <= TRAC: 
        try:
            Tb = minimize_scalar(RESULTANT,bounds=[Tpi,(Tpi+TRAC)/2], method='bounded',options={'xatol': 1e-5,'maxiter':100})     #[K], resulting bulk temperature    
            Tb = Tb["x"]
        except:
            Tb = (TRAC+Tpi)/2-0.001
            HTpo = SUB_TtoH.fun1(NISTP,limitsP,MMP,TRAC-0.001,Tpi)
            P66 = HTpo*mdot      
    else:
        Tb = minimize_scalar(RESULTANT,bounds=[(Tpi+TRAC)/2,Tpi], method='bounded',options={'xatol': 1e-5,'maxiter':100})     #[K], resulting bulk temperature   
        Tb = Tb["x"]
        
    P66 = P6(Tb)
    if abs(Tpo1(Tb)-Tpo2(Tb)) > 1.0:
        Tb = (TRAC+Tpi)/2-0.001
        HTpo = SUB_TtoH.fun1(NISTP,limitsP,MMP,TRAC-0.001,Tpi)
        P66 = HTpo*mdot                    

    return(Tpo2(Tb),P66,ReD(Tb),PrP(Tb),Tb)
