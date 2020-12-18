"""
P6
July 2020
A. Takken
"""

import numpy as np
from scipy.optimize import newton

import SUB_muPandkP
import SUB_cp
import SUB_HtoT

def fun1(Dch,DmeanM,Lch,Ach,mdotch,Tpi,TRAC,propellant,channellayout,NISTP,limitsP,MMP):
   
    mu = lambda Tb: SUB_muPandkP.fun1(propellant,Tb)[0]                         #[Pa s], dynamic viscosity
    k = lambda Tb: SUB_muPandkP.fun1(propellant,Tb)[1]                          #[W/m/K], thermal conductivity
    cp = lambda Tb: SUB_cp.fun1(Tb,NISTP,limitsP,MMP)                           #[J/kg/K], heat capacity at constant pressure
    Tpo1 = lambda Tb: SUB_HtoT.fun1(NISTP,limitsP,MMP,P6ch(Tb)/mdotch,Tpi,TRAC)      #[K], output temperature
    
    def Tpo2(Tb):
        return 2*Tb-Tpi                                                         #[K], output temperature
    
    def PrP(Tb):
        return mu(Tb)*cp(Tb)/k(Tb)                                              #[-], channel Prandtl number
    
    def ReD(Tb):
        return 4*mdotch/(np.pi*mu(Tb)*Dch)                                      #[-], channel Reynolds number
    
    def hP(Tb):
        return NuP(Tb)*k(Tb)/Dch                                                #[W/m2/K], convective heat transfer coefficient

    if Tpi == TRAC:
        def P6ch(Tb):
            return 0.0
    else:
        def P6ch(Tb):
            return hP(Tb)*Ach*(Tpo2(Tb)-Tpi)/np.log((TRAC-Tpi)/(TRAC-Tpo2(Tb))) #[W], channel propellant convection

    Rei = ReD(Tpi)                                                              #[-], inflow Reynolds number
    
    if channellayout == 0: #straight channels
        if Rei < 2300: #laminar flow                 
            def NuP(Tb):
                return (3.657+0.0677*(ReD(Tb)*PrP(Tb)*Dch/Lch)**1.33/           #[-], Nusselt number (Stephan)
                        (1+0.1*PrP(Tb)*(ReD(Tb)*Dch/Lch)**0.3))  
        elif Rei < 5.0e6: #turbulent flow
            def fDB(Tb):    
                return (0.790*np.log(ReD(Tb))-1.64)**-2                         #[-], friction factor (Bergman)
            def NuP(Tb):
                return (fDB(Tb)/8*(ReD(Tb)-1000)*PrP(Tb)/(1+12.7*(fDB(Tb)/8)    #[-], Nusselt number (Gnielinski (enhanced))
                           **(1/2)*(PrP(Tb)**(2/3)-1))*(1+(Dch/Lch)**(2/3)))    
     
    elif channellayout == 1: #spiral channels
        if Rei < 1.0e4: #laminar flow           
            def NuP(Tb):
                return 0.913*(ReD(Tb)*(Dch/DmeanM)**0.5)**0.476*PrP(Tb)**0.2    #[-], Nusselt number (Kalb & Seader)
        elif Rei < 1.0e5: #turbulent flow      
            def NuP(Tb):
                return 0.023*ReD(Tb)**0.85*PrP(Tb)**0.4*(Dch/DmeanM)**0.1       #[-], Nusselt number (Seban & McLaughlin)
            
    def RESULTANT(Tb):
        return Tpo1(Tb)-Tpo2(Tb)                                                #resulting newton formula, should approach 0
    
    if Tpi <= TRAC:
        
        Tb = newton(RESULTANT,(Tpi+TRAC)/2-0.10,tol=1.00e-3)                        #[K], resulting bulk temperature
    else:
        Tb = newton(RESULTANT,(Tpi+TRAC)/2+0.10,tol=1.00e-3)                        #[K], resulting bulk temperature
    
    return(Tpo1(Tb),P6ch(Tb),ReD(Tb),PrP(Tb),Tb)
