"""
P1, P2 & P3
July 2020
A. Takken
"""

import numpy as np
from scipy.optimize import newton
import math

def fun1(DouterA,DouterM,RACtype,ARACo,emO,TRAC,Tamb,pamb,g0,R_A,sigma,kI,insulation,LcavA):
    MMair = 28.9647e-3                                                                      #[kg/mol], molar mass air
    SuthC = [18.27e-6,120,291.15]                                                           #[-], Sutherland constants
        
    mu =    lambda Tf: SuthC[0]*(SuthC[2]+SuthC[1])/(Tf+SuthC[1])*(Tf/SuthC[2])**(3/2)      #[Pa s], dynamic viscosity
    cp =    lambda Tf: -8.0144e-08*Tf**3+2.1079e-04*Tf**2+2.0633e-02*Tf+9.8367e+02          #[J/kg/K], heat capacity at constant pressure
    kair =  lambda Tf: 6.25216e-05*Tf + 7.51105e-03                                         #[W/m/K], thermal conductivity air, from NIST & Excel
    rho =   lambda Tf: pamb/(R_A/MMair)/Tf                                                  #[kg/m3], density
    beta =  lambda Tf: 1/Tf                                                                 #[1/K], expansion coefficient
    Pr =    lambda Tf: mu(Tf)*cp(Tf)/kair(Tf)                                               #[-], Prandtl number
    nu =    lambda Tf: mu(Tf)/rho(Tf)                                                       #[m2/s], kinematic viscosity
    alpha = lambda Tf: nu(Tf)/Pr(Tf)                                                        #[m2/s], diffusivity
    Ra =    lambda Tf: g0*beta(Tf)*(Tinsu(Tf)-Tamb)*DouterA**3/nu(Tf)/alpha(Tf)             #[-], Rayleigh number

    if      RACtype == 0: #cone
                Nu =    lambda Tf: 0.7+0.35*abs(Ra(Tf))**0.125+0.51*abs(Ra(Tf))**0.25           #[-], cone Nusselt number
    elif    RACtype == 1: #cylinder
            if TRAC == Tamb:
                Nu =    lambda Tf: 0                                                            #[-], cylinder Nusselt number  
            else:
                Nu =    lambda Tf: (2/np.log(1+2/((0.518*abs(Ra(Tf))**(1/4)*(1+(0.559/Pr(Tf))   #[-], cylinder Nusselt number
                        **(3/5))**(-5/12))**15+(0.1*abs(Ra(Tf))**(1/3))**15)**(1/15)))                  

    h =     lambda Tf: Nu(Tf)*kair(Tf)/DouterA                              #[W/m2/K], convective heat transfer coefficient
    P1 =    lambda Tf: h(Tf)*ARACo*(Tinsu(Tf)-Tamb)                         #[W], power lost to outer wall convection
    P2 =    lambda Tf: emO*sigma*ARACo*(Tinsu(Tf)**4-Tamb**4)               #[W], power lost to outer wall radiation
                                  
    if insulation == 0: #no insulation
        Tf =    (TRAC+Tamb)/2               #[K], film temperature
        Tinsu = lambda Tf: TRAC             #[K], insulation temperature (equal to RAC temperature)
        P3 =    lambda Tf: 0.0              #[W], no conduction due to lack of insulation
    else: #insulation
        if insulation == 1: #Saffil M-Fil  
            kII =    lambda Tf: 0.0665*math.exp(0.0015*((TRAC+Tinsu(Tf))/2-273.15))         #[W/m/K], thermal conductivity
        else:
            kII =    lambda Tf: kI                                                          #[W/m/K], thermal conductivity
        Tinsu = lambda Tf: 2*Tf-Tamb                                                        #[K], insulation temperature
        P3 = lambda Tf: kII(Tf)*(2*np.pi*LcavA)/np.log((DouterA)/DouterM)*(TRAC-Tinsu(Tf))  #[W], conduction through insulation
        
        RESULTANT = lambda Tf: P3(Tf)-P1(Tf)-P2(Tf)                             #[-], resultant equation
        Tf = newton(RESULTANT,TRAC)                                             #[-], resulting film temperature

    return(P1(Tf),P2(Tf),P3(Tf),Tinsu(Tf),h(Tf))
