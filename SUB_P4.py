"""
P1, P2 & P3
July 2020
A. Takken
"""


def fun1(LsI,ARACi,TRAC,Tamb,pamb,g0,R_A):
  
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
    Ra =    lambda Tf: g0*beta(Tf)*(TRAC-Tamb)*LsI**3/nu(Tf)/alpha(Tf)                      #[-], Rayleigh number
    Nu =    lambda Tf: 0.00324*abs(Ra(Tf))**0.447                                           #[-], Nusselt number

    h =     lambda Tf: Nu(Tf)*kair(Tf)/LsI          #[W/m2/K], convective heat transfer coefficient
    P4 =    lambda Tf: h(Tf)*ARACi*(TRAC-Tamb)      #[W], power lost to outer wall convection

    Tf =    (TRAC+Tamb)/2                           #[K], film temperature
    
    return(P4(Tf),h(Tf))
