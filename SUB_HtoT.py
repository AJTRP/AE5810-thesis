"""
Calculate propellant temperature from enthalpy
July 2020
A. Takken
"""

from scipy.optimize import newton

def fun1(NISTP,limitsP,MMP,H6,Tpi,TRAC):

    def f(i,T):
        return ((NISTP[0][i]*T/1000+NISTP[1][i]*(T/1000)**2/2                               #[J/kg]
                +NISTP[2][i]*(T/1000)**3/3+NISTP[3][i]*(T/1000)**4/4
                -NISTP[4][i]/(T/1000)+NISTP[5][i]-NISTP[6][i])*1000/MMP)
    
    Hlimit0 = f(0,limitsP[0]) #[J/kg], limit 0
    Hlimit1 = f(1,limitsP[1]) #[J/kg], limit 1   

    ### Starting enthalpy of incoming propellant flow (0 at 298.15 K)
    if Tpi < limitsP[0]:
        Tpii = 0
    elif Tpi < limitsP[1]:
        Tpii = 1
    else:
        Tpii = 2
        
    HTpi = f(Tpii,Tpi)
    deltaH6 = H6+HTpi
    
    if deltaH6 < Hlimit0:
        ii = 0
    elif deltaH6 < Hlimit1:
        ii = 1
    else:
        ii = 2
        
    def f2(Tpo):
        return (deltaH6-f(ii,Tpo))
  
    Tpo = newton(f2,298.15)
    
    return(Tpo)
