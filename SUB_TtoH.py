"""
Calculate propellant temperature from enthalpy
July 2020
A. Takken
"""

def fun1(NISTP,limitsP,MMP,Tpo,Tpi):

    def f(i,T):
        return ((NISTP[0][i]*T/1000+NISTP[1][i]*(T/1000)**2/2                               #[J/kg]
                +NISTP[2][i]*(T/1000)**3/3+NISTP[3][i]*(T/1000)**4/4
                -NISTP[4][i]/(T/1000)+NISTP[5][i]-NISTP[6][i])*1000/MMP)
    
    if Tpo < limitsP[0]:
        Tpii = 0
    elif Tpo < limitsP[1]:
        Tpii = 1
    else:
        Tpii = 2
        
    if Tpi < limitsP[0]:
        Tpiii = 0
    elif Tpi < limitsP[1]:
        Tpiii = 1
    else:
        Tpiii = 2     
        
    H = f(Tpii,Tpo) - f(Tpiii,Tpi)
        

    
    return(H)
