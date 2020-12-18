"""
Sutherland's formula for propellant and air dynamic viscosity (muP), furthermore kP
July 2020
A. Takken
"""

#Crane Company (1988) Flow of fluids through valves,fittings, and pipe
#["Nitrogen","Water","Ammonia","Hydrogen"]

import math

def fun1(propellant,Tp):
    
    ### muP
    mu0 = [17.81e-6,0.0000,9.82e-6,8.76e-6]                     #Sutherland's constants
    Ts = [111,0.0000,370,72]                                    #Sutherland's constants
    T0 = [300.55,0.0000,293.15,293.85]                          #Sutherland's constants
       
    if propellant == 1: #Water
        if Tp <= 393.36: #at 2 bar
            muP = 0.014075241*math.exp(-0.016737*Tp)            #[Pa s], dynamic viscosity, by NIST & Excel
        else:           
            muP = 4.06056e-8*Tp-3.00963e-6                      #[Pa s], dynamic viscosity, by NIST & Excel
    elif propellant in [0,2,3]: #Nitrogen,-,Ammonia,Hydrogen
        muP = (mu0[propellant]*(T0[propellant]+Ts[propellant])  #[Pa s], dynamic viscosity, by Sutherland
        /(Tp+Ts[propellant])*(Tp/T0[propellant])**(3/2))      
    
    ### kP
    a = [5.7326e-5,[7.9792e-4,1.1479e-4],1.32783e-4,4.8422e-4]  #gradient, from NIST & Excel
    b = [0.0086744,[0.36934,-0.017661],-0.014539,0.040601]      #intersect, from NIST & Excel
    
    if propellant in [0,2,3]: #Nitrogen,-,Ammonia,Hydrogen
        kP = a[propellant]*Tp + b[propellant]                        
    elif propellant == 1:   #Water
        if Tp <= 393.36: #at 2 bar
            kP = a[propellant][0]*Tp + b[propellant][0]
        else:
            kP = a[propellant][1]*Tp + b[propellant][1]    
              
    return(muP,kP)
