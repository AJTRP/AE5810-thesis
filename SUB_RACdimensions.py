"""
SUB: Dimensions
July 2020
A. Takken
"""

import numpy as np

def fun1(RACtype,DinnerM,DouterM,DouterA,Dap,LcavM,rhoM):

    if RACtype == 0: #cone
        ARACi = (np.pi*((1/2*DinnerM)**2-(1/2*Dap)**2)                                          #[m2], cone inner surface area (including lid & aperture)
                +np.pi*(1/2*DinnerM)*np.sqrt(((DinnerM/DouterM*LcavM)**2+(1/2*DinnerM)**2)))     
        ARACo = (np.pi*((1/2*DouterA)**2-(1/2*Dap)**2)                                          #[m2], cone outer surface area (including lid & aperture)
                +np.pi*(1/2*DouterA)*np.sqrt((LcavM**2+(1/2*DouterA)**2)))
        VRAC = (np.pi/3*(((1/2*DouterM)**2*LcavM)-((1/2*DinnerM)**2*DinnerM/DouterM*LcavM))     #[m3], cone volume (including lid & aperture at 1 mm thickness)
                +np.pi*((1/2*DinnerM)**2-(1/2*Dap)**2)*1e-3)
    elif RACtype == 1: #cylinder
        ARACi = 2*np.pi*(1/2*DinnerM)**2-np.pi*(1/2*Dap)**2+np.pi*DinnerM*DinnerM/DouterM*LcavM #[m2], cylinder inner surface area (including lid & aperture)
        ARACo = 2*np.pi*(1/2*DouterA)**2-np.pi*(1/2*Dap)**2+np.pi*DouterA*LcavM                 #[m2], cylinder outer surface area (including lid & aperture)
        VRAC = (np.pi*(1/2*DouterM)**2*LcavM-np.pi*(1/2*DinnerM)**2*DinnerM/DouterM*LcavM       #[m3], cylinder volume (including lid & aperture at 1 mm thickness)
                +np.pi*((1/2*DinnerM)**2-(1/2*Dap)**2)*1e-3)

    MRAC = VRAC*rhoM
    
    return(ARACi,ARACo,MRAC)
