"""
pLoss
July 2020
A. Takken
"""

import numpy as np

def fun1(ReD,Dch,DmeanM,Lch,channellayout,R_A,Tb,mdotch,pIn,MMP):       
    if channellayout == 0: #straight channels
        if ReD < 2300: #laminar flow
            fDB = 64/ReD
        else: #turbulent flow
            fDB = (0.790*np.log(ReD)-1.64)**-2
    elif channellayout == 1: #
        if ReD < 2300: #laminar flow
            fDB = (1+0.14*(Dch/DmeanM)**0.97*ReD)**(1-0.644*(Dch/DmeanM)**0.312)*64/ReD
        else: #turbulent flow
            fDB = 1.216*ReD**-0.25 + 0.116*(Dch/DmeanM)**0.5
      
    pLoss = fDB*8*Lch*R_A*Tb*mdotch**2/(np.pi**2*pIn*MMP*Dch**5)     #[Pa], pressure loss
    
    return(pIn-pLoss)
    