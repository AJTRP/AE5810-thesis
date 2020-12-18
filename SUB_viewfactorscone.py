"""
View factors cone & radiation loss
July 2020
A. Takken
"""

import numpy as np

def fun1(r1,r3,h,absoIC):
    
    ### From:
    # https://web.engr.uky.edu/rtl/Catalog/

    ### There are three surfaces in the cone
    # Surface 1: cone interior
    # Surface 2: base (without aperture)
    # Surface 3: aperture
    
    # r1 [m]: cone base radius
    # r3 [m], aperture radius
    # h [m]: cone height
    # absoIC [-], absorptivity inner cavity
    
    ### Cone areas
    A1 = np.pi*2*r1*h           #[m2], inner wall surface area
    A2 = np.pi*(r1**2-r3**2)    #[m2], base surface area (without aperture)
    
    ### View factors
    #F1-1 (C-110)
    H = h/r1
    F11 = 1-1/(1+H**2)**(1/2)
    
    #F1-3 (C-117)
    R = r3/r1
    H2 = h/r1
    F13 = R**2/((1+H2**2)**(1/2))
    
    #F1-2 & F2-1 (by deduction)
    F12 = 1-F11-F13
    F21 = 1
    
    ### Calculate cone absorption loss
    Qout1 = 1*(1-absoIC) #all irradiation falls on surface 1
    Qout2 = 0
    QlostA = 0
    
    for i in range(0,20):
        QlostA += Qout1*F13
              
        Qreflin1 = Qout1*F11+Qout2*F21
        Qreflin2 = Qout1*F12

        Qout1 = Qreflin1*(1-absoIC)
        Qout2 = Qreflin2*(1-absoIC)
        
    ### Calculate cone emission loss
    Qout1 = A1/(A1+A2)
    Qout2 = A2/(A1+A2)
    QlostE = 0
    
    for i in range(0,20):
        QlostE += Qout1*F13
              
        Qreflin1 = Qout1*F11+Qout2*F21
        Qreflin2 = Qout1*F12

        Qout1 = Qreflin1*(1-absoIC)
        Qout2 = Qreflin2*(1-absoIC)
             
    return(QlostA,QlostE,F13)
