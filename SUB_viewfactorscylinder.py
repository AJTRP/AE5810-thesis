"""
View factors cylinder & radiation loss
July 2020
A. Takken
"""

import numpy as np

def fun1(r1,r4,h,absoIC):
      
    ### From:
    # https://web.engr.uky.edu/rtl/Catalog/
    
    ### There are four surfaces in the cylinder with aperture
    # Surface 1: bottom
    # Surface 2: side
    # Surface 3: top (without aperture)
    # Surface 4: aperture
    
    # r1 [m]: cylinder radius
    # r4 [m]: aperture radius
    # h [m]: cylinder height
    # absoIC [-], absorptivity inner cavity

    ### Cylinder areas
    A1 = np.pi*r1**2            #[m2], bottom surface area
    A2 = np.pi*2*r1*h           #[m2], side surface area
    A3 = np.pi*(r1**2-r4**2)    #[m2], top surface area (without aperture)
    
    ### View factors   
    #F1-2 (C-79)
    H = h/(2*r1)
    F12 = 2*H*((1+H**2)**(1/2)-H)
    
    #F1-4 (C-41)
    R1 = r1/h
    R2 = r4/h
    X = 1+(1+R2**2)/R1**2
    F14 = 1/2*(X-(X**2-4*(R2/R1)**2)**(1/2))
       
    #F2-2 (C-78)
    H = h/(2*r1)
    F22 = (1+H)-(1+H**2)**(1/2)
    
    #F2-1 (by deduction)
    F21 = 1/2*(1-F22)
    
    #F2-4 (C-82)
    R = r1/r4
    H2 = h/r4
    X2 = H2**2+R**2+1
    F24 = 1/(4*R*H2)*(-H2**2-R**2+1+(X2**2-4*R**2)**(1/2))
    
    #F3-2 (C-83)
    if r1 == r4:
        F32 = 0
    else:
        R = r1/r4
        H = h/r4
        
        F32 = 1/2*(1+1/(R**2-1)*(H*(4*R**2+H**2)**(1/2)-((1+R**2+H**2)**2-4*R**2)**(1/2)))
        
    #Other (by deduction)
    F13 = 1-F14-F12
    F23 = 1-F21-F22-F24
    F31 = 1-F32
    
    ### Calculate cylinder absorption loss  
    Qout1 = A1/(A1+A2)*(1-absoIC)
    Qout2 = A2/(A1+A2)*(1-absoIC)
    Qout3 = 0
    QlostA = 0
    for i in range(0,20):
        QlostA += (Qout1*F14+Qout2*F24)
        
        Qreflin1 = Qout2*F21+Qout3*F31
        Qreflin2 = Qout1*F12+Qout2*F22+Qout3*F32
        Qreflin3 = Qout1*F13+Qout2*F23
        
        Qout1 = Qreflin1*(1-absoIC)
        Qout2 = Qreflin2*(1-absoIC)
        Qout3 = Qreflin3*(1-absoIC)

    ### Calculate cylinder emission loss  
    Qout1 = A1/(A1+A2+A3)
    Qout2 = A2/(A1+A2+A3)
    Qout3 = A3/(A1+A2+A3)
    QlostE = 0
    for i in range(0,20):
        QlostE += (Qout1*F14+Qout2*F24)
        
        Qreflin1 = Qout2*F21+Qout3*F31
        Qreflin2 = Qout1*F12+Qout2*F22+Qout3*F32
        Qreflin3 = Qout1*F13+Qout2*F23
        
        Qout1 = Qreflin1*(1-absoIC)
        Qout2 = Qreflin2*(1-absoIC)
        Qout3 = Qreflin3*(1-absoIC)
        
    return(QlostA,QlostE)
