"""
Spiral calculation
July 2020
Source: https://mathworld.wolfram.com/ConicalSpiral.html
"""

import numpy as np
#import matplotlib.pyplot as plt

def fun1(RACtype,Lcav,D,pitch):

    if RACtype == 0: #conical
        hypot = np.sqrt(Lcav**2+(1/2*D)**2)         #[m], side length of cone
        N = hypot/pitch                             #[-], number of revolutions
        a = N*2.0*np.pi/Lcav                        #[1/m], angular frequency
        
        X,Y,Z = [],[],[] 
        arcL = 0                         
        n = 10000                                   #[-], number of cycles
        
        for i in range(0,n):
            z = Lcav*i/n
            x = (Lcav-z)/Lcav*D/2*np.cos(a*z)
            y = (Lcav-z)/Lcav*D/2*np.sin(a*z)
            
            if i >= 1:   
                arcL += np.sqrt((x-X[-1])**2+(y-Y[-1])**2+(z-Z[-1])**2)
            
            X.append(x)
            Y.append(y)
            Z.append(z)
    
        #plt.figure(1)
        #ax = plt.axes(projection='3d')
        #ax.plot3D(X,Y,Z, 'gray')
    
    if RACtype == 1: #cylindrical
        N = Lcav/pitch
        arcL = np.sqrt(4*np.pi**2*(1/2*D)**2+Lcav**2)*N
    
    return(arcL)
