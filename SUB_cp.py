"""
Calculate cp by NIST
July 2020
A. Takken
"""

def fun1(T,NIST,limits,MM):

    if T < limits[0]:
        ii = 0
    elif T < limits[1]:
        ii = 1
    else:
        ii = 2
           
    cp = (NIST[0][ii]+NIST[1][ii]*T/1000+NIST[2][ii]*(T/1000)**2+NIST[3][ii]*(T/1000)**3+NIST[4][ii]/(T/1000)**2)/MM   
                   
    return(cp)