"""
SUB: Insulation properties
July 2020
FROM JJ PREIJDE
"""

###NOTE: emissivity and specific heat capacity of Saffil M-FIL unknown!

def fun1(insulation):

    ###Insulation catalog and identifier
    #0: No insulation
    #1: Saffil M-FIL
    #2: MLI
    
    if insulation == 1:       
        kI = 0.22 #Thermal conductivity [W/mK] (at 1073.15 K), not used because an equation is found in Leenders
        emI = 0.09 #Emissivity [-], from the aluminium foil Leenders wrapped around the insulation
        T_maxI = 2273.15 #Maximum operating temperature [degrees Celsius]
        cI = 1000 #Specific heat capacity [J/kg K]
        rhoI = 100 #Density of the material [kg/m3]
        nameI = "Saffil M-Fil"
    elif insulation == 2:
        kI = 0.1E-3 #Thermal conductivity of possible MLI insulation [W/mK]
        emI = 0.04 #Emittance of possible MLI insulation [-]
        T_maxI = 5000 #Maximum operating temperature [degrees Celsius]
        cI = 1009 #Specific heat capacity [J/kg K]
        rhoI = 0.75+0.17 #Density of the material [kg/m3]
        nameI = "MLI"

    return(kI,emI,T_maxI,nameI)