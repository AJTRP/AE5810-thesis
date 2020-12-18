"""
SUB: Material properties
July 2020
FROM JJ PREIJDE
"""

def fun1(material):

    ###Material catalog and identifier
    #0: Copper
    #1: Tungsten
    #2: Molybdenum
    #3: Molybdenum with black paint coating
    #4: Molybdenum with white coating
    #5: Tungsten with black coating
    #6: Tungsten with white coating  
    
    if  material == 0:
        kM = 385 #Thermal conductivity [W/mK]
        emM = 0.65 #Emissivity [-] assuming oxide layer (0.70 by Leenders (although 0.77-0.87 in tables, the latter being oxidized copper), 0.63 by other sources, 0.65 by Leenders page 112)
        absoM = 0.98 #Absorptivity [-] assuming oxide layer
        T_maxM = 1084 #Maximum operating temperature [degrees Celsius]
        NISTM = [[17.72891],[28.09870],[-31.25289],[13.97243],[0.068611]] #Specific heat capacity matrix (NIST) [J/kg K]
        limitsM = [1358.0,1e10,1e10] #[K], NIST limits
        MMM = 63.546e-3 #Molar mass [kg/mol]
        rhoM = 8900 #Density of the material [kg/m3]
        sigmayM = 70E6 #Yield strength [Pa]
        nameM = "Copper"
    elif material == 1:
        kM = 173 #Thermal conductivity [W/mK]
        emM = 0.27 #Emissivity [-]
        absoM = 0.60 #Absorptivity [-]
        T_maxM = 3400 #Maximum operating temperature [degrees Celsius]
        NISTM = [[23.95930,-22.57640],[2.639680,90.27980],[1.257750,-44.27150],[-0.254642,7.176630],[-0.048407,-24.09740]] #Specific heat capacity matrix (NIST) [J/kg/K]
        limitsM = [1900.0,3680.0,1e10] #[K], NIST limits
        MMM = 183.84e-3 #Molar mass [kg/mol]
        rhoM = 19600 #Density of the material [kg/m3]
        sigmayM = 550E6 #Yield strength [Pa]
        nameM = "Tungsten"
    elif material == 2:
        kM = 140 #Thermal conductivity [W/mK]
        emM = 0.18 #Emissivity [-]
        absoM = 0.56 #Absorptivity [-]
        T_maxM = 2620 #Maximum operating temperature [degrees Celsius]
        NISTM = [[24.72736,1231.192],[3.960425,-963.4246],[-1.270706,283.7292],[1.153065,-28.04100],[-0.170246,-712.2047]] #Specific heat capacity matrix (NIST) [J/kg/K]
        limitsM = [1900.0,2896.0,1e10] #[K], NIST limits
        MMM = 95.96e-3 #Molar mass [kg/mol]
        rhoM = 10188 #Density of the material [kg/m3]
        sigmayM = 415E6 #Yield strength [Pa]
        nameM = "Molybdenum"
    elif material == 3:
        kM = 140 #Thermal conductivity [W/mK]
        emM = 0.86 #Emissivity [-]
        absoM = 0.96 #Absorptivity [-]
        T_maxM = 2620 #Maximum operating temperature [degrees Celsius]
        NISTM = [[24.72736,1231.192],[3.960425,-963.4246],[-1.270706,283.7292],[1.153065,-28.04100],[-0.170246,-712.2047]] #Specific heat capacity matrix (NIST) [J/kg/K]
        limitsM = [1900.0,2896.0,1e10] #[K], NIST limits
        MMM = 95.96e-3 #Molar mass [kg/mol]
        rhoM = 10188 #Density of the material [kg/m3]
        sigmayM = 415E6 #Yield strength [Pa]
        nameM = "Molybdenum with black paint coating"
    elif material == 4:
        kM = 140 #Thermal conductivity [W/mK]
        emM = 0.88 #Emissivity [-]
        absoM = 0.06 #Absorptivity [-]
        T_maxM = 2620 #Maximum operating temperature [degrees Celsius]
        NISTM = [[24.72736,1231.192],[3.960425,-963.4246],[-1.270706,283.7292],[1.153065,-28.04100],[-0.170246,-712.2047]] #Specific heat capacity matrix (NIST) [J/kg/K]
        limitsM = [1900.0,2896.0,1e10] #[K], NIST limits
        MMM = 95.96e-3 #Molar mass [kg/mol]
        rhoM = 10188 #Density of the material [kg/m3]
        sigmayM = 415E6 #Yield strength [Pa]
        nameM = "Molybdenum with white paint coating"
    elif material == 5:
        kM = 173 #Thermal conductivity [W/mK]
        emM = 0.86 #Emissivity [-]
        absoM = 0.96 #Absorptivity [-]
        T_maxM = 3400 #Maximum operating temperature [degrees Celsius]
        NISTM = [[23.95930,-22.57640],[2.639680,90.27980],[1.257750,-44.27150],[-0.254642,7.176630],[-0.048407,-24.09740]] #Specific heat capacity matrix (NIST) [J/kg/K]
        limitsM = [1900.0,3680.0,1e10] #[K], NIST limits
        MMM = 183.84e-3 #Molar mass [kg/mol]
        rhoM = 19600 #Density of the material [kg/m3]
        sigmayM = 550E6 #Yield strength [Pa]
        nameM = "Tungsten with black paint coating"
    elif material == 6:
        kM = 173 #Thermal conductivity [W/mK]
        emM = 0.88 #Emissivity [-]
        absoM = 0.06 #Absorptivity [-]
        T_maxM = 3400 #Maximum operating temperature [degrees Celsius]
        NISTM = [[23.95930,-22.57640],[2.639680,90.27980],[1.257750,-44.27150],[-0.254642,7.176630],[-0.048407,-24.09740]] #Specific heat capacity matrix (NIST) [J/kg/K]
        limitsM = [1900.0,3680.0,1e10] #[K], NIST limits
        MMM = 183.84e-3 #Molar mass [kg/mol]
        rhoM = 19600 #Density of the material [kg/m3]
        sigmayM = 550E6 #Yield strength [Pa]
        nameM = "Tungsten with white paint coating"

    return(emM,absoM,T_maxM,NISTM,limitsM,MMM,rhoM,nameM)