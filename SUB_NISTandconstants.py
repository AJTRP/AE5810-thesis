"""
NIST data & constants
July 2020
A. Takken
"""

def fun1(propellant):
    
    g0 =            9.81                    #[m/s2], gravitational acceleration Earth SL
    R_A =           8.314                   #[J/K/mol], universal gas constant
    sigma =         5.670e-8                #[W/m2/K4], Stefan Boltzmann constant

    PropNameMat =   ["Nitrogen","Water","Ammonia","Hydrogen"]                           #[-], investigated propellants
    NISTMat =       [[[28.98641,19.50583,35.51872],[1.853978,19.88705,1.128728],        #NIST thermal coefficients, excluding "G"!
                      [-9.647459,-8.598535,-0.196103],[16.63537,1.369784,0.014662],
                      [0.000117,0.527601,-4.553760],[-8.671914,-4.935202,-18.97091],
                      [0.0,0.0,0.0]],
                    [[-203.6060,30.09200,41.96426],[1523.290,6.832514,8.622053],
                      [-3196.413,6.793435,-1.499780],[2474.455,-2.534480,0.098119],
                      [3.855326,0.082139	,-11.15764],[-256.5478,-250.8810,-272.1797],
                      [-285.8304,-241.8264,-241.8264]],
                    [[19.99563,52.02427],[49.77119,18.48801],[-15.37599,-3.765128],
                      [1.921168,0.248541],[0.189174,-12.45799],[-53.30667,-85.53895],
                      [-45.89806,-45.89806]],
                    [[33.066178,18.563083,43.413560],[-11.363417,12.257357,-4.293079],
                      [11.432816,-2.859786,1.272428],[-2.772874,0.268238,-0.096876],
                      [-0.158558,1.977990,-20.533862],[-9.980797,-1.147438,-38.515158],
                      [0.0,0.0,0.0]]]
    NISTlimits =    [[500.0,2000.0],                                                    #[K], NIST thermal limits. Note that for ammonia, only 0 and 1 are valid (not 2)
                     [500.0,1700.0],
                     [1400.0,1.0e10],
                     [1000.0,2500.0]]
    MMMat =         [28.0134e-3,18.0153e-3,17.0305e-3,2.01588e-3]                       #[kg/mol], molar masses
    k_pMat =        [0.20,0.5918,0.507,0.1819]                                          #[W/m/K], thermal conductivities (not used)
    gammaPMat =     [1.40,1.33,1.32,1.41]                                               #[-], specific heat ratio
    
        
    return(g0,R_A,sigma,PropNameMat[propellant],NISTMat[propellant],NISTlimits[propellant],MMMat[propellant])
