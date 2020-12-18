"""
Transient theoretical model for STP (RAC & nozzle)
July 2020
A. Takken
"""
#######################Explanation & assumptions######################################################################
### Explanation
"""
1)  This theoretical model will calculate the heat transfer in a Receiver-Absorber Cavity (RAC), the pressure in the same area and 
    subsequent thrust generation in a nozzle.
2)  It is for transient use. Of course, when the input does not change for some time (mass flow, power input), the
    model will approach steady state.
3)  Four different propellants are up for selection: nitrogen (g), hydrogen (g), water (l) and ammonia (g). Note that the NIST
    website will be used for the relevant thermal data. Furthermore, possible evaporation of water will also be taken into account,
    although transitional two phase flow is not.
4)  There are two options for the RAC: a conical RAC (STT1, based on Leender's design) and a cylindrical RAC (STT2, based on Takken's
    design).
5)  It is possible to have insulation on the RAC outer wall, to minimize outer wall thermal losses.
6)  Python scripts called in this script are denoted by "SUB_"
"""

### Assumptions
"""
1)  The RAC has a constant temperature throughout the system, as only high-conductive materials (such as copper) and small
    thicknesses will be used. So, the outer wall of the RAC will be assumed to have a similar temperature as the inner wall 
    and the walls next to the channels.
2)  The surroundings of the RAC will be air which is at a windstill and at a uniform temperature.
3)  The power input will be at the inside of the RAC by incoming radiation. The power output will be in five forms: 
    P1. Outer wall convection
    P2. Outer wall radiation
    P3. Conduction through insulation (if applied)     
    P4. Inner wall convection
    P5. Inner wall radiation
    P6. Propellant convection (from RAC to propellant)
    P7. Heating of RAC
4)  When insulation is applied, its conduction (P3) will be taken into account. In that case, the outer wall heat flows will be
    from the insulation outer wall. Naturally, the inner insulation wall temperature will not be equal to the outer insulation
    wall temperature.
5)  The convection from RAC to propellant will be modelled as a serie of X (to be defined) segments, where the propellant has a higher
    temperature the further it goes down the channel(s).
6)  The cylindrical RAC (STT2) will have no inner wall losses, because there is only a small aperture where the irradiation falls
    through. Absorptivity is not also not applicable.
7)  No heat is lost to the surroundings in the transportation of the propellant from RAC to nozzle.
    The channel(s) have a uniform cross-section. 
8)  The nozzle is assumed to have ideal expansion, i.e. p_e = p_a.
9)  There are no heat losses in the nozzle.
10) The model is conservative in the sense that it does not take into account entrance and exit heating in the RAC.
11) The RAC is oriented horizontally, with the aperture opening to the side.
"""



#######################Imports#########################################################################################
### Import SUBscripts
import SUB_NISTandconstants
import SUB_materialproperties
import SUB_insulationproperties
import SUB_RACdimensions
import SUB_viewfactorscone
import SUB_viewfactorscylinder
import SUB_spiral
import SUB_cp
import SUB_nozzle
import SUB_pLoss
import SUB_P123
import SUB_P4
import SUB_P6

### Import Python packages
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

time0 = time.time()



#######################User inputs######################################################################################
### Inputs
PinR = 60.0                 #[W], incoming irradiance power
Peff = 0.50                 #[-], power input efficiency
n_t = 50.0/60               #[h], hours of running
n_i = [0/60.0,70/60.0]      #[h], hours of irradiation
t_step = 1.00               #[s], number of seconds per step

# Ambient properties (air is the surrounding medium)
Tamb = 298.15               #[K], ambient temperature
pamb = 1.01325e5            #[Pa], ambient pressure

# RAC properties
RACtype = 0                 #[-], RAC type. "0" for conical, "1" for cylindrical
material = 0                #[-], RAC material. "See SUB_materialproperties" for the options
TRAC = 298.15               #[K], starting RAC temperature
absoIC = 0.70               #[-], absorbtivity of inner cavity black paint on inner cone area

LcavM = 0.05014             #[m], cavity outer length/height, (excluding possible insulation)
DouterM = 0.031             #[m], cavity base outer diameter (excluding possible insulation)
DinnerM = 0.025             #[m], cavity base inner diameter
Dap = 0.025                 #[m], cavity aperture diameter
phi = 14*np.pi/180          #[rad], cavity half angle (0.0 in case of cylinder)

# Propellant & channel properties
propellant = 0              #[-], propellant. "0" for nitrogen (g), "1" for water (l), "2" for ammonia (g), "3" for hydrogen (g)
pIn = 2.0e5                 #[Pa], starting pressure propellant
Tpi = 298.15                #[K], starting temperature propellant
n_p = [[0,80]]
mdot = [300e-6]
#n_p = [[0,60]]                                         
#mdot = [100e-6]
#n_p = [[22,32],[35,44],[46,57],[60,71],[73,84],[86,96]]         #[min], begin and end time of flow
#mdot = [75e-6,100e-6,125e-6,150e-6,175e-6,200e-6]               #[kg/s], total mass flow  

channellayout = 0           #[-], channel layout. "0" for linear, "1" for spiral 
Dch = 0.0006                #[m], cross-sectional diameter of channel
nch = 8                     #[-], number of channels
pitch = 0.005               #[m], pitch between spiral cycles, BE AWARE THAT Dch*nch SHOULD BE LOWER THAN pitch

# insulation properties
insulation = 0              #[-], insulation. "0" for no insulation, "1" for Saffil M-FIL, "2" for MLI
tI = 40.0e-3                #[m], uniform thickness insulation



#######################Constants, NIST data and material properties#####################################################
g0,R_A,sigma,nameP,NISTP,limitsP,MMP = SUB_NISTandconstants.fun1(propellant)                #Various constants and propellant properties
emM,absoM,T_maxM,NISTM,limitsM,MMM,rhoM,nameM = SUB_materialproperties.fun1(material)       #RAC material properties
if insulation > 0: #if insulation
    kI,emI,T_maxI,nameI = SUB_insulationproperties.fun1(insulation)                         #Insulation properties
    emO = emI                                                                               #[-], outer wall emissivity
    DouterA = DouterM + tI*2                                                                #[m], overall outer diameter
    LcavA = LcavM +tI*2                                                                     #[m], overall length
else: #if no insulation
    nameI = "No insulation"
    kI = 0.0                                                                                #[-]
    emO = emM                                                                               #[-], outer wall emissivity
    DouterA = DouterM                                                                       #[m], overall outer diameter
    LcavA = LcavM                                                                           #[m], overall length



#######################RAC heat transfer#################################################################################
### One-time calculations
# Dimensions RAC
ARACi,ARACo,MRAC = SUB_RACdimensions.fun1(RACtype,DinnerM,DouterM,DouterA,Dap,LcavM,rhoM)   #[-], calculated dimensions
DmeanM = (DouterM+DinnerM)/2                                                                #[m], mean RAC diameter (without insulation)

if channellayout == 0:      #for linear channels
    if RACtype == 0: #cone
        Lch = np.sqrt(LcavM**2+DmeanM**2)                                   #[m], linear channel length (cone)
    elif RACtype == 1: #cylinder
        Lch = LcavM                                                         #[m], linear channel length (cylinder)    
elif channellayout == 1:    #for spiral channels
    Lch = SUB_spiral.fun1(RACtype,LcavM,DmeanM,pitch)                       #[m], spiral channel length
    
print(Lch)

# P4. Inner wall convection
if RACtype == 0: #cone
    Dav = DinnerM/2                                             #[m], cavity average diameter
elif RACtype == 1: #cylinder
    phi = 0.0                                                   #[rad], cavity half angle
    Dav = DinnerM                                               #[m], cavity average diameter
      
LsI = ((4.79*np.cos(phi)**4.43 - 0.37*np.sin(phi)**0.719)*Dav                   #[m], characteristic length for inner convection
      +(1.06*np.cos(phi)**3.24-0.0462*np.sin(phi)**0.286)*Dap
      +(7.07*np.cos(phi)**5.31+0.221*np.sin(phi)**2.43)*LcavM)           

# P5. Inner wall radiation
if RACtype == 0: #for cone
    RlossA,RlossE,F13 = SUB_viewfactorscone.fun1(DinnerM/2,Dap/2,LcavM,absoIC)      #Loss factors for inner wall radiation (cone))
elif RACtype == 1: #for cylinder
    RlossA,RlossE = SUB_viewfactorscylinder.fun1(DinnerM/2,Dap/2,LcavM,absoIC)  #Loss factors for inner wall radiation (cylinder))

# P7. Heating of RAC
Pin = PinR*Peff*(1-RlossA)                  #[W], absorbed incoming radiation
Ach = Lch*np.pi*Dch                         #[m2], segment inner wall area per channel

# OVERRIDES for Leenders' cases
MRAC = 0.066                                                                        #[kg], RAC mass (without insulation)
ARACi = 2.267e-3                                                                    #[m2], inner area RAC
if insulation == 0:
    ARACo = 3.897e-3                                                                #[m2], outer area RAC
elif insulation == 1:
    RACtype = 1                                                                     #[-], cylindrical
    DouterA = 0.11                                                                  #[m], outer diameter with insulation
    LcavA = 0.10                                                                    #[m], overall length with insulation
    ARACo = 2*np.pi*DouterA**2/4+np.pi*DouterA*LcavA-np.pi*(Dap/2)**2               #[m2], outer area RAC insulation

### Loop
PMatrix,TMatrix,iMatrix,pLMatrix,FMatrix,IspMatrix,vRMatrix,ReDMatrix,PrPMatrix = [],[],[],[],[],[],[],[],[]  #starting empty matrices
NM = 0
for i in range(0,int(n_t*3600/t_step)):    
    
    ### P1, P2 and P3. Outer (insulation) wall convection & radiation, combined in conduction (with insulation)
    P1,P2,P3,Tinsu,h123 = SUB_P123.fun1(DouterA,DouterM,RACtype,ARACo,emO,TRAC,Tamb,pamb,g0,R_A,sigma,kI,insulation,LcavA)

    ### P4. Inner wall convection
    P4,h4 = SUB_P4.fun1(LsI,ARACi,TRAC,Tamb,pamb,g0,R_A)
    
    '''     
    MRAC = 0.087
    ARACi = 0.002346
    ARACo = 0.005124
    P4 = P4*3.29/1.71
    P1 = P1*10.06/13.41
    '''
    


    ### P5. Inner wall radiation    
    P5 = emM*sigma*ARACi*(TRAC**4-Tamb**4)*RlossE       

    ### P6. Propellant convection (from RAC to propellant)  
    if i >= n_p[NM][0]/60.0*3600/t_step and i < n_p[NM][1]/60.0*3600/t_step: #propellant flow   
        mdotch = mdot[NM]/nch
        Tpo,P6ch,ReD,PrP,Tb = SUB_P6.fun1(Dch,DmeanM,Lch,Ach,mdotch,Tpi,TRAC,propellant,channellayout,NISTP,limitsP,MMP)                                  

        ### pLoss, F & Isp
        pp = SUB_pLoss.fun1(ReD,Dch,DmeanM,Lch,channellayout,R_A,Tb,mdotch,pIn,MMP)                     #[Pa], pressure loss
        F,Isp,ReT,At,Ae = SUB_nozzle.fun1(pp,mdot[NM],R_A,MMP,Tpo,pamb,g0,propellant,NISTP,limitsP)     #nozzle outputs

        ### Velocity check      
        v = mdotch/(pIn/(R_A/MMP*Tpi))/(1/4*np.pi*Dch**2)                       #[m/s], propellant end velocity
        vmax = 175*(1/(pIn/(R_A/MMP*Tpi)))**0.43                                #[m/s], propellant maximum velocity  
   
    elif i > n_p[NM][1]/60.0*3600/t_step and NM != len(mdot)-1 and i < n_p[NM+1][1]/60.0*3600/t_step: #go to next mass flow
        NM += 1
    else: #no propellant flow
        P6ch,F,Isp,Tpo,pp,v,vmax,ReD,PrP = 0.0,0.0,0.0,Tpi,pIn,0.0,1234.0,0.0,0.0
     
    ### P7. Heating of RAC
    if i < n_i[0]*3600/t_step or i >= n_i[1]*3600/t_step: #no heating
        PinL = 0    
    else:
        PinL = Pin
    P7 = PinL - (P1+P2+P4+P5+P6ch*nch)              #[W], power to heat the RAC
    cpM = SUB_cp.fun1(TRAC,NISTM,limitsM,MMM)       #[J/kg/K], specific heat coefficient at constant pressure (material)
    TRAC = TRAC + P7*t_step/cpM/MRAC                #[K], resulting RAC temperature

    ### Matrix saves    
    iMatrix.append(i*t_step/60.0)
    PMatrix.append([Pin,P1,P2,P3,P4,P5,P6ch*nch,P7])
    TMatrix.append([TRAC,Tinsu,Tpo])
    pLMatrix.append(pp)
    FMatrix.append(F)
    IspMatrix.append(Isp)
    vRMatrix.append(v/vmax)
    ReDMatrix.append(ReD)
    PrPMatrix.append(PrP)
    if i % 20 == 0:
        print(TRAC)



#######################Outputs############################################################################################
plt.figure(1)
plt.plot(iMatrix,[item[0] for item in PMatrix],label="Pin")
plt.plot(iMatrix,[item[1] for item in PMatrix],label="P1")
plt.plot(iMatrix,[item[2] for item in PMatrix],label="P2")
plt.plot(iMatrix,[item[3] for item in PMatrix],label="P3")
plt.plot(iMatrix,[item[4] for item in PMatrix],label="P4")
plt.plot(iMatrix,[item[5] for item in PMatrix],label="P5")
plt.plot(iMatrix,[item[6] for item in PMatrix],label="P6")
plt.plot(iMatrix,[item[7] for item in PMatrix],label="P7")
plt.xlabel("Time [min]")
plt.ylabel("P [W]")
plt.legend() 

plt.figure(2)
plt.plot(iMatrix,[item[0] for item in TMatrix],label="T_RAC from Python")
plt.plot(iMatrix,[item[2] for item in TMatrix],label="Tp")
plt.grid()
axes = plt.gca()
axes.set_xlim([0,120])
axes.set_ylim([250,750])
axes.xaxis.set_major_locator(MultipleLocator(10))
axes.yaxis.set_major_locator(MultipleLocator(50))
axes.xaxis.set_minor_locator(AutoMinorLocator(2))
axes.yaxis.set_minor_locator(AutoMinorLocator(2))
axes.grid(which='major', color='#CCCCCC', linestyle='--')
axes.grid(which='minor', color='#CCCCCC', linestyle=':')
plt.xlabel("Time [min]")
plt.ylabel("T [K]")
plt.legend()  

plt.figure(3)
plt.plot(iMatrix,pLMatrix,label="pp")
plt.xlabel("Time [min]")
plt.ylabel("pp [Pa]")
plt.legend()     

plt.figure(4)
plt.plot(iMatrix,FMatrix,label="F")
plt.xlabel("Time [min]")
plt.ylabel("Thrust [N]")
plt.legend()     

plt.figure(5)
plt.plot(iMatrix,IspMatrix,label="Isp")
plt.xlabel("Time [min]")
plt.ylabel("Isp [s]")
plt.legend()   

plt.show()

print("---------------")
print("RAC type:",RACtype)
print("RAC material:",material)
print("Insulation:",insulation)
print("Propellant:",propellant)
print("Power input:",PinR,"[W] at an efficiency of",Peff)
print("Max RAC temperature:","%.1f" % max([item[0] for item in TMatrix]),"[K]")  
print("Max prop temperature:","%.1f" % max([item[2] for item in TMatrix]),"[K]")
print("Min v/vmax:","%.3f" % min(vRMatrix),"[-], max v/vmax:","%.3f" % max(vRMatrix),"[-]")
print("Min ReD:","%.1f" % min(ReDMatrix),"[-], max ReD:","%.1f" % max(ReDMatrix),"[-]")
print("Min PrP:","%.3f" % min(PrPMatrix),"[-], max PrP:","%.3f" % max(PrPMatrix),"[-]")
#print("Equilibrium throat area:","%.5f" % (At*1e6),"[mm2], equilibrium exit area:","%.5f" % (Ae*1e6),"[mm2]")
print("---------------")



#######################Links############################################################################################
#http://www-eng.lbl.gov/~dw/projects/DW4229_LHC_detector_analysis/calculations/emissivity2.pdf
#https://www.engineeringtoolbox.com/solar-radiation-absorbed-materials-d_1568.html
#http://www.solarmirror.com/fom/fom-serve/cache/43.html
