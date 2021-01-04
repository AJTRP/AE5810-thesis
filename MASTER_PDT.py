"""
Preliminary Design Tool (PDT) for STP
December 2020
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
7)  For further explanation, see thesis document at https://github.com/AJTRP/AE5810-thesis/tree/AE5810-documents
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
# Power & run time
PinR = 250.0                #[W], incoming irradiance power
Peff = 1.00                 #[-], power input efficiency
n_t = 95.0/60               #[h], hours of running
n_i = [0.0,95.0/60]         #[h], hours of irradiation [begin,end]
t_step = 1.00               #[s], number of seconds per step

# Ambient properties (air or vacuum is the surrounding medium)
Tamb = 298.15               #[K], ambient temperature
pamb = 1.01325e5            #[Pa], ambient pressure. If set to 0 Pa, no convection losses are assumed

# RAC properties
RACtype = 1                 #[-], RAC type. "0" for conical, "1" for cylindrical
material = 0                #[-], RAC material. "See SUB_materialproperties" for the options
TRAC = 298.15               #[K], starting RAC temperature
absoIC = 0.70               #[-], absorbtivity of inner cavity black paint on inner area

LcavC = 0.028               #[m], channel length along RAC centerline (for channel length calculation)
LcavI = 0.034               #[m], cavity inner length/height (for inner convection loss)
LcavA = 0.038               #[m], overall length RAC
DinnerM = 0.0118            #[m], cavity base inner diameter
DouterM = 0.019             #[m], cavity base outer diameter (excluding possible insulation)
DmeanM = 0.0124             #[m], cavity channel diameter (distance centerline cavity to midpoint channel, times two)
Dap = 0.004                 #[m], cavity aperture diameter
phi = 0.0*np.pi/180         #[rad], cavity half angle (unused in case of cylinder)

# Insulation properties
insulation = 0              #[-], insulation. "0" for no insulation, "1" for Saffil M-FIL, "2" for MLI
tI = 40.0e-3                #[m], uniform thickness insulation

# Propellant properties
propellant = 0              #[-], propellant. "0" for nitrogen (g), "1" for water (l), "2" for ammonia (g), "3" for hydrogen (g)
pIn = 8.16e5                #[Pa], starting pressure propellant
Tpi = 298.15                #[K], starting temperature propellant
mdot = [300e-6]             #[kg/s], total mass flow
n_p = [[0,95]]              #[min], begin and end time of flow

# Channel properties
channellayout = 1           #[-], channel layout. "0" for linear, "1" for spiral 
Dh = 0.0006                 #[m], cross-sectional diameter of channel
nch = 12                    #[-], number of channels
pitch = 0.0016*nch          #[m], pitch between spiral cycles

# Nozzle properties
ksiF = 0.96                 #[-], nozzle quality/efficiency, average of 0.92-1.00 by Sutton
pe_min = 100                #[Pa], minimum nozzle exit pressure



#######################Constants, NIST data and material properties#####################################################
g0,R_A,sigma,nameP,NISTP,limitsP,MMP = SUB_NISTandconstants.fun1(propellant)                #Various constants and propellant properties
emM,absoM,T_maxM,NISTM,limitsM,MMM,rhoM,nameM = SUB_materialproperties.fun1(material)       #RAC material properties
if insulation > 0: #if insulation
    kI,emI,T_maxI,nameI = SUB_insulationproperties.fun1(insulation)                         #Insulation properties
    emO = emI                                                                               #[-], outer wall emissivity
    DouterA = DouterM + tI*2                                                                #[m], overall outer diameter
    LcavA = LcavA +tI*2                                                                     #[m], overall length
else: #if no insulation
    nameI = "No insulation"
    kI = 0.0                                                                                #[-]
    emO = emM                                                                               #[-], outer wall emissivity
    DouterA = DouterM                                                                       #[m], overall outer diameter
    LcavA = LcavA                                                                           #[m], overall length



#######################One-time calculations#################################################################################
# Dimensions RAC
ARACi,ARACo,MRAC = SUB_RACdimensions.fun1(RACtype,DinnerM,DouterM,DouterA,Dap,LcavA,rhoM)   #[-], calculated dimensions

if channellayout == 0:      #for linear channels
    if RACtype == 0: #cone
        Lch = LcavC/np.cos(phi)                                             #[m], linear channel length (cone)
    elif RACtype == 1: #cylinder
        Lch = LcavC                                                         #[m], linear channel length (cylinder)    
elif channellayout == 1:    #for spiral channels
    Lch = SUB_spiral.fun1(RACtype,LcavC,DmeanM,pitch)                       #[m], spiral channel length
Aheat = Lch*np.pi*Dh*nch                    #[m2], heated wall area
Acs = 1/4*np.pi*Dh**2                       #[m2], cross-sectional channel area
    
if RACtype == 0: #cone
    Dav = DinnerM/2                                             #[m], cavity average diameter
elif RACtype == 1: #cylinder
    phi = 0.0                                                   #[rad], cavity half angle
    Dav = DinnerM                                               #[m], cavity average diameter   
LsI = ((4.79*np.cos(phi)**4.43 - 0.37*np.sin(phi)**0.719)*Dav                   #[m], characteristic length for inner convection
      +(1.06*np.cos(phi)**3.24-0.0462*np.sin(phi)**0.286)*Dap
      +(7.07*np.cos(phi)**5.31+0.221*np.sin(phi)**2.43)*LcavI)          

# Radiation loss factors & power input
if RACtype == 0: #for cone
    RlossA,RlossE,F13 = SUB_viewfactorscone.fun1(DinnerM/2,Dap/2,LcavI,absoIC)  #Loss factors for inner wall radiation (cone)
elif RACtype == 1: #for cylinder
    RlossA,RlossE = SUB_viewfactorscylinder.fun1(DinnerM/2,Dap/2,LcavI,absoIC)  #Loss factors for inner wall radiation (cylinder)
Pin = PinR*Peff*(1-RlossA)                  #[W], absorbed incoming radiation

# Nozzle discharge coefficient relation (from [Johnson1998])
Cda = (0.937-0.968)/(0.016-0.008)           #[-], Cd relation slope
Cdb = 0.968-Cda*0.008                       #[-], Cd relation intercept

#######################Override options######################################################################################
ARACi = 0.00121      
ARACo = 0.00274
MRAC = 0.0429
Aheat = 0.0023929



##########################Loop###############################################################################################
PMatrix,TMatrix,iMatrix,pcMatrix,FMatrix,IspMatrix,vRMatrix,ReDMatrix,PrPMatrix = [],[],[],[],[],[],[],[],[]  #starting empty matrices
NM = 0
for i in range(0,int(n_t*3600/t_step)):   

    ### P1, P2 and P3. Outer (insulation) wall convection & radiation, combined in conduction (with insulation)
    P1,P2,P3,Tinsu,h123 = SUB_P123.fun1(DouterA,DouterM,RACtype,ARACo,emO,TRAC,Tamb,pamb,g0,R_A,sigma,kI,insulation,LcavA)

    ### P4. Inner wall convection
    if pamb < 0.5: #if vacuum
        P4,h4, = 0.0,0.0
    else:
        P4,h4 = SUB_P4.fun1(LsI,ARACi,TRAC,Tamb,pamb,g0,R_A)

    ### P5. Inner wall radiation    
    P5 = emM*sigma*ARACi*(TRAC**4-Tamb**4)*RlossE       

    ### Propellant flow
    if i >= n_p[NM][0]/60.0*3600/t_step and i < n_p[NM][1]/60.0*3600/t_step: 
        
        ### P6. Propellant convection (from RAC to propellant)
        mdotch = mdot[NM]/nch   
        Tpo,P6,ReD,PrP,Tb = SUB_P6.fun1(Dh,DmeanM,Lch,Aheat,mdot[NM],mdotch,Tpi,TRAC,propellant,channellayout,NISTP,limitsP,MMP,Acs)                     

        ### pc, F & Isp
        pc = SUB_pLoss.fun1(ReD,Dh,DmeanM,Lch,channellayout,R_A,Tb,mdotch,pIn,MMP)                                          #[Pa], pressure after pressure loss is applied
        F,Isp,ReT,At,Ae,Cd = SUB_nozzle.fun1(pc,mdot[NM],R_A,MMP,Tpo,pamb,g0,propellant,NISTP,limitsP,pe_min,ksiF,Cda,Cdb)  #[-], nozzle outputs

        ### Velocity check      
        v = mdotch/(pIn/(R_A/MMP*Tpi))/Acs                                      #[m/s], propellant end velocity
        vmax = 175*(1/(pIn/(R_A/MMP*Tpi)))**0.43                                #[m/s], propellant maximum velocity  
   
    elif i > n_p[NM][1]/60.0*3600/t_step and NM != len(mdot)-1 and i < n_p[NM+1][1]/60.0*3600/t_step: #go to next mass flow
        NM += 1
    else: #no propellant flow
        P6,F,Isp,Tpo,pc,v,vmax,ReD,PrP,ReT,Cd = 0.0,0.0,0.0,Tpi,pIn,0.0,1234.0,0.0,0.0,0.0,0.0
     
    ### P7. Heating of RAC
    if i < n_i[0]*3600/t_step or i >= n_i[1]*3600/t_step: #no heating
        PinL = 0    
    else:
        PinL = Pin
    P7 = PinL - (P1+P2+P4+P5+P6)                    #[W], power to heat the RAC
    cpM = SUB_cp.fun1(TRAC,NISTM,limitsM,MMM)       #[J/kg/K], specific heat coefficient at constant pressure (material)
    TRAC = TRAC + P7*t_step/cpM/MRAC                #[K], resulting RAC temperature

    ### Matrix saves    
    iMatrix.append(i*t_step/60.0)
    PMatrix.append([Pin,P1,P2,P3,P4,P5,P6,P7])
    TMatrix.append([TRAC,Tinsu,Tpo])
    pcMatrix.append(pc)
    FMatrix.append(F)
    IspMatrix.append(Isp)
    vRMatrix.append(v/vmax)
    ReDMatrix.append(ReD)
    PrPMatrix.append(PrP)

    ### Exceeding material melting temperature
    if TRAC > T_maxM:
        print("RAC temperature exceeds material melting point!")
        break

print(Pin)
print(P1,P2,P4,P5,P6,P7)
print(h123,h4,NUU)

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
plt.grid()
axes = plt.gca()
axes.xaxis.set_major_locator(MultipleLocator(4))
axes.yaxis.set_major_locator(MultipleLocator(25))
axes.xaxis.set_minor_locator(AutoMinorLocator(2))
axes.yaxis.set_minor_locator(AutoMinorLocator(2))
axes.grid(which='major', color='#CCCCCC', linestyle='--')
axes.grid(which='minor', color='#CCCCCC', linestyle=':')
plt.xlabel("Time [min]")
plt.ylabel("P [W]")
plt.legend() 

plt.figure(2)
plt.plot(iMatrix,[item[0] for item in TMatrix],label="T_RAC from Python")
if max(IspMatrix) > 0.0:
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

if max(IspMatrix) > 0.0: 
    plt.figure(3)
    plt.plot(iMatrix,FMatrix,label="F")
    plt.xlabel("Time [min]")
    plt.ylabel("Thrust [N]")
    plt.legend()     
    
    plt.figure(4)
    plt.plot(iMatrix,IspMatrix,label="Isp")
    plt.xlabel("Time [min]")
    plt.ylabel("Isp [s]")
    plt.legend()   
    
    plt.figure(5)
    plt.plot(iMatrix,pcMatrix,label="pc")
    plt.xlabel("Time [min]")
    plt.ylabel("Chamber pressure [Pa]")
    plt.legend()     
plt.show()

print("---------------")
print("RAC type:",RACtype)
print("RAC material:",material)
print("Insulation:",insulation)
print("Propellant:",propellant)
print("Power input:",PinR,"[W] at an efficiency of",Peff*100.0,"[%]")
print("Max RAC temperature:","%.1f" % max([item[0] for item in TMatrix]),"[K]")  
if max(IspMatrix) > 0.0: 
    print("Max prop temperature:","%.1f" % max([item[2] for item in TMatrix]),"[K]")
    print("Thermal efficiency:","%.1f" % (P6/PinR*100.0),"[%]")
    print("Max Isp:","%.1f" % max(IspMatrix),"[s]")
    print("Max thrust:","%.3f" % max(FMatrix),"[N]")
    print("Max pressure loss:","%.1f" % (pIn-min(pcMatrix)),"[Pa]")
    print("Max v/vmax:","%.3f" % max(vRMatrix),"[-]")
    print("Min PrP:","%.3f" % min(PrPMatrix),"[-], max PrP:","%.3f" % max(PrPMatrix),"[-]")
    print("Min channel ReD:","%.1f" % min(ReDMatrix),"[-], max channel ReD:","%.1f" % max(ReDMatrix),"[-]")
    print("Nozzle ReT:","%.1f" % ReT,"[-], giving a discharge coefficient of:","%.3f" % Cd,"[-]")
    print("Throat diameter:","%.3f" % np.sqrt(At*1e6*4/np.pi),"[mm], Exit diameter:","%.3f" % np.sqrt(Ae*1e6*4/np.pi),"[mm]")
print("---------------")
