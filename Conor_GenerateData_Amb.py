import numpy as np
#import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math
from scipy.optimize import minimize
#from openpyxl import Workbook, load_workbook
import pandas as pd
import scipy.stats as stats
import os
import openpyxl
import time
import pyDOE

exp_conditions_gendata=[120, 54.9, 1.55, 825, 0.0962]
KnownParameters_gendata=[16.64, 6.41, 0.3]
Sigma_gendata=[0.03, 0.0165]

def kinetic_model(c,W,u,theta):
    CBA = c[0] #mol/L
    CEtOH = c[1]
    CEB = c[2] #mol/L
    CW = c[3]
    
    tempC = u[0] #oC
    flowuLmin = u[1] #uL/min
    InletC = u[2] #mol/L
    dSphere = u[3]*10**(-6) #um
    
    KP1 = theta[0]
    KP2 = theta[1]
    KW  = theta[2]
    R = 8.314
    pi = math.pi
    Tort = 3.745

    flowLs = flowuLmin*10**(-6)/60 #L/s
    TK=tempC+273.15 #K
    TM=((140+70)/2)+273.15 #K
    rSphere = dSphere/2 #m
    
    TubeDiameter = 1000*10**(-6) #m
    TubeArea = pi*(TubeDiameter/2)**(2) #m2
    SuperVelocity = (flowLs/1000)/TubeArea #m/s
    
    kpermass = math.exp(-KP1-KP2*10000*(1/TK-1/TM)/R) #L/g s
    rhocat = 770 #kg/m3
    porosity = 0.32
    kc = 9.7*10**(-5)*SuperVelocity**0.0125 #m/s
    SAsphere = 4*pi*rSphere**2 #m2
    Voloftube1sphere = TubeArea*dSphere #m3
    ac = SAsphere/Voloftube1sphere #m2/m3
    visEtOH = math.exp(-7.3714+2770/(74.6787+TK)) #mPa s 
    D = (7.4*10**(-8)*TK*(46.07)**0.5)/(100*100*visEtOH*92.5**0.6) #m2/s
    De = D*porosity/Tort #m2/s
    Thiele = rSphere*math.sqrt(kpermass*rhocat/De)
    eta = (3/(Thiele**2))*(Thiele/math.tanh(Thiele)-1)
    Omega = eta/(1+eta*kpermass*rhocat/(kc*ac))  
                   
    rate=(kpermass*CBA*CEtOH)/((1+KW*CW)**2)
    #rate=Omega*kpermass*CBA
    
    dCBAdW = -rate/flowLs
    dCEtOHdW = -rate/flowLs
    dCEBdW = rate/flowLs 
    dCWdW = rate/flowLs 
         
    return [dCBAdW, dCEtOHdW, dCEBdW, dCWdW]
        
#Function to create model predicted outlet concentrations for a given set of experimental conditions and a given set of Kinetic parameters    
def generatepredictions(ExpCond, theta):
    ConcPredicted=np.zeros(2)
    soln = odeint(kinetic_model,[ExpCond[2],17.09-1.6824*ExpCond[2],0,0],[0,ExpCond[4]], mxstep = 3000, args=(ExpCond,theta))
    CBA = soln[1, 0]
    CEB = soln[1, 2]
    ConcPredicted = [CBA, CEB]
    return ConcPredicted
#CheckGen_gendata=generatepredictions(exp_conditions_gendata, KnownParameters_gendata)

def ModelplusError(ExpCond, theta, stdev):
    ModelConc=generatepredictions(ExpCond, theta)
    Error=np.zeros(2)
    for i in range (0, 2):
        Error[i]=np.random.normal(0,stdev[i])
    SimulatedConc=ModelConc+Error
    return SimulatedConc
SimulatedExpResults_gendata=ModelplusError(exp_conditions_gendata, KnownParameters_gendata, Sigma_gendata)
print SimulatedExpResults_gendata