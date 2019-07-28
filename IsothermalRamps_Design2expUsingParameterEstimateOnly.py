# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 18:31:50 2018

@author: Conor
"""

import numpy as np
import scipy.integrate as integrate
import matplotlib
import time
import os
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import minimize
import math
import scipy.stats as stats
import pyDOE
import openpyxl

start=time.time()

FactorialParameters= [9.1152 , 7.98194]   #for paper I used 9.1152 and 7.98194
GuessParameters= [9,8]
ReactorVol = 98.175 #ul
DeadVol = 44.2 #ul
testtime=5869 #s

sigma = [0.03, 0.0165]
n_y=2

def ImportData(filename):
    os.chdir('C:/Users/Conor/OneDrive - University College London/Closed Loop System/Transient - Simple ODEs/Sensitivity to inital guess/Ramp F')
    exp_data = pd.read_excel(filename)
    data = exp_data[['TM', 'CBA', 'CEB', 'FeedConc', 'vo', 'alphav', 'temp0', 'alphaT', 'TL', 'TE', 'Res']].values   # experimental conditions  
    datamatrix=np.zeros([len(data),len(data[0,:])])
    for i in range(len(data)):
        datamatrix[i,0] = data[i,0]
        datamatrix[i,1] = data[i,1]
        datamatrix[i,2] = data[i,2]
        datamatrix[i,3] = data[i,3]
        datamatrix[i,4] = data[i,4]
        datamatrix[i,5] = data[i,5]
        datamatrix[i,6] = data[i,6]
        datamatrix[i,7] = data[i,7]
        datamatrix[i,8] = data[i,8]
        datamatrix[i,9] = data[i,9]
        datamatrix[i,10] = data[i,10]
    return datamatrix  
#SulfuricAcidData=ImportData('Dec 14C Python Format Corrected.xlsx')   #this is later used to make predictions for different data sets with different parameter estiamtes

def ImportManyFiles(filename1, filename2):
    Data1=ImportData(filename1)
    Data2=ImportData(filename2)
    CombinedDatafile=np.concatenate((Data1,Data2), axis=0)
    return CombinedDatafile
#SulfuricAcidData=ImportManyFiles('ImportIsothermal120.xlsx', 'ImportIsothermal140.xlsx')
    
def Times(vo, alphaV, Vr, Vd, tmSec): #Code to calculate the time the sample left and entered the reactor for a given measurment time
    #units are uL for Volume, uL/min for flowrate and SECONDS for time, they are later converted to minutes in the code
    tmMin=tmSec/60.0 #need decimal to avoid truncation errors
    tLeftmin=(vo-np.sqrt(vo*vo-2*alphaV*(vo*tmMin-0.5*alphaV*tmMin*tmMin-Vd)))/alphaV
    tentermin=(vo-np.sqrt(vo*vo-2*alphaV*(vo*tmMin-0.5*alphaV*tmMin*tmMin-(Vr+Vd))))/alphaV  
    tLeftsec=tLeftmin*60
    tentersec=tentermin*60
    restimesec=tLeftsec-tentersec
    return tLeftsec, tentersec, restimesec

def IsothermalBatch(C, t, theta, T):
    KP1 = theta[0]
    KP2 = theta[1]
    rateconstant=np.exp(-KP1-KP2*10000*(1/T-1/378.15)/8.314)
    dCdTow=-rateconstant*C
    return dCdTow

def SimulateSingleIsothermalDataPoint(theta, FeedConc, T, vo, alphaV, tmSec, Vr, Vd):
    residencetime=Times(vo, alphaV, Vr, Vd, tmSec)[2]
    timepoints=np.linspace(0,residencetime,100)
    SimualtedCBA=integrate.odeint(IsothermalBatch, FeedConc, timepoints, args=(theta,T)) 
    CBenzoicAcid=SimualtedCBA[-1]
    CEthylBenzoate=FeedConc-CBenzoicAcid
    return CBenzoicAcid, CEthylBenzoate
#check=SimulateSingleIsothermalDataPoint(GuessParameters, 1.556148, 392.15, 29.8, 0.253, testtime, ReactorVol, DeadVol)

def generatepredictions(ExpCond, theta, Vr, Vd):
    ConcPredicted=np.zeros([len(ExpCond),2])
    for i in range (0, len(ExpCond)):
        soln = SimulateSingleIsothermalDataPoint(theta, ExpCond[i,3], ExpCond[i,6], ExpCond[i,4], ExpCond[i,5], ExpCond[i,0], Vr, Vd)
        CBA = soln[0]
        CEB = soln[1]
        ConcPredicted[i,:] = [CBA, CEB]
    return ConcPredicted
#This code can be used to make prediction profiles for any actually conducted experiment, with any parameter set we like
#CheckGeneratedSol=generatepredictions(SulfuricAcidData, GuessParameters, ReactorVol, DeadVol)

def SaveFile(filenameNoExtension):
    wbwrite=openpyxl.Workbook()
    wswrite=wbwrite.active #opens the first sheet in the wb
    
    #Mainipulating Files
    wswrite['A'+str(1)]="TM" #using the excel numbering system    
    wswrite['B'+str(1)]="Predicted CBA" #using the excel numbering system 
    wswrite['C'+str(1)]="Predicted CEB" #using the excel numbering system 
    wswrite['D'+str(1)]="FeedConc" #using the excel numbering system 
    wswrite['E'+str(1)]="vo" #using the excel numbering system 
    wswrite['F'+str(1)]="alphav" #using the excel numbering system 
    wswrite['G'+str(1)]="temp0" #using the excel numbering system 
    wswrite['H'+str(1)]="alphaT" #using the excel numbering system 
    wswrite['I'+str(1)]="KP1" #using the excel numbering system 
    wswrite['J'+str(1)]="KP2" #using the excel numbering system 
    for i in range (0, len(CheckGeneratedSol)):
        wswrite['A'+str(i+2)]=SulfuricAcidData[i,0] #using the excel numbering system    
        wswrite['B'+str(i+2)]=CheckGeneratedSol[i,0] #using the excel numbering system 
        wswrite['C'+str(i+2)]=CheckGeneratedSol[i,1] #using the excel numbering system 
        wswrite['D'+str(i+2)]=SulfuricAcidData[i,3] #using the excel numbering system 
        wswrite['E'+str(i+2)]=SulfuricAcidData[i,4] #using the excel numbering system 
        wswrite['F'+str(i+2)]=SulfuricAcidData[i,5] #using the excel numbering system    
        wswrite['G'+str(i+2)]=SulfuricAcidData[i,6] #using the excel numbering system 
        wswrite['H'+str(i+2)]=SulfuricAcidData[i,7] #using the excel numbering system 
        wswrite['I'+str(i+2)]=SimulatingParams[0] #using the excel numbering system 
        wswrite['J'+str(i+2)]=SimulatingParams[1] #using the excel numbering system 
    

    #Saving File
    wbwrite.save(filenameNoExtension+".xlsx")# overwrites without warning. So be careful
    return
#SaveFile('PredictionDec14_139oC_RampF')

def loglikelihood(theta, ExpCond, Vr, Vd):
    ModelPredictions=generatepredictions(ExpCond, theta, Vr, Vd)
    BAPredict=ModelPredictions[:,0]
    EBPredict=ModelPredictions[:,1]
    BAmeasured=ExpCond[:,1]
    EBmeasured=ExpCond[:,2]
    rho_1 = BAmeasured - BAPredict
    rho_2 = EBmeasured - EBPredict
    rho = (rho_1/sigma[0])**2+(rho_2/sigma[1])**2
    residuals = np.sum(rho)
    neg_loglikelihood = math.log(2*math.pi) + 0.5 * (math.log(sigma[0]**2) + math.log(sigma[1]**2)) + 0.5 * residuals
    obj_fun = 0.5 * residuals
    return obj_fun
#CheckLog=loglikelihood(KnownParameters, SulfuricAcidData, ReactorVol, DeadVol)

def parameter_estimation(thetaguess, ExpCond, Vr, Vd):
    new_estimate = minimize(loglikelihood, thetaguess, method = 'Nelder-Mead', options = {'maxiter':2000}, args=(ExpCond, Vr, Vd,))
    #print "preformed minimisation of log liklihood"
    return new_estimate
#minimisedlog = parameter_estimation(GuessParameters, SulfuricAcidData, ReactorVol, DeadVol)
#params = minimisedlog.x
#wt_residuals = 2*minimisedlog.fun
#print minimisedlog

## Adequacy test/chisquare #
alpha = 0.05
confidence_level = 1 - alpha
#dofreedom = (((len(SulfuricAcidData)) * n_y) - len(GuessParameters))

def chisquare_test(conf_level, dof):
    ref_chisquare = stats.chi2.ppf((conf_level),dof)
    return ref_chisquare
#chisq_ref = chisquare_test(confidence_level, dofreedom)

Disturbance = 0.01
def perturbation(epsilon,TrueParams):      
    perturbated_matrix = np.zeros([len(TrueParams)+1,len(TrueParams)])
    for j in range(len(TrueParams)):
        for k in range(len(TrueParams)):
            if j==k:
                perturbated_matrix[j,k] = TrueParams[j] * (1 + epsilon)
            else:
                perturbated_matrix[j,k] = TrueParams[k]
    for j in range(len(TrueParams)):
        perturbated_matrix[-1,j] = TrueParams[j]
    return perturbated_matrix
#PerturbedParameterMatrix=perturbation(Disturbance,params)

def sensitivity(OneExp, epsilon, TrueParams, Vr, Vd):
    KPMatrix=perturbation(epsilon,TrueParams)
    PredictedMeasurable= np.zeros([len(KPMatrix),2])
    for i in range(len(KPMatrix)):
        #Solution = odeint(kinetic_model,[OneExp[2], 17.09-1.6824*OneExp[2], 0, 0], [0,9.8174*10**(-5)], mxstep = 3000, args=(OneExp, KPMatrix[i,:]))
        Solution = SimulateSingleIsothermalDataPoint(KPMatrix[i,:], OneExp[3], OneExp[6], OneExp[4], OneExp[5], OneExp[0], Vr, Vd)
        PredictedMeasurable[i,0]=Solution[0]
        PredictedMeasurable[i,1]=Solution[1]
    sensitivity_matrix = np.zeros([len(TrueParams),n_y])
    for j in range(len(TrueParams)):
        for k in range(n_y):
            sensitivity_matrix[j,k] = ((PredictedMeasurable[j,k] - PredictedMeasurable[-1,k])/(epsilon*TrueParams[j]))  #divide by eepslion*theta?
    return sensitivity_matrix  
#testsens=sensitivity(SulfuricAcidData[0,:], Disturbance, params, ReactorVol, DeadVol)

#Make information matrix for a single experiment
def information(OneExp,TrueParams, epsilon, Vr, Vd):
    Fisher=np.zeros([len(TrueParams),len(TrueParams)])
    for j in range(n_y):
        sens=sensitivity(OneExp, epsilon, TrueParams, Vr, Vd)[:,j]
        Fisher = Fisher + (1/(sigma[j]**2)) * np.outer(sens,sens)
    return Fisher          
#testFisher=information(SulfuricAcidData[0,:], params, Disturbance, ReactorVol, DeadVol)
                
#Here we get Fisher for all N Experiments
def obs_Fisher(ExpCond,epsilon,TrueParams, Vr, Vd):
    obs_information = np.zeros([len(ExpCond),len(TrueParams),len(TrueParams)])
    for j in range(len(ExpCond)):
        obs_information[j,:,:] = information(ExpCond[j,:], TrueParams, epsilon, Vr, Vd)
    overall_obs_Fisher = np.zeros([len(TrueParams),len(TrueParams)])
    for j in range(len(ExpCond)):
        overall_obs_Fisher = overall_obs_Fisher + obs_information[j,:,:]
    return overall_obs_Fisher 
#testobsFisher=obs_Fisher(SulfuricAcidData,Disturbance,params, ReactorVol, DeadVol) 

def obs_covariance(ExpCond,epsilon,TrueParams, Vr, Vd):
    obs_variance_matrix = np.linalg.inv(obs_Fisher(ExpCond,epsilon,TrueParams, Vr, Vd))
    return obs_variance_matrix
#Cov=obs_covariance(SulfuricAcidData,Disturbance,params, ReactorVol, DeadVol)

def correlation(Covariance):
    correlationmatrix = np.zeros([len(Covariance),len(Covariance)])
    for i in range(len(Covariance)):
        for j in range(len(Covariance)):
            correlationmatrix[i,j] = Covariance[i,j]/(np.sqrt(Covariance[i,i] * Covariance[j,j]))
    return correlationmatrix
#Corr=correlation(Cov)     

def t_test(Covariance, TrueParams, conf_level, dof):  
    t_values = np.zeros(len(TrueParams))
    conf_interval = np.zeros(len(TrueParams))
    for j in range(len(TrueParams)):
        conf_interval[j] = np.sqrt(Covariance[j,j]) * stats.t.ppf((1 - ((1-conf_level)/2)), dof) 
        t_values[j] = TrueParams[j]/(conf_interval[j])
    t_ref = stats.t.ppf((1-(1-conf_level)),dof)
    return conf_interval,t_values,t_ref  
#ttests=t_test(Cov, params,confidence_level,dofreedom)                          

####################################################################################################################################
#Here I start MBDoE for isothermal flowrate ramps
####################################################################################################################################

#For testing the code I create an experimental design matrix
FixInitialFlowA = 50
FixAlphaVA      = 2
FixFeedConcA    = 1.25
FixTempA        = 130+273.15
FixInitialFlowB = 100
FixAlphaVB      = 7
FixFeedConcB    = 1.5
FixTempB        = 110+273.15


MBDoEDesign    =[FixAlphaVA, FixInitialFlowA, FixFeedConcA, FixTempA, FixAlphaVB, FixInitialFlowB, FixFeedConcB, FixTempB]
#Here we define the variable range for our designed variables.
Variableranges=np.zeros([8,2])
Variableranges[0,:]=[0.00001, 1.25] #For ramp flowrate, 
Variableranges[1,:]=[7.5, 100] #inital flowrates
Variableranges[2,:]=[0.9, 1.55] #inlet feed concentration
Variableranges[3,:]=[343.15, 413.15] #constant temperature in Kelvin
Variableranges[4,:]=Variableranges[0,:] #For ramp flowrate, 
Variableranges[5,:]=Variableranges[1,:] #inital flowrates
Variableranges[6,:]=Variableranges[2,:] #inlet feed concentration
Variableranges[7,:]=Variableranges[3,:] #constant temperature in Kelvin

HPLCTime=7.0*60 
#Now I need a function to simulate a 100 min experiment, sampling every 7 minutes.
#I need to make make this experiment end early if the flowrate reaches 0
def Simulate100minExp(Design, Vr, Vd, SampleTime, theta):
    alphaV=Design[0]
    vo=Design[1]
    FeedConc=Design[2]
    T=Design[3]
    alphaVB=Design[4]
    voB=Design[5]
    FeedConcB=Design[6]
    TB=Design[7]    
    
    vfinal=5
    TimetoEnd=(vo-vfinal)/alphaV
    TimetoEndB=(voB-vfinal)/alphaVB
    
    ExpTime=100*60
    ExpTimeB=100*60
    #If the flowrate reaches 0 before 100 min we need to end the experiment early
    if TimetoEnd<100:
        ExpTime=TimetoEnd*60        
    NumberofSamples=np.trunc(ExpTime/SampleTime)              #Problem this line only works for floats, not np.ndarrays if you use math.trunc...
    MeasurmentTimes=np.linspace(SampleTime, SampleTime*NumberofSamples, int(NumberofSamples)) #I need the Number of Samples to be an integer not a float  
    SimulatedDataSet=np.zeros([len(MeasurmentTimes),11])
    SimulatedDataSet[:,0]=MeasurmentTimes
    timecal=np.zeros([len(MeasurmentTimes),3])
    for i in range(int(NumberofSamples)):
        SimulatedDataSet[i,3]=FeedConc
        SimulatedDataSet[i,4]=vo
        SimulatedDataSet[i,5]=alphaV
        SimulatedDataSet[i,6]=T
        timecal[i,:]=Times(vo, alphaV, Vr, Vd, MeasurmentTimes[i])
        SimulatedDataSet[i,8]=timecal[i,0]
        SimulatedDataSet[i,9]=timecal[i,1]
        SimulatedDataSet[i,10]=timecal[i,2]
    SimulatedResults=generatepredictions(SimulatedDataSet, theta, Vr, Vd)
    for i in range(int(NumberofSamples)):
        SimulatedDataSet[i,1]=SimulatedResults[i,0]
        SimulatedDataSet[i,2]=SimulatedResults[i,1]
        
    if TimetoEndB<100:
        ExpTimeB=TimetoEndB*60
    NumberofSamplesB=np.trunc(ExpTimeB/SampleTime)              #Problem this line only works for floats, not np.ndarrays if you use math.trunc...
    MeasurmentTimesB=np.linspace(SampleTime, SampleTime*NumberofSamplesB, int(NumberofSamplesB)) #I need the Number of Samples to be an integer not a float
    SimulatedDataSetB=np.zeros([len(MeasurmentTimesB),11])
    SimulatedDataSetB[:,0]=MeasurmentTimesB
    timecalB=np.zeros([len(MeasurmentTimesB),3])
    for i in range(int(NumberofSamplesB)):
        SimulatedDataSetB[i,3]=FeedConcB
        SimulatedDataSetB[i,4]=voB
        SimulatedDataSetB[i,5]=alphaVB
        SimulatedDataSetB[i,6]=TB
        timecalB[i,:]=Times(voB, alphaVB, Vr, Vd, MeasurmentTimesB[i])
        SimulatedDataSetB[i,8]=timecalB[i,0]
        SimulatedDataSetB[i,9]=timecalB[i,1]
        SimulatedDataSetB[i,10]=timecalB[i,2]
    SimulatedResultsB=generatepredictions(SimulatedDataSetB, theta, Vr, Vd)
    for i in range(int(NumberofSamplesB)):
        SimulatedDataSetB[i,1]=SimulatedResultsB[i,0]
        SimulatedDataSetB[i,2]=SimulatedResultsB[i,1]
        
    AllSimulatedDataSets=np.concatenate((SimulatedDataSet,SimulatedDataSetB), axis=0)
    return AllSimulatedDataSets
#testexp=Simulate100minExp(MBDoEDesign, ReactorVol, DeadVol, HPLCTime, FactorialParameters)

#We stop using previous information for design
def ExpectedCovFunctionforMultipleDataPoints(Design, Vr, Vd, SampleTime, TrueParams, epsilon ):
    #PastFisher=np.linalg.inv(PastCovariance)
    FutureExpData=Simulate100minExp(Design, Vr, Vd, SampleTime, TrueParams)
    FisherofDesignDataPoint=np.zeros([len(FutureExpData),len(TrueParams),len(TrueParams)])
    TotalNewFisher=np.zeros([len(TrueParams),len(TrueParams)])
    for i in range(len(FutureExpData)):
        FisherofDesignDataPoint[i,:,:]=information(FutureExpData[i], TrueParams, epsilon, Vr, Vd)
        TotalNewFisher=TotalNewFisher+FisherofDesignDataPoint[i,:,:]
    TotalExpectedFisher=TotalNewFisher #+PastFisher we do not include past info here
    #I need a try check loop here for cases where the experimental design is so poor that we get non inversible information matrix
    try:
        ExpectedCov=np.linalg.inv(TotalExpectedFisher)
    except np.linalg.LinAlgError as err:
            if 'Singular matrix' in str(err):
                ExpectedCov=np.zeros([len(TotalExpectedFisher),len(TotalExpectedFisher)]) #here give very bad result so the solver knows this is a poor experiment
                ExpectedCov[0,:]=[0.9,0.56]
                ExpectedCov[1,:]=[0.56,4.37]
                print ('Singular Matrix error occured')
            else:
                print ('some other error occured')
                print (i)
                raise
    return ExpectedCov
#TestExpCov=ExpectedCovFunctionforMultipleDataPoints(MBDoEDesign, ReactorVol, DeadVol, HPLCTime, FactorialParameters, Disturbance )

def ObjFunction(Design, Vr, Vd, SampleTime, TrueParams, epsilon, Criteria):
    if Criteria == "A":
        Obj=np.trace(ExpectedCovFunctionforMultipleDataPoints(Design, Vr, Vd, SampleTime, TrueParams, epsilon))
    elif Criteria =="D":
        Obj=np.log(np.linalg.det(ExpectedCovFunctionforMultipleDataPoints(Design, Vr, Vd, SampleTime, TrueParams, epsilon))) #we take log of determinant
    elif Criteria ==1:
        Obj=ExpectedCovFunctionforMultipleDataPoints(Design, Vr, Vd, SampleTime, TrueParams, epsilon)[0,0]
    elif Criteria ==2:
        Obj=ExpectedCovFunctionforMultipleDataPoints(Design, Vr, Vd, SampleTime, TrueParams, epsilon)[1,1]    
    else:
        e,v=np.linalg.eig(ExpectedCovFunctionforMultipleDataPoints(Design, Vr, Vd, SampleTime, TrueParams, epsilon))    
        Obj=np.max(e)
    return Obj
#TestObj=ObjFunction(MBDoEDesign, ReactorVol, DeadVol, HPLCTime, FactorialParameters, Disturbance, 'D' )

#need to set a seed so I get the same results everytime so the work is reproducible
def LatinGenerator(N_variables, N_sample, VarRanges):
    np.random.seed(1) #This sets the seed for my random functions like the latin generator to make them give the same results everytime.
    sampling=pyDOE.lhs(N_variables, N_sample)
    Designs=np.zeros([N_sample,N_variables])
    for i in range (0, N_variables):
        Designs[:,i]=VarRanges[i,0]+sampling[:,i]*(VarRanges[i,1]-VarRanges[i,0])
    return Designs
NumberofExpInScreen = 10000
Screening=LatinGenerator(len(MBDoEDesign),NumberofExpInScreen,Variableranges)

def FindBestGuess(ScreeningDesigns, Vr, Vd, SampleTime, TrueParams, epsilon, Criteria):
    BestObjective=100
    BestDesign=0
    for i in range (0, len(ScreeningDesigns)): 
        CriteriaValue=ObjFunction(ScreeningDesigns[i,:], Vr, Vd, SampleTime, TrueParams, epsilon, Criteria) #Problem in this line...
        if CriteriaValue < BestObjective:
            #print "We found a better guess"
            BestObjective = CriteriaValue
            BestDesign = i
    BestGuessDesign=ScreeningDesigns[BestDesign,:]
    print ("Finished screening, now starting MBDoE optimisation")
    return BestGuessDesign        
BestDesignFromScreening=FindBestGuess(Screening, ReactorVol, DeadVol, HPLCTime, FactorialParameters, Disturbance, 'D')

def MBDoE(NewExp, Vr, Vd, SampleTime, TrueParams, epsilon, Criteria):
    new_design=minimize(ObjFunction, NewExp, method = 'SLSQP', bounds = ([Variableranges[0,0],Variableranges[0,1]],[Variableranges[1,0],Variableranges[1,1]],[Variableranges[2,0],Variableranges[2,1]],[Variableranges[3,0],Variableranges[3,1]],[Variableranges[4,0],Variableranges[4,1]],[Variableranges[5,0],Variableranges[5,1]],[Variableranges[6,0],Variableranges[6,1]],[Variableranges[7,0],Variableranges[7,1]]), options = {'maxiter':10000, 'ftol':1e-20}, args = (Vr, Vd, SampleTime, TrueParams, epsilon, Criteria,))
    return new_design    
MinimisedMBDoE=MBDoE(BestDesignFromScreening, ReactorVol, DeadVol, HPLCTime, FactorialParameters, Disturbance, 'D')
NewDesign=MinimisedMBDoE.x
Objvalue=MinimisedMBDoE.fun

def CreateSeperateExperimentMatrix(Design, Vr, Vd, SampleTime, theta):
    alphaV=Design[0]
    vo=Design[1]
    FeedConc=Design[2]
    T=Design[3]
    TimetoEnd=vo/alphaV
    ExpTime=100*60
    #If the flowrate reaches 0 before 100 min we need to end the experiment early
    if TimetoEnd<100:
        ExpTime=TimetoEnd*60
    NumberofSamples=np.trunc(ExpTime/SampleTime)              #Problem this line only works for floats, not np.ndarrays if you use math.trunc...
    MeasurmentTimes=np.linspace(SampleTime, SampleTime*NumberofSamples, int(NumberofSamples)) #I need the Number of Samples to be an integer not a float
    SimulatedDataSet=np.zeros([len(MeasurmentTimes),11])
    timecal=np.zeros([len(MeasurmentTimes),3])
    SimulatedDataSet[:,0]=MeasurmentTimes
    for i in range(int(NumberofSamples)):
        SimulatedDataSet[i,3]=FeedConc
        SimulatedDataSet[i,4]=vo
        SimulatedDataSet[i,5]=alphaV
        SimulatedDataSet[i,6]=T
        timecal[i,:]=Times(vo, alphaV, Vr, Vd, MeasurmentTimes[i])
        SimulatedDataSet[i,8]=timecal[i,0]
        SimulatedDataSet[i,9]=timecal[i,1]
        SimulatedDataSet[i,10]=timecal[i,2]
    SimulatedResults=generatepredictions(SimulatedDataSet, theta, Vr, Vd)
    for i in range(int(NumberofSamples)):
        SimulatedDataSet[i,1]=SimulatedResults[i,0]
        SimulatedDataSet[i,2]=SimulatedResults[i,1]
    return SimulatedDataSet

NewExperimentA=CreateSeperateExperimentMatrix(NewDesign[0:4], ReactorVol, DeadVol, HPLCTime, FactorialParameters)
NewExperimentB=CreateSeperateExperimentMatrix(NewDesign[4:8], ReactorVol, DeadVol, HPLCTime, FactorialParameters)


#rivalDesignA=[0.01, 5, 1.55, 378]
#rivalDesignB=[ 0.01, 10, 1.55, 413]
#NewExperimentA=CreateSeperateExperimentMatrix(rivalDesignA, ReactorVol, DeadVol, HPLCTime, FactorialParameters)
#NewExperimentB=CreateSeperateExperimentMatrix(rivalDesignB, ReactorVol, DeadVol, HPLCTime, FactorialParameters)

def SaveFile(filenameNoExtension, experimentdataset):
    wbwrite=openpyxl.Workbook()
    wswrite=wbwrite.active #opens the first sheet in the wb
    
    #Mainipulating Files
    wswrite['A'+str(1)]="TM" #using the excel numbering system    
    wswrite['B'+str(1)]="CBA" #using the excel numbering system 
    wswrite['C'+str(1)]="CEB" #using the excel numbering system 
    wswrite['D'+str(1)]="FeedConc" #using the excel numbering system 
    wswrite['E'+str(1)]="vo" #using the excel numbering system 
    wswrite['F'+str(1)]="alphav" #using the excel numbering system 
    wswrite['G'+str(1)]="temp0" #using the excel numbering system 
    wswrite['H'+str(1)]="alphaT" #using the excel numbering system 
    wswrite['I'+str(1)]="TL" #using the excel numbering system 
    wswrite['J'+str(1)]="TE" #using the excel numbering system 
    wswrite['K'+str(1)]="Res" #using the excel numbering system 
    for i in range (0, len(experimentdataset)):
        wswrite['A'+str(i+2)]=experimentdataset[i,0] #using the excel numbering system    
        wswrite['B'+str(i+2)]=experimentdataset[i,1] #using the excel numbering system 
        wswrite['C'+str(i+2)]=experimentdataset[i,2] #using the excel numbering system 
        wswrite['D'+str(i+2)]=experimentdataset[i,3] #using the excel numbering system 
        wswrite['E'+str(i+2)]=experimentdataset[i,4] #using the excel numbering system 
        wswrite['F'+str(i+2)]=experimentdataset[i,5] #using the excel numbering system    
        wswrite['G'+str(i+2)]=experimentdataset[i,6] #using the excel numbering system 
        wswrite['H'+str(i+2)]=experimentdataset[i,7] #using the excel numbering system 
        wswrite['I'+str(i+2)]=experimentdataset[i,8] #using the excel numbering system 
        wswrite['J'+str(i+2)]=experimentdataset[i,9] #using the excel numbering system 
        wswrite['K'+str(i+2)]=experimentdataset[i,10] #using the excel numbering system 
    

    #Saving File
    wbwrite.save(filenameNoExtension+".xlsx")# overwrites without warning. So be careful
    return
SaveFile('papercheck A', NewExperimentA)
SaveFile('papercheck B', NewExperimentB)


end=time.time()
runtime=end-start
print(runtime)