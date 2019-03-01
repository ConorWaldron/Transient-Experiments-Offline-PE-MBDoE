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



def kinetic_model_A(c,W,u,theta):
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
    #KW  = theta[2]
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
                   
    #rate=Omega*(kpermass*CBA*CEtOH)/((1+KW*CW)**2)
    rate=Omega*kpermass*CBA*CEtOH
    
    dCBAdW = -rate/flowLs
    dCEtOHdW = -rate/flowLs
    dCEBdW = rate/flowLs 
    dCWdW = rate/flowLs 
         
    return [dCBAdW, dCEtOHdW, dCEBdW, dCWdW]
    
def kinetic_model_B(c,W,u,theta):
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
def SingleExpGeneratePredictions(SinlgeExp, theta, model):
    ConcPredicted=np.zeros(2)
    soln = odeint(model,[SinlgeExp[2],17.09-1.6824*SinlgeExp[2],0,0],[0,SinlgeExp[4]], mxstep = 3000, args=(SinlgeExp,theta))
    CBA = soln[1, 0]
    CEB = soln[1, 2]
    ConcPredicted = [CBA, CEB]
    return ConcPredicted

def generatepredictions(ExpCond, theta, model):
    ConcPredicted=np.zeros([len(ExpCond),2])
    for i in range (0, len(ExpCond)):
        soln = odeint(model,[ExpCond[i,2],17.09-1.6824*ExpCond[i,2],0,0],[0,ExpCond[i,4]], mxstep = 3000, args=(ExpCond[i,:],theta))
        CBA = soln[1, 0]
        CEB = soln[1, 2]
        ConcPredicted[i,:] = [CBA, CEB]
    return ConcPredicted

def loglikelihood(theta, Measurments, ExpCond, stdev, model):
    ModelPredictions=generatepredictions(ExpCond, theta, model)
    BAPredict=ModelPredictions[:,0]
    EBPredict=ModelPredictions[:,1]
    BAmeasured=Measurments[:,0]
    EBmeasured=Measurments[:,1]
    rho_1 = BAmeasured - BAPredict
    rho_2 = EBmeasured - EBPredict
    rho = (rho_1/stdev[0])**2+(rho_2/stdev[1])**2
    residuals = np.sum(rho)
    neg_loglikelihood = math.log(2*math.pi) + 0.5 * (math.log(stdev[0]**2) + math.log(stdev[1]**2)) + 0.5 * residuals
    obj_fun = 0.5 * residuals
    return obj_fun

def parameter_estimation(thetaguess, measureddata, inputconditions, stdev, model):
    new_estimate = minimize(loglikelihood, thetaguess, method = 'Nelder-Mead', options = {'maxiter':2000}, args=(measureddata, inputconditions, stdev, model))
    #print "preformed minimisation of log liklihood"
    return new_estimate

def chisquare_test(conf_level, dof):
    ref_chisquare = stats.chi2.ppf((conf_level),dof)
    return ref_chisquare

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

def sensitivity(OneExp, epsilon, TrueParams, model):
    KPMatrix=perturbation(epsilon,TrueParams)
    PredictedValues= np.zeros([len(KPMatrix),4])
    PredictedMeasurable= np.zeros([len(KPMatrix),2])
    for i in range(len(KPMatrix)):
        Solution = odeint(model,[OneExp[2], 17.09-1.6824*OneExp[2], 0, 0], [0,OneExp[4]], mxstep = 3000, args=(OneExp, KPMatrix[i,:]))
        PredictedValues[i,:] = Solution[1,:]
        PredictedMeasurable[i,0]=Solution[1,0]
        PredictedMeasurable[i,1]=Solution[1,2]
    sensitivity_matrix = np.zeros([len(TrueParams),2])
    for j in range(len(TrueParams)):
        for k in range(2):
            sensitivity_matrix[j,k] = ((PredictedMeasurable[j,k] - PredictedMeasurable[-1,k])/(epsilon*TrueParams[j]))  #divide by eepslion*theta?
    return sensitivity_matrix 

def information(OneExp,TrueParams, epsilon, stdev, model):
    Fisher=np.zeros([len(TrueParams),len(TrueParams)])
    for j in range(len(stdev)):
        sens=sensitivity(OneExp, epsilon, TrueParams, model)[:,j]
        Fisher = Fisher + (1/(stdev[j]**2)) * np.outer(sens,sens)
    return Fisher 

def obs_Fisher(ExpCond,epsilon,TrueParams, stdev, model):
    obs_information = np.zeros([len(ExpCond),len(TrueParams),len(TrueParams)])
    for j in range(len(ExpCond)):
        obs_information[j,:,:] = information(ExpCond[j,:], TrueParams, epsilon, stdev, model)
    overall_obs_Fisher = np.zeros([len(TrueParams),len(TrueParams)])
    for j in range(len(ExpCond)):
        overall_obs_Fisher = overall_obs_Fisher + obs_information[j,:,:]
    return overall_obs_Fisher                            

def obs_covariance(ExpCond,epsilon,TrueParams, stdev, model):
    obs_variance_matrix = np.linalg.inv(obs_Fisher(ExpCond,epsilon,TrueParams, stdev, model))
    return obs_variance_matrix

def correlation(Covariance):
    correlationmatrix = np.zeros([len(Covariance),len(Covariance)])
    for i in range(len(Covariance)):
        for j in range(len(Covariance)):
            correlationmatrix[i,j] = Covariance[i,j]/(np.sqrt(Covariance[i,i] * Covariance[j,j]))
    return correlationmatrix

def t_test(Covariance, TrueParams, conf_level, dof):  
    t_values = np.zeros(len(TrueParams))
    conf_interval = np.zeros(len(TrueParams))
    for j in range(len(TrueParams)):
        conf_interval[j] = np.sqrt(Covariance[j,j]) * stats.t.ppf((1 - ((1-conf_level)/2)), dof) 
        t_values[j] = TrueParams[j]/(conf_interval[j])
    t_ref = stats.t.ppf((1-(1-conf_level)),dof)
    return conf_interval,t_values,t_ref  

def ExpectedModelDeviation(SinlgeExp, thetaA, thetaB, modelA, modelB):
    ResponsesModelA=SingleExpGeneratePredictions(SinlgeExp, thetaA, modelA)
    ResponsesModelB=SingleExpGeneratePredictions(SinlgeExp, thetaB, modelB)
    Deviation=np.subtract(ResponsesModelA, ResponsesModelB)
    return Deviation
    
def ExpectedVariance(SinlgeExp, thetaA, thetaB, modelA, modelB, MeasrumentVar, epsilon, CovarA, CovarB):      
    SensAatgivenExp=sensitivity(SinlgeExp, epsilon, thetaA, modelA)
    SensBatgivenExp=sensitivity(SinlgeExp, epsilon, thetaB, modelB)
    ModelAVar=np.matmul(np.matmul(np.transpose(SensAatgivenExp),CovarA), (SensAatgivenExp))
    ModelBVar=np.matmul(np.matmul(np.transpose(SensBatgivenExp),CovarB), (SensBatgivenExp))
    ExpVar=2*MeasrumentVar+ModelAVar+ModelBVar
    return ExpVar
                        
def MB_2Models_ObjFunction(SinlgeExp, thetaA, thetaB, modelA, modelB, MeasrumentVar, epsilon, CovarA, CovarB):
    ModelDeviation=ExpectedModelDeviation(SinlgeExp, thetaA, thetaB, modelA, modelB)
    TotalVar=ExpectedVariance(SinlgeExp, thetaA, thetaB, modelA, modelB, MeasrumentVar, epsilon, CovarA, CovarB)
    Objective=np.matmul(np.matmul(np.transpose(ModelDeviation),np.linalg.inv(TotalVar)), ModelDeviation)
    return -Objective

def LatinGenerator(N_variables, N_sample, VarRanges, DiameterUsed, CatMassUsed):
    np.random.seed(1) #This sets the seed for my random functions like the latin generator to make them give the same results everytime.
    sampling=pyDOE.lhs(N_variables, N_sample)
    Designs=np.zeros([N_sample,5])
    Designs[:,3]=DiameterUsed
    Designs[:,4]=CatMassUsed
    for i in range (0, N_variables):
        Designs[:,i]=VarRanges[i,0]+sampling[:,i]*(VarRanges[i,1]-VarRanges[i,0])
    return Designs

def FindBestGuess(Designs, thetaA, thetaB, modelA, modelB, MeasrumentVar, epsilon, CovarA, CovarB):
    BestObjective=100
    BestDesign=0
    for i in range (0, len(Designs)):
        CriteriaValue=MB_2Models_ObjFunction(Designs[i], thetaA, thetaB, modelA, modelB, MeasrumentVar, epsilon, CovarA, CovarB)
        if CriteriaValue < BestObjective:
            #print "We found a better guess"
            BestObjective = CriteriaValue
            BestDesign = i
    BestGuessDesign=Designs[BestDesign,:]
    return BestGuessDesign        

def MD_Optimisation(NewExp, thetaA, thetaB, VarRanges, DiameterUsed, CatMassUsed, modelA, modelB, MeasrumentVar, epsilon, CovarA, CovarB):
    new_design=minimize(MB_2Models_ObjFunction, NewExp, method = 'SLSQP', bounds = ([VarRanges[0,0],VarRanges[0,1]],[VarRanges[1,0],VarRanges[1,1]],[VarRanges[2,0],VarRanges[2,1]],[DiameterUsed,DiameterUsed],[CatMassUsed,CatMassUsed]), options = {'maxiter':10000, 'ftol':1e-20}, args = (thetaA, thetaB, modelA, modelB, MeasrumentVar, epsilon, CovarA, CovarB,))
    return new_design

############################################################################################################################################
############################################################################################################################################

                    #Start of MD code here
                    
############################################################################################################################################
############################################################################################################################################                   
def MD_Buzzi(filename, DiameterUsed, CatMassUsed):
    os.chdir('C:\Users\Conor\OneDrive - University College London\Closed Loop System\Conor Python')
    start=time.time()
    #load data
    exp_data = pd.read_excel(filename)
    performed_exp = exp_data[['Temperature', 'Flowrate', 'Concentration', 'dsphere', 'cat mass']].values   # experimental conditions
    y_m = exp_data[['BAC', 'EBC']].values   # experiment measurements    
    
    #We use a function to make a matrix of experimental conditions
    def Matrix_conditions(performed_exp):
        data = np.zeros([len(performed_exp),len(performed_exp[0,:])])
        for i in range(len(performed_exp)):
            data[i,0] = performed_exp[i,0]
            data[i,1] = performed_exp[i,1]
            data[i,2] = performed_exp[i,2]
            data[i,3] = performed_exp[i,3]
            data[i,4] = performed_exp[i,4]
        return data
    exp_conditions = Matrix_conditions(performed_exp)
    
    #This code only works for two measurments per experiment
    sigma_BAC = 0.03 # estimate of variance
    sigma_EBC = 0.0165
    sigma = [sigma_BAC, sigma_EBC]
    
    ErrorVariance = np.zeros([len(sigma),len(sigma)])
    ErrorVariance[0,0]=sigma_BAC*sigma_BAC
    ErrorVariance[1,1]=sigma_EBC*sigma_EBC
    
    GuessParameters_modelA = [14, 6]
    GuessParameters_modelB = [16, 7, 0.5]
    
    confidence_level=0.95
    dofreedom_A=(((len(y_m)) * 2) - len(GuessParameters_modelA))
    dofreedom_B=(((len(y_m)) * 2) - len(GuessParameters_modelB))
    chisq_ref_A_AfterFact = chisquare_test(confidence_level, dofreedom_A)
    chisq_ref_B_AfterFact = chisquare_test(confidence_level, dofreedom_B)
    
    #Here we preform parameter estimation on the previously collected data
    minimisedlogA=parameter_estimation(GuessParameters_modelA, y_m, exp_conditions, sigma, kinetic_model_A)
    paramsA_AfterFact = minimisedlogA.x
    wt_residualsA_AfterFact = 2*minimisedlogA.fun
    
    minimisedlogB=parameter_estimation(GuessParameters_modelB, y_m, exp_conditions, sigma, kinetic_model_B)
    paramsB_AfterFact = minimisedlogB.x
    wt_residualsB_AfterFact = 2*minimisedlogB.fun
    
    Disturbance=0.01
    Cov_modelA_AfterFact=obs_covariance(exp_conditions, Disturbance, paramsA_AfterFact, sigma, kinetic_model_A)
    Cov_modelB_AfterFact=obs_covariance(exp_conditions, Disturbance, paramsB_AfterFact, sigma, kinetic_model_B)
    
    Corr_modelA_AfrerFact=correlation(Cov_modelA_AfterFact)
    Corr_modelB_AfrerFact=correlation(Cov_modelB_AfterFact)
    
    conf_interval_A_AfterFact,t_values_A_AfterFact,t_ref_A_AfterFact = t_test(Cov_modelA_AfterFact, paramsA_AfterFact, confidence_level, dofreedom_A)
    conf_interval_B_AfterFact,t_values_B_AfterFact,t_ref_B_AfterFact = t_test(Cov_modelB_AfterFact, paramsB_AfterFact, confidence_level, dofreedom_B)
    
    #Here we begin the design of the next experiment
    Variableranges=np.zeros([3,2])
    Variableranges[0,:]=[80, 120]
    Variableranges[1,:]=[15, 60]
    Variableranges[2,:]=[0.9, 1.55]    
    
    #TestDeviation=ExpectedModelDeviation([120,40,1.5,850,0.01], paramsA_AfterFact, paramsB_AfterFact, kinetic_model_A, kinetic_model_B)
    #testSumofVariances=ExpectedVariance([120,40,1.5,850,0.01], paramsA_AfterFact, paramsB_AfterFact, kinetic_model_A, kinetic_model_B, ErrorVariance, Disturbance, Cov_modelA_AfterFact, Cov_modelB_AfterFact)
    #testobj=MB_2Models_ObjFunction([120,40,1.5,850,0.01], paramsA_AfterFact, paramsB_AfterFact, kinetic_model_A, kinetic_model_B, ErrorVariance, Disturbance, Cov_modelA_AfterFact, Cov_modelB_AfterFact)
        
    NumberofExpInScreen = 10000
    Screening=LatinGenerator(3,NumberofExpInScreen,Variableranges, DiameterUsed, CatMassUsed)    
    BestDesignFromScreening=FindBestGuess(Screening, paramsA_AfterFact, paramsB_AfterFact, kinetic_model_A, kinetic_model_B, ErrorVariance, Disturbance, Cov_modelA_AfterFact, Cov_modelB_AfterFact)
                            
    
    MinimisedMBDoE=MD_Optimisation(BestDesignFromScreening, paramsA_AfterFact, paramsB_AfterFact, Variableranges, DiameterUsed, CatMassUsed, kinetic_model_A, kinetic_model_B, ErrorVariance, Disturbance, Cov_modelA_AfterFact, Cov_modelB_AfterFact)
    NewDesign=MinimisedMBDoE.x
    Objvalue=MinimisedMBDoE.fun
    
    #########################################################################################################
    #########################################################################################################
    
    #testing
    
    #########################################################################################################
    #########################################################################################################
    
    
    #testobj_HunterBest=MB_2Models_ObjFunction([120, 15, 1.55, DiameterUsed, CatMassUsed], paramsA_AfterFact, paramsB_AfterFact, kinetic_model_A, kinetic_model_B, ErrorVariance, Disturbance, Cov_modelA_AfterFact, Cov_modelB_AfterFact)
    #testobj_BuzziBest=MB_2Models_ObjFunction(NewDesign, paramsA_AfterFact, paramsB_AfterFact, kinetic_model_A, kinetic_model_B, ErrorVariance, Disturbance, Cov_modelA_AfterFact, Cov_modelB_AfterFact)    
    
    #MinimisedMBDoETestHunter=MD_Optimisation([120, 15, 1.55, DiameterUsed, CatMassUsed], paramsA_AfterFact, paramsB_AfterFact, Variableranges, DiameterUsed, CatMassUsed, kinetic_model_A, kinetic_model_B, ErrorVariance, Disturbance, Cov_modelA_AfterFact, Cov_modelB_AfterFact)
    #NewDesignTestHunter=MinimisedMBDoE.x
    #ObjvalueTestHunter=MinimisedMBDoE.fun
    
    end=time.time()
    runtime=end-start
    print(runtime)
    
    return paramsA_AfterFact, wt_residualsA_AfterFact, chisq_ref_A_AfterFact, Cov_modelA_AfterFact, Corr_modelA_AfrerFact, conf_interval_A_AfterFact, t_values_A_AfterFact, t_ref_A_AfterFact, paramsB_AfterFact, wt_residualsB_AfterFact, chisq_ref_B_AfterFact, Cov_modelB_AfterFact, Corr_modelB_AfrerFact, conf_interval_B_AfterFact, t_values_B_AfterFact, t_ref_B_AfterFact, BestDesignFromScreening, NewDesign    
