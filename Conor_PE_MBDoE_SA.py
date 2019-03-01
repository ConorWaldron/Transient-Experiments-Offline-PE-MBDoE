#Code to simulate Sulfuric Acid RXR
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


#This code has been validated for PE of three parameters against gPROMs
#os.chdir('C:/Users/User/Documents/LabVIEW Data/2015 Projects/Combined')
os.chdir('C:\Users\Conor\OneDrive - University College London\Closed Loop System\Conor Python')

def PE_N_paramThree(filename, ObjectiveCriteria):
    start=time.time()
    exp_data = pd.read_excel(filename)
    #performed_exp = exp_data[['Temperature', 'Flowrate', 'Concentration', 'Dsphere', 'Mass']].values   # experimental conditions
    performed_exp = exp_data[['Temperature', 'Flowrate', 'Concentration']].values   # experimental conditions
    y_m = exp_data[['BA', 'EB']].values   # experiment measurements    
    
    #We use a function to make a matrix of experimental conditions
    def Matrix_conditions(performed_exp):
        data = np.zeros([len(performed_exp),len(performed_exp[0,:])])
        for i in range(len(performed_exp)):
            data[i,0] = performed_exp[i,0]
            data[i,1] = performed_exp[i,1]
            data[i,2] = performed_exp[i,2]
        return data
    exp_conditions = Matrix_conditions(performed_exp)
    
    #This code only works for two measurments per experiment
    sigma_BAC = 0.03 # estimate of variance
    sigma_EBC = 0.0165
    sigma = [sigma_BAC, sigma_EBC]
    
    R = 8.314 # gas constant
    n_theta = 2                            # number of unknown model parameters
    n_y = 2                                # number of measured responses
    n_phi = 3                              # number of design variables
    ParametersGuess = [8.9, 7.9] 
    
    def kinetic_model(c,W,u,theta):
        CBA = c[0] #mol/L
        CEtOH = c[1]
        CEB = c[2] #mol/L
        CW = c[3]
        
        tempC = u[0] #oC
        flowuLmin = u[1] #uL/min
        InletC = u[2] #mol/L
        
        KP1 = theta[0]
        KP2 = theta[1]
        R = 8.314
    
        flowLs = flowuLmin*10**(-6)/60 #L/s
        TK=tempC+273.15 #K
        TM=((140+70)/2)+273.15 #K
       
        kpermass = math.exp(-KP1-KP2*10000*(1/TK-1/TM)/R) #L/g s
        rate=kpermass*CBA
           
        dCBAdW = -rate/flowLs
        dCEtOHdW = -rate/flowLs
        dCEBdW = rate/flowLs 
        dCWdW = rate/flowLs 
             
        return [dCBAdW, dCEtOHdW, dCEBdW, dCWdW]
        
    #Function to create model predicted outlet concentrations for a given set of experimental conditions and a given set of Kinetic parameters    
    def generatepredictions(ExpCond, theta):
        ConcPredicted=np.zeros([len(ExpCond),2])
        for i in range (0, len(ExpCond)):
            soln = odeint(kinetic_model,[ExpCond[i,2],17.09-1.6824*ExpCond[i,2],0,0],[0,9.8174*10**(-5)], mxstep = 3000, args=(ExpCond[i,:],theta))
            CBA = soln[1, 0]
            CEB = soln[1, 2]
            ConcPredicted[i,:] = [CBA, CEB]
        return ConcPredicted
    #CheckGen=generatepredictions(exp_conditions, ParametersGuess)
    
     #Function to evaluate the objective function (sum of the residuals divided by st dev for a given set of conducted experiments and a given set of kinetic parameters             
    def loglikelihood(theta, Measurments, ExpCond):
        ModelPredictions=generatepredictions(ExpCond, theta)
        BAPredict=ModelPredictions[:,0]
        EBPredict=ModelPredictions[:,1]
        BAmeasured=y_m[:,0]
        EBmeasured=y_m[:,1]
        rho_1 = BAmeasured - BAPredict
        rho_2 = EBmeasured - EBPredict
        rho = (rho_1/sigma[0])**2+(rho_2/sigma[1])**2
        residuals = np.sum(rho)
        neg_loglikelihood = math.log(2*math.pi) + 0.5 * (math.log(sigma[0]**2) + math.log(sigma[1]**2)) + 0.5 * residuals
        obj_fun = 0.5 * residuals
        return obj_fun
    #CheckLog=loglikelihood(ParametersGuess, y_m, exp_conditions)
    
        #Function to get parameter estimate
    def parameter_estimation(thetaguess, measureddata, inputconditions):
        new_estimate = minimize(loglikelihood, thetaguess, method = 'Nelder-Mead', options = {'maxiter':2000}, args=(measureddata, inputconditions,))
        #print "preformed minimisation of log liklihood"
        return new_estimate
    minimisedlog = parameter_estimation(ParametersGuess, y_m, exp_conditions)
    params = minimisedlog.x
    wt_residuals = 2*minimisedlog.fun
    #print minimisedlog
        
        ## Adequacy test/chisquare #
    alpha = 0.05
    confidence_level = 1 - alpha
    dofreedom = (((len(performed_exp)) * n_y) - len(ParametersGuess))
    def chisquare_test(conf_level, dof):
        ref_chisquare = stats.chi2.ppf((conf_level),dof)
        return ref_chisquare
    chisq_ref = chisquare_test(confidence_level, dofreedom)
    
     ##Create matrix of perturbed parameters to be used later for making information matrix
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
    #PerturbedParameterMatrix=perturbation(Disturbance,estimate)
    
    #Here you need to decide where in the experimental design space to do the sensitivity analysis
    Examplexp=np.zeros(7)
    Examplexp[0]=120
    Examplexp[1]=20
    Examplexp[2]=1.5
    Examplexp[3]=825
    Examplexp[4]=0.0982
    def sensitivity(OneExp, epsilon, TrueParams):
        KPMatrix=perturbation(epsilon,TrueParams)
        PredictedValues= np.zeros([len(KPMatrix),4])
        PredictedMeasurable= np.zeros([len(KPMatrix),2])
        for i in range(len(KPMatrix)):
            Solution = odeint(kinetic_model,[OneExp[2], 17.09-1.6824*OneExp[2], 0, 0], [0,9.8174*10**(-5)], mxstep = 3000, args=(OneExp, KPMatrix[i,:]))
            PredictedValues[i,:] = Solution[1,:]
            PredictedMeasurable[i,0]=Solution[1,0]
            PredictedMeasurable[i,1]=Solution[1,2]
        sensitivity_matrix = np.zeros([len(TrueParams),n_y])
        for j in range(len(TrueParams)):
            for k in range(n_y):
                sensitivity_matrix[j,k] = ((PredictedMeasurable[j,k] - PredictedMeasurable[-1,k])/(epsilon*TrueParams[j]))  #divide by eepslion*theta?
        return sensitivity_matrix  
    #testsens=sensitivity(Examplexp, Disturbance, params)
    
    #Make information matrix for a single experiment
    def information(OneExp,TrueParams, epsilon):
        Fisher=np.zeros([len(TrueParams),len(TrueParams)])
        for j in range(n_y):
            sens=sensitivity(OneExp, epsilon, TrueParams)[:,j]
            Fisher = Fisher + (1/(sigma[j]**2)) * np.outer(sens,sens)
        return Fisher
    #testFisher=information(Examplexp, params, Disturbance)    #Why do I have negative information here!!!
    
    #Here we get Fisher for all N Experiments
    def obs_Fisher(ExpCond,epsilon,TrueParams):
        obs_information = np.zeros([len(ExpCond),len(TrueParams),len(TrueParams)])
        for j in range(len(ExpCond)):
            obs_information[j,:,:] = information(ExpCond[j,:], TrueParams, epsilon)
        overall_obs_Fisher = np.zeros([len(TrueParams),len(TrueParams)])
        for j in range(len(ExpCond)):
            overall_obs_Fisher = overall_obs_Fisher + obs_information[j,:,:]
        return overall_obs_Fisher 
    #testobsFisher=obs_Fisher(performed_exp,Disturbance,params)   #Why do I have negative information here!!!
    
    def obs_covariance(ExpCond,epsilon,TrueParams):
        obs_variance_matrix = np.linalg.inv(obs_Fisher(ExpCond,epsilon,TrueParams))
        return obs_variance_matrix
    Cov=obs_covariance(performed_exp,Disturbance,params)
    
    def correlation(Covariance):
        correlationmatrix = np.zeros([len(Covariance),len(Covariance)])
        for i in range(len(Covariance)):
            for j in range(len(Covariance)):
                correlationmatrix[i,j] = Covariance[i,j]/(np.sqrt(Covariance[i,i] * Covariance[j,j]))
        return correlationmatrix
    Corr=correlation(Cov)
    
    def t_test(Covariance, TrueParams, conf_level, dof):  
        t_values = np.zeros(len(TrueParams))
        conf_interval = np.zeros(len(TrueParams))
        for j in range(len(TrueParams)):
            conf_interval[j] = np.sqrt(Covariance[j,j]) * stats.t.ppf((1 - ((1-conf_level)/2)), dof) 
            t_values[j] = TrueParams[j]/(conf_interval[j])
        t_ref = stats.t.ppf((1-(1-conf_level)),dof)
        return conf_interval,t_values,t_ref  
    ttests=t_test(Cov, params,confidence_level,dofreedom)
    

    #Here we start MBDoE
        
    #expected Fisher with 1 new experiment
    def ExpectedCovFunction(PastCovariance, NewExp, TrueParams, epsilon):
        PastFisher=np.linalg.inv(PastCovariance)
        FisherofDesign=information(NewExp, TrueParams, epsilon)
        TotalExpectedFisher=PastFisher+FisherofDesign
        ExpectedCov=np.linalg.inv(TotalExpectedFisher)
        return ExpectedCov
    #TestExpCov=ExpectedCovFunction(Cov, InitialGuessDesign, params, Disturbance)
        
    def ObjFunction(NewExp, PastCovariance, TrueParams, epsilon, Criteria):
        if Criteria == "A":
            Obj=np.trace(ExpectedCovFunction(PastCovariance, NewExp, TrueParams, epsilon))
        elif Criteria =="D":
            Obj=np.log(np.linalg.det(ExpectedCovFunction(PastCovariance, NewExp, TrueParams, epsilon))) #we take log of determinant
        elif Criteria ==1:
            Obj=ExpectedCovFunction(PastCovariance, NewExp, TrueParams, epsilon)[0,0]
        elif Criteria ==2:
            Obj=ExpectedCovFunction(PastCovariance, NewExp, TrueParams, epsilon)[1,1]    
        else:
            e,v=np.linalg.eig(ExpectedCovFunction(PastCovariance, NewExp, TrueParams, epsilon))    
            Obj=np.max(e)
        return Obj
    #MarcosExp=[140, 7.954527, 1.55]
    #MarcosObFunctionValue=ObjFunction(MarcosExp, Cov, params, Disturbance, ObjectiveCriteria)
    
    Variableranges=np.zeros([3,2])
    Variableranges[0,:]=[80, 140]
    Variableranges[1,:]=[7.5, 30]
    Variableranges[2,:]=[0.9, 1.55]
    
    #need to set a seed so I get the same results everytime so the work is reproducible
    def LatinGenerator(N_variables, N_sample, VarRanges):
        np.random.seed(1) #This sets the seed for my random functions like the latin generator to make them give the same results everytime.
        sampling=pyDOE.lhs(N_variables, N_sample)
        Designs=np.zeros([N_sample,N_variables])
        for i in range (0, N_variables):
            Designs[:,i]=VarRanges[i,0]+sampling[:,i]*(VarRanges[i,1]-VarRanges[i,0])
        return Designs
    NumberofExpInScreen = 10000
    Screening=LatinGenerator(3,NumberofExpInScreen,Variableranges)
            
    def FindBestGuess(Designs, PastCovariance, TrueParams, epsilon, Criteria):
        BestObjective=100
        BestDesign=0
        for i in range (0, len(Designs)): 
            CriteriaValue=ObjFunction(Designs[i,:], PastCovariance, TrueParams, epsilon, Criteria)
            if CriteriaValue < BestObjective:
                #print "We found a better guess"
                BestObjective = CriteriaValue
                BestDesign = i
        BestGuessDesign=Designs[BestDesign,:]
        return BestGuessDesign        
    BestDesignFromScreening=FindBestGuess(Screening, Cov, params, Disturbance, ObjectiveCriteria)
                            
    def MBDoE(NewExp, PastCovariance, TrueParams, epsilon, Criteria):
        new_design=minimize(ObjFunction, NewExp, method = 'SLSQP', bounds = ([Variableranges[0,0],Variableranges[0,1]],[Variableranges[1,0],Variableranges[1,1]],[Variableranges[2,0],Variableranges[2,1]]), options = {'maxiter':10000, 'ftol':1e-20}, args = (PastCovariance, TrueParams, epsilon, Criteria,))
        return new_design    
    MinimisedMBDoE=MBDoE(BestDesignFromScreening, Cov, params, Disturbance, ObjectiveCriteria)
    NewDesign=MinimisedMBDoE.x
    Objvalue=MinimisedMBDoE.fun
        
        
    end=time.time()
    runtime=end-start
    print(runtime)    
    return params, wt_residuals, chisq_ref, Cov, Corr, ttests, NewDesign
    #return Screening, ScreeningResultsCheck, BestDesignFromScreening