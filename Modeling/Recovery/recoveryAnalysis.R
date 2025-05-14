# fitting simulated data with maaDDM (1phi and 2phi)
# creator: CCT
# source of script: aDDM.Rmd from JM


#clear working environment
rm(list=ls())

#clear all plots
if(!is.null(dev.list())) dev.off()

#load required libraries
library(rtdists)
library(dfoptim)
library(readxl)
library(tidyr)
library(dplyr)
library(zoo)
library(tibble)
library(readr)
pacman::p_load(tidyverse, ez)
#parallel computing stuff
library(parallel)
library(doParallel)
library(foreach)
numCores <- detectCores()
registerDoParallel(cores=numCores)
#JAGS packages
library(R2jags) #should be put at the start but keep it here for the moment...
library(rtdists) #to be on the safe side when loading workspace (and not executing the first chunk above)

# model fitting information
modelName <- c("M1","M2","M2_phiAlarger","M2_phiBlarger") #maaDDM1phi and maaDDM2phi
nModels = length(modelName)
parametersList.M1= c('mu_bound','sigma_bound','mu_ndt','sigma_ndt',
                     'mu_drift','sigma_drift','mu_weight1','sigma_weight1',
                     'mu_theta','sigma_theta','mu_phy','sigma_phy',
                     'bound','ndt','drift','weight1','theta','phy')

parametersList.M2 = c('mu_bound','sigma_bound','mu_ndt','sigma_ndt',
                      'mu_drift','sigma_drift','mu_weight1','sigma_weight1',
                      'mu_theta','sigma_theta','mu_phy','sigma_phy', 'mu_phy2','sigma_phy2', 
                      'bound','ndt','drift','weight1','theta','phy', 'phy2')

parametersList.M2_phiAlarger <- parametersList.M2 # same list for other M2 (in which the difference between phiA and phiB is controlled)
parametersList.M2_phiBlarger <- parametersList.M2 # same list for other M2 (in which the difference between phiA and phiB is controlled)
modelfile = c("BayesModel_maaDDM.txt", "BayesModel_maaDDM_2Phi.txt", "BayesModel_maaDDM_2Phi.txt", "BayesModel_maaDDM_2Phi.txt")


# Define column names for later modeling results
column_names <- c("simNum", "generatedModel", "fittingModel",  "DIC", "Rhats",
                  "d","w","theta","phiA","phiB","a","t0")

# Initialize an empty data frame with the specified column names
tblResults <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(tblResults) <- column_names


# load data and run model fitting for each model in each simulation
cur_folder <- getwd()
nSims = 10
nFitting = 2

for (s in 1:nSims) {
  foldername <- paste(cur_folder, "/sim_data_v3/", sep = "")
  filename <- paste0(foldername, "sim_", s, ".Rdata",  sep = "")
  load(filename)
  
  DIC <- data.frame(matrix(ncol = nModels, nrow = nModels))
  DICcomparison<- matrix(NA, ncol = nModels, nrow = 1)
  fittingResults <- vector("list", 0)
  
  for (m in 1:nModels) { # In each simulation, we have data simulated by M1, M2, M2_phiAlarger, M2_phiBlarger
    data <- sim_data[[modelName[m]]]
    nSubj <- length(unique(data$subject))
    
    # only include valid trials
    validT <-
      which((is.na(data$F_XA) == 0) |
              (is.na(data$F_YA) == 0) |
              (is.na(data$F_XB) == 0) | 
              (is.na(data$F_YB) == 0))  

    data <- data[validT,]
    validT <- which(abs(data$RT) > 0.35)
    data<-data[validT,]
    N <- length(data$RT)
    S <- length(unique(data$subject))
    
    # get data and initial values together and specify the model
    ddmData <-
      list(
        'N' = N,
        'S' = S,
        'P' = data$subject,
        'attribute1A' = data$X_A,
        'attribute1B' = data$Y_A,
        'attribute2A' = data$X_B,
        'attribute2B' = data$Y_B,
        'RT' = data$RT,
        'fixProp1' = data$F_XA,
        'fixProp2' = data$F_YA,
        'fixProp3' = data$F_XB,
        'fixProp4' = data$F_YB
      )
    
    # initialize space for estimated parameters
    fittingResults[[modelName[m]]] <- list()
    
    for (fm in 1:nFitting) { # fitted model
      if (fm == 1) {
        estPara <- parametersList.M1
      } else if (fm == 2) {
        estPara <- parametersList.M2
      }
      
      #starting values
      T1<-Sys.time()
      nChains <- 5 # or 8
      
      maaDDMresults <-
        jags.parallel(
          ddmData,
          inits = NULL,
          parameters.to.save = estPara,
          
          model.file = modelfile[fm],
          
          working.directory = '../BayesModels_noHungry',
          
          # n.chains = nChains, n.iter = 60000, n.burnin = 30000,n.thin = 12,
          n.chains = nChains,
          n.iter = 10000, #7000
          n.burnin = 1000, #500
          n.thin = 2,
          DIC = TRUE,
          jags.module = c("glm", "dic", "wiener")
        )
      
      # Check duration of computing
      T2 <- Sys.time() #T2-T1
      
      
      # check convergence
      ddmRhats <-
        maaDDMresults$BUGSoutput$summary[, 8] #check with max(Rhats), which should ideally be < 1.01 (1.05 would also be okay)
      conv <- max(ddmRhats)
      
      # check parameter values (group-level)
      DIC[m,fm] <- maaDDMresults$BUGSoutput$DIC
      boundSep <- mean(log(1+exp(maaDDMresults$BUGSoutput$sims.list$mu_bound)))
      ndt <- mean(log(1+exp(maaDDMresults$BUGSoutput$sims.list$mu_ndt)))
      
      weight1 <- mean(pnorm(maaDDMresults$BUGSoutput$sims.list$mu_weight1))
      drift <- mean(log(1+exp(maaDDMresults$BUGSoutput$sims.list$mu_drift)))
      
      theta <- mean(maaDDMresults$BUGSoutput$sims.list$mu_theta)
      phyA <- mean(maaDDMresults$BUGSoutput$sims.list$mu_phy)
      
      if (fm == 1) {
        phyB <- NaN
      } else {
        phyB <- mean(maaDDMresults$BUGSoutput$sims.list$mu_phy2)
      }
      
      # store results
      fittingResults[[modelName[m]]][[modelName[fm]]] <- maaDDMresults
      
      
      # Create a data frame row with the generated variables
      new_row <- data.frame(simNum = s, generatedModel = m, fittingModel = fm, DIC = DIC[m,fm], Rhats = conv, 
                            d = drift,w = weight1, theta = theta, phiA = phyA, phiB = phyB, a = boundSep, t0 = ndt)
      
      # Append the new row to the data frame
      tblResults <- rbind(tblResults, new_row)
    }
    
    # model comparison
    if (m == 1) {
      DICcomparison[m] <- DIC[m,m] < DIC[m,2]   
    } else if (m == 2) {
      DICcomparison[m] <- DIC[m,m] < DIC[m,1]   
    }
  }
  
  res_foldername <- paste(cur_folder, "/sim_result_v3/", sep="")
  res_filename <- paste0(res_foldername, "sim_",s, ".Rdata", sep = "")
  save.image(file = res_filename)
}

