# simulation of mama data with maaDDM (1phi and 2phi)
# creator: CCT
# source of script: simulation_maaDDM2.R from SG

# Subject-level manipulations:
# 1. w
# 2. theta
# 3. phi1
# 4. phi2
# 5. NDT
# 6. a (threshold)
#
# Trial-level manipulations:
# 1. X_A, X_B, Y_A, Y_B (option values for option X and Y and attribute A and B)
# 2. F_X_A, X_B, Y_A, Y_B

#clear working environment
rm(list=ls())

# change working directory
setwd("/Users/chih-chungting/Desktop/WorkingProjects/JM_maaDDM/Hunger_across_domains/Modeling/Recovery")
cur_folder <- getwd()

#load required libraries
library(rtdists)

# task parameters
nTrials = 200
nSubjs = 50
nSims = 10
modelName <- c("M1","M2","M2_phiAlarger","M2_phiBlarger")
nModels = length(modelName) #to test the impact of phi differences

#basic function of maaDDM
maaDDM <- function(n,d,w,theta,phi,a,t0,X_A,Y_A,X_B,Y_B,F_X_A,F_Y_A,F_X_B,F_Y_B){
  
  #define drift rate d
  v <- d* (F_X_A * (w * (X_A - theta * Y_A) + (1-w) * phi * (X_B - theta * Y_B)) + 
    F_Y_A * (w * (theta * X_A - Y_A) + (1-w) * phi * (theta * X_B - Y_B)) + 
    F_X_B * (w * phi * (X_A - theta * Y_A) + (1-w) * (X_B - theta * Y_B)) + 
    F_Y_B * (w * phi * (theta * X_A - Y_A) + (1-w) * (theta * X_B - Y_B)))
  
  #generate decision
  y <- rdiffusion(n,a,v,t0)
  return(y)
}
#basic function of maaDDM with 2 phi
maaDDM_phis <- function(n,d,w,theta,phiA,phiB,a,t0,X_A,Y_A,X_B,Y_B,F_X_A,F_Y_A,F_X_B,F_Y_B){
  
  #define drift rate d
  v <- d* (F_X_A * (w * (X_A - theta * Y_A) + (1-w) * phiB * (X_B - theta * Y_B)) + 
    F_Y_A * (w * (theta * X_A - Y_A) + (1-w) * phiB * (theta * X_B - Y_B)) + 
    F_X_B * (w * phiA * (X_A - theta * Y_A) + (1-w) * (X_B - theta * Y_B)) + 
    F_Y_B * (w * phiA * (theta * X_A - Y_A) + (1-w) * (theta * X_B - Y_B)))
  
  #generate decision
  y <- rdiffusion(n,a,v,t0)
  return(y)
}


# Initialize a listfor the current dataset
sim_data <- list()
sim_param <-list()

# simulation and data saving 
for (s in 11:15) {# number of simulation 
  for (m in 1:nModels) {#simulated model
    print(paste0('current simulation is ',s, ' out of ', nSims, '; model is ',m, ' out of', nModels))
    # Initialize a sublist for the current dataset
    sim_data[[modelName[m]]] <- list()
    sim_param[[modelName[m]]] <- list()

    # Define column names
    column_names <- c("subject", "trial","RT", "X_A", "Y_A", "X_B", "Y_B", 
                      "F_XA","F_YA","F_XB","F_YB",
                      "d","w","theta","phiA","phiB","a","t0")
    
    # Initialize an empty data frame with the specified column names
    tblData <- data.frame(matrix(ncol = 18, nrow = 0))
    colnames(tblData) <- column_names
    
    # range of parameters (similar to Yang and Krajbich, 2022)
    d <- runif(nSubjs,0.2, 0.35) 
    w <- runif(nSubjs,0.55,0.8) 
    theta <- runif(nSubjs,0.2,0.8)
    if (m == 1) { # M1 and M2 is phiA>phiB or phiA<phiB
      phiA <- runif(nSubjs,0.2,0.8)
    } else if (m >=2) {
      phiA <- runif(nSubjs,0.2,0.4) # M2, M3 and M4 is phiA>phiB or phiA<phiB
    }
    phiB <- runif(nSubjs,0.6,0.8)
    a <- abs(rnorm(nSubjs,2.5,0.2))
    t0 <- runif(nSubjs,0.3,0.6)
    
    for (p in 1:nSubjs) {
      
      #attribute values (A/B refer to attributes A/B, X/Y refer to options X/Y)
      X_A = round(runif(nTrials,0,10))
      Y_A = 10 - X_A
      X_B = round(runif(nTrials,0,10))
      Y_B = 10 - X_B
      
      # correct the phiA and phiB for each model
      if (m == 1) {
        phiB = matrix(NaN,1,nSubjs)
      } else if (m == 2 & (p<(nSubjs/2))) { # half participants showed higher phiA and phiB are fully random
          tempPhi <- phiA[p]
          phiA[p] <- phiB[p]
          phiB[p] <- tempPhi
      } else if (m == 3) { #phiA is larger than phiB 
        if (phiA[p]<phiB[p]) {
          tempPhi <- phiA[p]
          phiA[p] <- phiB[p]
          phiB[p] <- tempPhi
        }
      } else if (m == 4) { #phiA is smaller than phiB 
        if (phiA[p]>phiB[p]) {
          tempPhi <- phiA[p]
          phiA[p] <- phiB[p]
          phiB[p] <- tempPhi
        }
      }
      
      for (n in 1:nTrials){
        # define fixations
        F_X_A <- abs(rnorm(1,0.28,0.20))
        F_Y_A <- abs(rnorm(1,0.25,0.20))
        F_X_B <- abs(rnorm(1,0.24,0.20))
        while (any((F_X_A+F_Y_A+F_X_B)>1)) {
          F_X_A <- abs(rnorm(1,0.28,0.20))
          F_Y_A <- abs(rnorm(1,0.25,0.20))
          F_X_B <- abs(rnorm(1,0.24,0.20))
        }
        F_Y_B <- 1 - (F_X_A+F_Y_A+F_X_B)
        # make sure one attribute is better in X and another is better in Y
        if ((X_A[n] > Y_A[n]) & (X_B[n] > Y_B[n])) {
          X_B[n] <- Y_B[n]
          Y_B[n] <- 10-X_B[n]
        } else if ((X_A[n] < Y_A[n]) & (X_B[n] < Y_B[n])) {
          X_B[n] <- Y_B[n]
          Y_B[n] <- 10-X_B[n]
        }
        
        #simulate
        if (m == 1) {
          ch <- maaDDM(1,d[p], w[p],theta[p],phiA[p],a[p],t0[p],X_A[n],Y_A[n],X_B[n],Y_B[n],F_X_A,F_Y_A,F_X_B,F_Y_B)
        } else if (m > 1) { 
          ch <- maaDDM_phis(1,d[p], w[p],theta[p],phiA[p],phiB[p],a[p],t0[p],X_A[n],Y_A[n],X_B[n],Y_B[n],F_X_A,F_Y_A,F_X_B,F_Y_B)
        }
        
        # main DVs
        Choice <- (ch$response=='upper') *1 + (ch$response=='lower') *-1  #upper is choosing option X
        RT_X <- ch$rt*Choice
        
        
        # Create a data frame row with the generated variables
        new_row <- data.frame(subject = p, trial = n, RT = RT_X, X_A = X_A[n], Y_A = Y_A[n], X_B = X_B[n], Y_B = Y_B[n], 
                              F_XA = F_X_A, F_YA = F_Y_A,F_XB = F_X_B, F_YB = F_Y_B,
                              d = d[p], w = w[p], theta = theta[p], phiA = phiA[p], phiB = phiB[p], a = a[p], t0 = t0[p])
        
        # Append the new row to the data frame
        tblData <- rbind(tblData, new_row)
      }
    }
    sim_data[[modelName[m]]] <- tblData
    
    # save parameters
    sim_param[[modelName[m]]] <- list(
      bound = matrix(NA, nrow = 1, ncol = nSubjs),
      ndt = matrix(NA, nrow = 1, ncol = nSubjs),
      weight = matrix(NA, nrow = 1, ncol = nSubjs),
      drift = matrix(NA, nrow = 1, ncol = nSubjs),
      theta = matrix(NA, nrow = 1, ncol = nSubjs),
      phiA = matrix(NA, nrow = 1, ncol = nSubjs),
      phiB = matrix(NA, nrow = 1, ncol = nSubjs)
    )
    
    sim_param[[modelName[m]]]$bound <- a
    sim_param[[modelName[m]]]$ndt   <- t0
    sim_param[[modelName[m]]]$weight <- w
    sim_param[[modelName[m]]]$drift <- d
    sim_param[[modelName[m]]]$theta <- theta
    sim_param[[modelName[m]]]$phiA <- phiA
    sim_param[[modelName[m]]]$phiB <- phiB
  }
  
  foldername <- paste(cur_folder, "/sim_data_v3/", sep="")
  filename <- paste0(foldername, "sim_",s, ".Rdata", sep = "")
  save(sim_data, sim_param, file = filename)
}

