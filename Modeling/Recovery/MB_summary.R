# Summarize recovery analysis
# creator: CCT
# source of script: aDDM.Rmd from JM

# main purposes:
# 1. model recovery matrix
# 2. parameter recovery (correlation analysis)

#clear working environment
rm(list=ls())
setwd("/Users/chih-chungting/Desktop/WorkingProjects/JM_maaDDM/Hunger_across_domains/Modeling/Recovery")
currentFolder<-getwd()

#clear all plots
if(!is.null(dev.list())) dev.off()


## Model recovery
# load
nSims <- 15
BestModel <-  matrix(0, nrow = 4, ncol = 2)
compM <-c(2,1)
for (s in 1:nSims) {
  res_foldername <- paste(currentFolder, "/sim_result_v3/", sep="")
  res_filename <- paste0(res_foldername, "sim_",s, ".Rdata", sep = "")
  load(res_filename)
  
  for (m in 1:4) { # simulated model
    for (fm in 1:2) { # fitted model
      BestModel[m,fm] <- BestModel[m,fm] + as.numeric(DIC[m,fm] < DIC[m,compM[fm]])
    }
  }
}
save(BestModel,file = paste0(currentFolder, "/recoveryResults_v3.RData"))


# visualization
# confusion matrix
load(file = paste0(currentFolder, "/recoveryResults_v3.RData"))
library(ggplot2)
for (m_int in 2:4) { # model of interest
  conf_matrix_table <- as.table(BestModel)
  Simulation <- factor(c("M1", "M1", "M2", "M2"))
  Fitting <- factor(c("M1", "M2", "M1", "M2"))
  Y <- c(BestModel[1,1]/nSims, BestModel[1,2]/nSims, BestModel[m_int,1]/nSims, BestModel[m_int,2]/nSims)
  df <- data.frame(Simulation, Fitting, Y)
  
  ggplot(data =  df, mapping = aes(x = Simulation, y = Fitting)) +
    geom_tile(aes(fill = Y*100), colour = "white") +
    geom_text(aes(label = sprintf("%1.f", Y*100)), vjust = 1, size = 20) +
    scale_fill_gradient(low = "white", high = "green", guide = "colourbar",name = "Recovery (%)") +
    theme_bw() + 
    theme(legend.position = "right", 
          axis.text = element_text(size = 25),        # Increase axis text size
          axis.title = element_text(size = 25),       # Increase axis title size
          plot.title = element_text(size = 25))       # Increase plot title size)
  
  # inverse matrix
  Y <- c(BestModel[1,1]/(BestModel[1,1]+BestModel[m_int,1]), 
         BestModel[1,2]/(BestModel[1,2]+BestModel[m_int,2]),
         BestModel[m_int,1]/(BestModel[1,1]+BestModel[m_int,1]),
         BestModel[m_int,2]/(BestModel[1,2]+BestModel[m_int,2]))
  df <- data.frame(Simulation, Fitting, Y)
  ggplot(data =  df, mapping = aes(x = Simulation, y = Fitting)) +
    geom_tile(aes(fill = Y*100), colour = "white") +
    geom_text(aes(label = sprintf("%1.f", Y*100)), vjust = 1, size = 20) +
    scale_fill_gradient(low = "white", high = "green", guide = "colourbar",name = "Recovery (%)") +
    theme_bw() + 
    theme(legend.position = "right", 
          axis.text = element_text(size = 25),        # Increase axis text size
          axis.title = element_text(size = 25),       # Increase axis title size
          plot.title = element_text(size = 25))       # Increase plot title size)
}


# check if parameters are converge 
max_rhat <-list()
for (s in 1:nSims) {
  nSims <- 15
  print(s)
  print(nSims)
  res_foldername <- paste(currentFolder, "/sim_result_v3/", sep="")
  res_filename <- paste0(res_foldername, "sim_",s, ".Rdata", sep = "")
  load(res_filename)
  
  for (m in 1:4) {  # simulated model
    if (s == 1) {
      max_rhat[[modelName[m]]]<-matrix(NA,nSims,2)
    }
    for (fm in 1:2) { # fitted model
      tempPara <- sim_param[[modelName[m]]]
      tempFittingRes <- fittingResults[[modelName[m]]][[modelName[fm]]]
      max_rhat[[modelName[m]]][s,fm]<-max(tempFittingRes$BUGSoutput$summary[,8]) #check with max(Rhats), which should ideally be < 1.01 (1.05 would also be okay)
    }
  }
  nSims <- 15
}

# compute the probability of convergence
P_converge <- matrix(NA,4,2);
for (m in 1:4) { 
  for (fm in 1:2) {
    P_converge[m,fm] = sum(max_rhat[[modelName[m]]][,fm]<1.1)/nSims;
  }
}
row.names(P_converge) <- c('M1', 'M2', 'M2_phiAlarger', 'M2_phiBlarger')

save(BestModel, max_rhat, file = paste0(currentFolder, "/recoveryResults_v3.RData"))


# parameter recovery
corrResult <- list()
ParVal_all <- list()

## load parameters
for (m in 1:4) {
  if (m == 1){
    sim_paraName <- names(sim_param$M1)
    est_paraName <- parametersList.M1[13:18]
  } else {
    sim_paraName <- names(sim_param$M2)
    est_paraName <- parametersList.M2[15:21]
  }
  
  corrResult[[modelName[m]]]<-list()
  ParVal_all[[modelName[m]]] <- list() 
  for (p in 1:7) {
    corrResult[[modelName[m]]][[sim_paraName[p]]]<-matrix(NaN,nSims,3)
    ParVal_all[[modelName[m]]][[sim_paraName[p]]]<-matrix(NaN,nSims*50,2)
  }
}

## transform parameters and run correlation analysis (between estimated and actual)
for (s in 1:nSims) {
  res_foldername <- paste(currentFolder, "/sim_result_v3/", sep="")
  res_filename <- paste0(res_foldername, "sim_",s, ".Rdata", sep = "")
  load(res_filename)
  for (m in 1:4) {
    if (m == 1) {
      fm <-1
      nPara <-6
    } else {
      fm <-2
      nPara <-7
      }
    tempPara <- sim_param[[modelName[m]]]
    tempFittingRes <- fittingResults[[modelName[m]]][[modelName[fm]]]
    
    for (p in 1:nPara) {
      
      if (p <=3) {
        estParVal <- log(1+exp(tempFittingRes$BUGSoutput$mean[[est_paraName[p]]])) # drift
        estParVal_median <- log(1+exp(tempFittingRes$BUGSoutput$median[[est_paraName[p]]])) # drift
      } else if (p == 4) {
        estParVal <- pnorm(tempFittingRes$BUGSoutput$mean[[est_paraName[p]]]) #weight
        estParVal_median <- pnorm(tempFittingRes$BUGSoutput$median[[est_paraName[p]]]) #weight
      } else {
        estParVal <- tempFittingRes$BUGSoutput$mean[[est_paraName[p]]]
        estParVal_median <- tempFittingRes$BUGSoutput$median[[est_paraName[p]]]
      }
      
      if (p ==3) { # due to mismatched order of parameters between simParam and estParam
        simParVal <-tempPara[[sim_paraName[4]]]  # drift
      } else if (p ==4) {
        simParVal <-tempPara[[sim_paraName[3]]] # weight
      } else {
        simParVal <-tempPara[[sim_paraName[p]]]
      }
      
      tempCorr <- cor.test(estParVal,simParVal)
      tempResults <- c(as.numeric(tempCorr$estimate),as.numeric(tempCorr$statistic),tempCorr$p.value)
      corrResult[[modelName[m]]][[sim_paraName[p]]][s,]<-tempResults
      
      
      # store all parameters
      tempParList <- data.frame(estParVal,simParVal)
      if (s > 1) {
        tempParList_all <- ParVal_all[[modelName[m]]][[est_paraName[p]]]
        tempParList_all <- rbind(tempParList_all,tempParList)
      } else {
        tempParList_all <- tempParList
      }
      ParVal_all[[modelName[m]]][[est_paraName[p]]]<-tempParList_all
    }
  }
}

save(BestModel, max_rhat, ParVal_all, corrResult,file = paste0(currentFolder, "/recoveryResults_v3.RData"))
        
# are phiA and phiB correlated? 
cor.test(ParVal_all[[modelName[m]]][[est_paraName[6]]][,1], ParVal_all[[modelName[m]]][[est_paraName[7]]][,1], method = "pearson")

# visualization (one figure/parameter for each model: plot individual dots for simulated and estimated param. values; add correlation results)

load(file = paste0(currentFolder, "/recoveryResults_v3.RData"))


# Create the base plot
m_int = 2 # which model (used to generated data) you are interested in
para_plot <- list()
for (p in 1:7) {
  x <- ParVal_all[[modelName[m_int]]][[est_paraName[p]]]$simParVal
  y <- ParVal_all[[modelName[m_int]]][[est_paraName[p]]]$estParVal
  combPara <- data.frame(x = x, y = y) # Combine x and y into a data frame
  
  # add statistical results
  tempCorr <- cor.test(x,y)
  pdf(file = paste0("correlation_",est_paraName[p], ".pdf"),width = 10,height = 7)
  
  para_plot[[est_paraName[p]]] <- ggplot(combPara, aes(x = x, y = y)) +
    geom_point(color = "blue") +  # Scatter plot
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # Correlation line
    geom_abline(intercept = 0, slope = 1, color = "black") +  # 45-degree dashed line
    labs(title = paste("parameter : ", est_paraName[p]),
         x = "Sim. para.", y = "Est. para.") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14),  # Change title font size
      axis.title = element_text(size = 12),  # Change axis titles font size
      axis.text = element_text(size = 15)    # Change axis text font size
    ) + coord_fixed(ratio = 1, xlim = c(min(combPara),max(combPara)), ylim = c(min(combPara),max(combPara))) +
    annotate("text", x = min(combPara), y=max(combPara), label = paste0("r = ", round(as.numeric(tempCorr$estimate),2), ", p = ", sprintf("%.4f",as.numeric(tempCorr$p.value))), 
             hjust = 0, vjust = 1, size = 5, color = "black")
  print(para_plot[[est_paraName[p]]])
  dev.off()
}

#