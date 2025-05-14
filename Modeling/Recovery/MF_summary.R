# Relationships between fixation (dwell-time advantage) and choices
# creator: CCT

#clear working environment and set up workind directory
rm(list=ls())
setwd("~/Desktop/WorkingProjects/JM_maaDDM/Hunger_across_domains/Modeling/Recovery")

#clear all plots
if(!is.null(dev.list())) dev.off()

#load required libraries
LibraryList <- c("readxl", "ggpubr", "rstatix", "readr", 
                 "tidyr","dplyr","zoo", "ggplot2", "tibble",
                 "tidyverse","ez","lme4", "lmerTest")
invisible(lapply(LibraryList, library, character.only = TRUE))


# load an random simulation to see modelNames and nModel
load('sim_data_v3/sim_1.Rdata')
modelNames <- names(sim_data)
nModel <- length(modelNames)

# functions
## compute bins
calculate_bins <- function(left_time, right_time) {
  fix_sum <- left_time + right_time 
  fix_diff <- left_time - right_time
  if (sum(fix_sum != 0) > 0) {
    bin_sep <- seq(min(fix_diff), max(fix_diff), length.out = 6)
    bins <- cut(fix_diff, breaks = bin_sep, labels = 1:5, include.lowest = TRUE)
    bindata <- bins
  } else {
    bindata <- rep(NaN, length(left_time))
  }
  return(bindata)
}

## add variables: zscores, bins
add_variables <- function(data, s, m) { #data, simulation number, generative model

  # --------------------------------------------------#
  # dwell difference (left vs. right)
  data$XY_dwelldiff <- (data$F_XA+data$F_XB) - (data$F_YA+data$F_YB) 
  data$ΑΒ_dwelldiff <- (data$F_XA+data$F_YA) - (data$F_XB+data$F_YB) 
  data$XY_A_dwelldiff <- data$F_XA - data$F_YA
  data$XY_B_dwelldiff <- data$F_XB - data$F_YB
  
  # vd of each attribute (X vs. Y)
  data$A_vd <- data$X_A - data$Y_A
  data$B_vd <- data$X_B - data$Y_B
  
  # Choice data (ch_X: 0 = choosing y, 1= choosing x; ch_betterA: 1: choosing option with better A)
  data$ch_X <- as.numeric(data$RT > 0)
  data$ch_betterA <- as.numeric((data$RT > 0 & (data$X_A >= data$Y_A)) | (data$RT < 0 & (data$X_A < data$Y_A)))
  # --------------------------------------------------#
  # Initialize z_data
  z_data <- tibble()
  # Get unique subjects
  subjlist <- unique(data$subject)
  # main loop to add variables for each P
  for (n in subjlist) {
    # Create a subset for the current subject
    dat_sub <- data[data$subject == n & data$absRT>0.25,]
    
    # --------------------------------------------------#
    # relabel subject numbers from 1 to N
    dat_sub$subject <- rep(n+100*m+1000*s, length(dat_sub$RT))
    dat_sub$simulatedM <- rep(m, length(dat_sub$RT))
    dat_sub$Nsimulation <- rep(s, length(dat_sub$RT)) 
    
    # --------------------------------------------------#
    # compute the bins
    dat_sub$attA_diff_bin <- # VD (OptXA-OptYA)
      calculate_bins(
        dat_sub$X_A,
        dat_sub$Y_A
      )
    
    dat_sub$attB_diff_bin <- # VD (OptXB-OptYB)
      calculate_bins(
        dat_sub$X_B,
        dat_sub$Y_B
      )
    
    dat_sub$option_fixdiff_bin <- # dwell-time advantage (OptX-OptY)
      calculate_bins(
        dat_sub$F_XA + dat_sub$F_XB,
        dat_sub$F_YA + dat_sub$F_YB
      )
    
    
    dat_sub$attribute_fixdiff_bin <- # dwell-time advantage (AttA-AttB)
      calculate_bins(
        dat_sub$F_XA + dat_sub$F_YA,
        dat_sub$F_XB + dat_sub$F_YB
      )
    
    dat_sub$attributeA_fixdiff_bin <- # dwell-time advantage (AttX_A-AttY_A)
      calculate_bins(
        dat_sub$F_XA, 
        dat_sub$F_YA
      )
    
    dat_sub$attributeB_fixdiff_bin <- # dwell-time advantage (AttX_B-AttY_B)
      calculate_bins(
        dat_sub$F_XB, 
        dat_sub$F_YB
      )
    
    # --------------------------------------------------#
    # Calculate z-scores for each variable
    dat_sub$z_A_vd <- scale(dat_sub$A_vd)[, 1]
    dat_sub$z_B_vd <- scale(dat_sub$B_vd)[, 1]
    dat_sub$z_XY_dwelldiff <- scale(dat_sub$XY_dwelldiff)[, 1]
    dat_sub$z_ΑΒ_dwelldiff <- scale(dat_sub$ΑΒ_dwelldiff)[, 1]
    dat_sub$z_XY_A_dwelldiff <- scale(dat_sub$XY_A_dwelldiff)[, 1]
    dat_sub$z_XY_B_dwelldiff <- scale(dat_sub$XY_B_dwelldiff)[, 1]
    z_data <- rbind(z_data, dat_sub)
  }
  return(z_data)
}


# Call functions to compute dwell-time advantage and zscores for each P
data_bin <- list()

for (m in 1:nModel) { #simulated model
  # Initialize binned data
  data_bin[[modelNames[m]]]<-list()
  temp_bin <- tibble()
  
  for (s in 1:10) {
    # load data
    load(paste0('sim_data_v3/sim_',s, '.Rdata'))
    dat <- sim_data[[modelNames[m]]]
    dat$absRT <- abs(dat$RT)
    par <- sim_param[[modelNames[m]]]
    subjList <- unique(dat$subject)
    
    # add variables
    temp_bin<-add_variables(dat,s,m)
    
    # group data from all simulated data (50 subject *nSims = total subjects)
    data_bin[[modelNames[m]]] <- rbind(data_bin[[modelNames[m]]], temp_bin) 
  }
}
save(modelNames, nModel, data_bin,file = 'MF_summary.RData')


# general assessment (RT is largest when VD is smallest)
load('MF_summary.RData')

# 1. Create an empty data frame to store all lines data
all_bin_results <- data.frame()

# 2. name of x axis 
model_situation <- c('1phi', '2phi', 
                     '2phi_phiA>phiB', '2phi_phiA<phiB') 
x_name <- c('value difference (X vs. Y)','Attention on option X vs. Y', 'Attention on attribute A vs. B', 
            'Attention on X_A vs. Y_A', 'Attention on X_B vs. Y_B') 
var_name <- c('opt_diff_bin', #1
              'option_fixdiff_bin', #2 
              'attribute_fixdiff_bin', #3
              'attributeA_fixdiff_bin', #4
              'attributeB_fixdiff_bin') #5

var_int = 5;
# 3. main loop
for (m in c(1:4)) {
  loadData <- data_bin[[modelNames[m]]]
  
  if (var_int == 3) {
    temp_data <- loadData%>%
      group_by(.data[[var_name[var_int]]]) %>%
      summarise(
        mean_ch = mean(ch_betterA, na.rm = TRUE)) 
  } else {
    temp_data <- loadData%>%
      group_by(.data[[var_name[var_int]]]) %>%
      summarise(
        mean_ch = mean(ch_X, na.rm = TRUE)) 
  }

  
  # Create a temporary data frame
  temp_bin_results <- data.frame(x = c(1:5), y = temp_data$mean_ch, line_id = factor(m))
  
  # Combine with the main data frame
  all_bin_results <- rbind(all_bin_results, temp_bin_results)
}

# Plot dwell-time advantage on choice
p <- ggplot(all_bin_results, aes(x = x, y = y, color = line_id)) +
  geom_line(size = 1) + 
  labs(title = ,
       x = x_name[var_int],
       y = "P(choose X or better A)") +
  scale_color_discrete(name = "simulated model", labels = model_situation) +
  theme_minimal()+
  coord_cartesian(xlim = c(1, 5), ylim = c(0.2, 0.8))

print(p)


# GLM analysis
## Define the models as a list of formulas
GLM_model_formulas <- list(
  GLM1 = ch_X ~ z_A_vd + z_B_vd + z_XY_dwelldiff,
  
  GLM2 = ch_X ~ z_A_vd + z_B_vd + z_XY_A_dwelldiff + z_XY_B_dwelldiff,
  
  GLM3 = ch_betterA ~ z_A_vd + z_B_vd + z_ΑΒ_dwelldiff
  )

## Initialize a list to store all glm results by dataset and model type
glmm_results <- list()
aic_glmm <- data.frame(dataset = character(), glmm = character(), aic = numeric(), stringsAsFactors = FALSE)

## min. number of trials
minTrial <- 15

## Loop through each dataset
for (m in 1:4) {
  # Retrieve the dataset using `get`
  dataset <- data_bin[[modelNames[m]]]
  dataset[((dataset$F_XA+ dataset$F_XB) > 0 |
             (dataset$F_YA + dataset$F_YB) > 0),] # only consider the datapoint where both fixations are evaluated
  
  # Initialize a sublist for the current dataset
  glmm_results[[modelNames[m]]] <- list()
  
  # Loop through each model formula
  for (n_glmm in names(GLM_model_formulas)) {
    formula <- as.formula(paste(deparse(GLM_model_formulas[[n_glmm]]),"+ (1 | subject)"))
    glmm_results[[modelNames[m]]][[n_glmm]] <- list()
    
    # Run the GLMM model and store the result
    glmm_results[[modelNames[m]]][[n_glmm]] <- glmer(
      formula,
      data = dataset,
      family = binomial(link = "logit"),
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 5e5))
    )
    
    # Extract the AIC value
    aic_value <- AIC(glmm_results[[modelNames[m]]][[n_glmm]])
    
    # Append the AIC value to the results data frame
    aic_glmm <- rbind(aic_glmm, data.frame(dataset = modelNames[m], glmm = n_glmm, aic = aic_value))
  }
}

temp_GLMM_coef<-summary(glmm_results[[modelNames[m]]][[n_glmm]])
z_val = temp_GLMM_coef$coefficients[,"z value"]
p = temp_GLMM_coef$coefficients[,"Pr(>|z|)"]


# GLM analysis
## Initialize a list to store all glm results by dataset and model type
aic_glm_ind <- list()
glm_ind_results <- list()
beta_list <- list()
coef_ttest <- list()
for (m in 1:4) {
  beta_list[[modelNames[m]]] <- list()
}
## Loop through each model formula
for (n_GLM in names(GLM_model_formulas)) {
  formula <- GLM_model_formulas[[n_GLM]]
  
  # Loop through each dataset
  for (m in 1:4) {
    glm_ind_results[[modelNames[m]]][[n_GLM]] <- list()
    
    # Retrieve the dataset using `get`
    dataset <- data_bin[[modelNames[m]]]
    dataset[((dataset$F_XA+ dataset$F_XB) > 0 |
               (dataset$F_YA + dataset$F_YB) > 0),] # only consider the datapoint where both fixations are evaluated
    subjlist <- unique(dataset$subject)
    beta_list[[modelNames[m]]][[n_GLM]] <- matrix(NA,length(subjlist),5)
    count <-0
    
    for (n in subjlist) {
      subjData <- dataset[dataset$subject == n, ]
      count <- count+1
      
      # Initialize a sublist for the current dataset
      glm_ind_results[[modelNames[m]]][[n_GLM]][[as.character(s)]]  <- list()
        
      if (sum((subjData$F_XA+subjData$F_XB)>0)>minTrial & 
          sum((subjData$F_YA+subjData$F_YB)>0)>minTrial){
        # Run the GLM model and store the result
        glm_ind_results[[modelNames[m]]][[n_GLM]][[as.character(n)]]  <-
          glm(formula,
              data = subjData,
              family = binomial(link = "logit")
          )
        
        # Extract the AIC value for checking purpose (make sure it looks "normal")
        aic_value <-
          AIC(glm_ind_results[[modelNames[m]]][[n_GLM]][[as.character(n)]])
        
        # store data from GLM2, we care about z_XY_A_dwelldiff and z_XY_B_dwelldiff
        if (n_GLM == "GLM2") {
          coef_ttest[[modelNames[m]]] <- matrix(NA,1,6)
          temp_beta <- glm_ind_results[[modelNames[m]]][[n_GLM]][[as.character(n)]]
          if (temp_beta$converged == 0) {
            beta_list[[modelNames[m]]][[n_GLM]][count,] <- c(NA,NA,NA,NA,NA)
          } else {
            beta_list[[modelNames[m]]][[n_GLM]][count,] <- as.numeric(temp_beta$coefficients)#c(temp_beta$coefficients[["z_XY_A_dwelldiff"]],temp_beta$coefficients[["z_XY_B_dwelldiff"]])
          }
        }
        
        # Append the AIC value to the results data frame
        aic_glm_ind <-
          rbind(aic_glm_ind,
                data.frame(
                  dataset = modelNames[m],
                  model = n_GLM,
                  subjID = s,
                  aic = aic_value
                ))
      } else {
        glm_ind_results[[modelNames[m]]][[n_GLM]][[as.character(n)]]  <- NA
        aic_glm_ind <-
          rbind(aic_glm_ind,
                data.frame(
                  dataset = modelNames[m],
                  model = n_GLM,
                  subjID = n,
                  aic = NA
                ))
      }
    }
    if (n_GLM == "GLM2"){
      valid_coef <- beta_list[[modelNames[m]]][[n_GLM]]
      t_test_comp_dwelldiff <-
        t.test(valid_coef[, 4], valid_coef[, 5], paired = TRUE)
      coef_ttest[[modelNames[m]]] <-
        c(
          as.numeric(t_test_comp_dwelldiff$estimate),
          as.numeric(t_test_comp_dwelldiff$conf.int[1]),
          as.numeric(t_test_comp_dwelldiff$conf.int[2]),
          as.numeric(t_test_comp_dwelldiff$statistic),
          as.numeric(t_test_comp_dwelldiff$parameter),
          t_test_comp_dwelldiff$p.value
        )
    }
  }
}
## Save GLM results
file_name <- paste0("GLMresults_minTrial_", minTrial, ".Rdata")
save(data_bin, GLM_model_formulas, glm_ind_results, aic_glm_ind, glmm_results, beta_list, aic_glmm, coef_ttest, file = file_name)

# test differences between 
