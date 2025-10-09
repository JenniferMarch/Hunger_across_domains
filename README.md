# Domain Specific Effects of Hunger on Attention and Choice

This is the repository where I am storing the code used to generate **1 behavioural** and **2 modeling** analyses for the manuscript: "Domain Specific Effects of Hunger on Attention and Choice"

## 1 Behavioual Analyses

All analyses are carried out in R using RMarkdown (Version 4.4.1). The analyses are structured using the suffix A_ to D_ and at the beginning of each file further information about the analyses is provided.

**A_Preprocess.Rmd** 
*Input:* to run this file, the original .csv files are needed, which can be requested from jennifer.march@uni-hamburg.de
*Output:* creates four .RData files which are used for all subsequent behavioural ("had_data.RData") and modeling analyses ("food_modeling_data.RData", "discount_modeling_data.RData" and "social_modeling_data.RData")

**B_GLMMs.Rmd** 
*Input:* had_data.RData
*Output:* GLMM results, Fig 2, Fig S2-S5, had_data2.RData (recoded dwell time information relevant for D_EyeTracking.Rmd)

**C_Graphs.Rmd** 
*Input:* had_data.RData
*Output:* Fig S1

**D_EyeTracking.Rmd**
*Input:* had_data.RData, food_modeling_data.RData, discount_modeling_data.RData and social_modeling_data.RData
*Output:* eyetracking analysis, Fig 3, Fig S6

## 2 Modeling Analyses

**Requirements:** JAGS, dwiener

All analyses are carried out in R using RMarkdown (Version 4.4.1) and JAGS (see Folder BayesModels). The analyses are named after the respective model, except for the preprocessing, which essentially transforms the food_modeling_data.RData into the the files 
*data_prep.RData*, *data_discount_prep.RData* and *data_social_prep.RData* used for modelling analyses across tasks. At the beginning of each file further information about the analyses is provided. 

**Note** we used parallel computing to speed up model fit, yet, this only works for mac and linux, windows uses a requires a different function to do parallel computing. Nevertheless, the analyses also work on windows just without parallel computing and thus takes longer.

**A_Prep_Modeling_###.Rmd** 
*Input:* food_modeling_data.RData/ discount_modeling_data.RData/ social_modeling_data.RData
*Output:* data_prep.RData/ data_discount_prep.RData/ data_social_prep.RData

**B_maaDDM2phi_###.Rmd** 
*Input:* data_prep.RData/ data_discount_prep.RData/ data_social_prep.RData
*Output:* posterior parameter distributions (maaDDM2phisp_###.RData), posterior predictive checks (ppc_maaDDM2phisp_###.RData), Fig S8/9/10

**C_parameter_recovery_2phi_###.Rmd** 
*Input:* data_prep.RData/ data_discount_prep.RData/ data_social_prep.RData and results from respective model fit (maaDDM2phisp_###.RData)
*Output:* parameter recovery (param_recov_2phi_sp_###.RData), Fig S11/12/13

**Model_Plots_HAD.RData.Rmd**
*Input:* maaDDM2phisp_food.RData, maaDDM2phisp_discount.RData, maaDDM2phisp_social.RData
*Output:* Fig 4, Fig S7


### Recommended Folder Structure
*…to speed up reproducibility*

1.	Create folder (e.g. “food_modeling”)
2.	Create R.proj for this folder 
3.	Copy R data file into folder (e.g. “food_modeling_data.RData”)
4.	Copy all modeling scripts (.Rmd files) into folder
5.	Dowload Folder “BayesModels” (.txt files for the modeling in JAGS) and copy into folder, such that “BayesModels” is a folder in the folder “food_modeling”
**Note** If you rename the folder BayesModels, or the models in that folder, you have to rename them when calling the models in the .Rmd files!
6.	Run analyses!
