The main scripts from CCT are listed below:
- Behaviour/Meta_modelfree.Rmd
- Modeling/Cloth/aDDM.Rmd (for clothing data)
- Modeling/Cloth/aDDM_food.Rmd (for food data)
- Modeling/Cloth/aDDM_social.Rmd (for social data)
- Modeling/Cloth/aDDM_TD.Rmd (for temporal discounting data)

Everything organized here so far (Nov. 15. 2024.) is slightly different from Jenna's original setting. That is, the choice is encoded as left (1) and right (0) rather than choosing option 1 with attribute A1>attribute A2 (1) or with attribute A1<attribute A2 (0).

Note: The scripts used the data from Jenna (food, social, temporal discounting) and Yang et al., 2022 (clothing). Since the script relies on the data organized by other scripts (e.g., Behaviour/A_Preprocess.Rmd), WE HAVE TO RERUN THIS SCRIPT if any change in those files is made. 
