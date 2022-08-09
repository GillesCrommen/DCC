# DCC
R code for "A Gaussian model for survival data subject to dependent censoring and confounding".

This repository contains all the R code used for the simulations and data application in the paper "A Gaussian model for survival data subject to dependent censoring and confounding" (http://arxiv.org/abs/2208.04184).

The Functions.R file is called in almost all the other files and contains the likelihood functions, asymptotic variance calculations and confidence interval construction among other things.

The JTPA clean data.R file starts from the full JTPA data set that can be found at https://www.upjohn.org/data-tools/employment-research-data-center/national-jtpa-study and cleans it. We extract variables such as race, gender, age, treatment status, whether they have a high school diploma, whether they are married, whether they participated in JTPA programs, the censoring indicator and the unemployment duration. This file compiles data from 2 follow-up surveys. The cleaned data set is included in the repository as clean_dataset_JTPA.csv.

The Simulations.R and JTPA application.R files can be used to recreate the results shown in sections 4 and 5 of the paper respectively. 
