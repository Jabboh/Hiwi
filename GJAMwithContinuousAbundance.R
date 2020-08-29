################################Analysis of Continuous Abundances of Several Mosquito Species###
#Running a JSDM with annual continous abundances with the model GJAM:
#1. Data preparation
#2. Fitting the most complex model (all environmental covariates 
#+ according interaction and quadratic terms)
#3. Internal Model validation of this models by analyzing its residuals with DHARMa
#4. Selecting the "best" model by step-wise selection based on the DIC 
#>> Finding the "most appropriate" specifications of covariates
#5. Redoing the internal validation for the "most appropriate" model.
#6. Results
#6.1. In-sample (not on test data) or maybe model exploration?:
#a. Comparing the coefficients: Size and credibility intervals
#b. Correlations between the responses / Residual Correlation Parameter of gjam
#c. Response Curves 
#d. Variable importance
#6.2. Out-of-Sample
#a. Conditional Predictions of gjam vs. Unconditional Predictions of 
#gjam 
#b. Comparing the uncertainty of the different "prediction"-types
rm(list=ls())
setwd("C:\\Users\\jakob\\Documents\\Hiwi\\Data")
#install.packages("gridExtra")

#loading packages
library(readxl) #to read in the data
library(rstanarm) #Doing Bayesian probit GLMs for single species
library(gjam) #Doing the joint estimation with the GJAM modelling aproach
library(ggplot2) #for plotting
library(DHARMa) # for checking in-sample validity of the models
library(dplyr) # for simplified syntax and neater code
