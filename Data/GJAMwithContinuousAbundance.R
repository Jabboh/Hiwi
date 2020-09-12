################################Analysis of Continuous Abundances of Several Mosquito Species###
#Running a JSDM with annual continous abundances with the model GJAM:
#1. Data preparation
#2. Fitting the most complex model (all environmental covariates 
#+ according interaction and quadratic terms)
#3. Internal Model checking of this model by analyzing its residuals with DHARMa
#4. Selecting the "best" model by step-wise selection based on DIC 
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
#install.packages("mgsub")
#loading packages
library(readxl) #to read in the data
library(rstanarm) #Doing Bayesian probit GLMs for single species
library(gjam) #Doing the joint estimation with the GJAM modelling aproach
library(ggplot2) #for plotting
library(DHARMa) # for checking in-sample validity of the models
library(dplyr) # for simplified syntax and neater code
library(here) #for creating paths and opening/storing files in R
library(gjam) #for modelling with gjam
library(tidyr)
##################################### 1.Data Preperation #######################################
#read in the data (Monthly species PA data for seven different mosquito species
#and according environmental covariates)
df <- read_excel(here("Data/RoizDaten.xlsx"))

#Checking the rough structure of the data set
str(df)

#delete last row, bc it does not contain raw data
df <- df[-nrow(df),]

#set NAs to 0, because these correspond to captures of 0 of the mosquito species at that trap
df[is.na(df)] <- 0

#Reading in the spatial coordinates of the different trap locations
coords <- read_excel(here("Data/Traps_coordenadas_geograficas.xls"))
#adding lon-lat column to the data frame df
df <- merge(df, coords[, c("trap", "Norte", "Oeste")], by="trap", all.x= T, sort = F)

#Select our three response variables (a) Culex modestus, (b) Culex perexiguus and 
#(c) Culex pipiens
y <- as_tibble(df[,c("Cx.pipiens", "Cx.modestus", "Cx.perexiguus")])

#Normalizing continuous Covariates: Otherwise, interaction terms hardly interpretable and skewed
df[,15:18] <- scale(df[,15:18])
df[,24:31] <- scale(df[,15:18])

#make landuse variables to factors
df$Unit <- as.factor(df$Unit)
levels(df$Unit) <- c("Rice", "Cultives", "Marshland", "Scrubland", "SandDunes", "Fishponds")
#Split the data set in training (70 %) and test (30%) set
train_size <- floor(0.7 * nrow(y))
#set a random seed for replicability
set.seed(333)
#sample the training IDs
train_id <- sample(seq_len(nrow(y)), size = train_size)
#partition data into train and test set
train <- df[train_id, ]
test <- df[-train_id, ]
y_train <- y[train_id, ]
y_test <- y[-train_id, ]

##########################2. Fitting the most complex model########################################
#The most complex model includes the following covariates as well as the interactions and
#quadratic terms for the continuous variables:
#1. Unit Code
#2. Hydroperiod(500 m buffer): I chose that buffer, bc Roiz also uses this buffer for negatative
#binomial models of Culex modestus and pipiens
#3. NDVI (2000 m buffer): Roiz uses this buffer in his negative-binomial model for Culex perexiguus
#4. Distance to next rice paddy
#5. Distance to next human settlement
#These covariates yield the following model formula
form <- ~ Unit + (HIDRO_MEAN_500 + NDVI_2000 + DIST_ARROZ + DIST_urbano)^2 +
  I(HIDRO_MEAN_500^2) + I(NDVI_2000^2) + I(DIST_ARROZ^2) + I(DIST_urbano^2)

#Make a dataframe of the covariates in the form gjam requires (separate dummy columns per factor
#for factor variables) and separate columns for interaction and quadratic terms (gjamPredict
#requires this)
model_xdata <- stats::model.matrix(form, train)    
x_train <- as.data.frame(model_xdata)
#make the model formula in the form gjam requires
form_gj <- paste(names(x_train)[-1], collapse ="+")
form_gj <- paste0("~", form_gj)

#Define gjam model settings
#variable type: continuous abundance, because we have mean abundance values (Abundance per trap 
#night)
types <- c("CA", "CA", "CA")
#define model/algorithm parameters: 4000 gibbs steps + burnin of 1000
ml   <- list(ng = 4000, burnin = 1000, typeNames = types)
#run the gjam model
m_com <- gjam(form_gj,
              ydata = y_train, xdata = x_train, modelList = ml)
summary(m_com)
#The model performs very bad (RMSPE = 833), BUT I just want to give an example how a comprehensive
#analysis of mosquito abundance data could look like in gjam.

####3. Internal Model Checking of the model by analyzing its residuals with DHARMa####
##We check the unconditional predictions of our model and throw predictions of all species into
#one pot >> This saves us a lot of work; for a more rigourous analysis, different prediction
#types and different species can be tested individually.
#To get a solid methadological background see: https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
#Or my summary of the method "XXX"

#Create a function that produces a DHARMa-object needed for the analysis
dharm_obj <- function(model, xdata, ydata){
  stopifnot(nrow(xdata) == nrow(ydata))
  #We need simulations (in-sample predictions)
  #define the modelling setting
  newdata <- list(xdata = as.data.frame(xdata), nsim = 4000)
  #calculating the in-sample predictions (simulations)
  sim <- gjamPredict(output = model, newdata = newdata)
  
  #Create a dharma object. For this, I specify the following arguments:
  #1. 4000 simulations of fitted responses per observation
  #(simulated response),(The single draws from the predictive posterior are stored in ychains.)
  # 2. the observed responses (observedResponse),
  # 3. the mean of the 4000 simulated fitted responses which is the "expected" 
  #value of the predicted y of an observation (fittedPredictedResponse; stored in 
  #gjamObject$sdList$yMu),
  #and other arguments concerning the specific handling of the "scaled residuals".
  dharm_gj_un <- createDHARMa(simulatedResponse = t(sim$ychains),
                              observedResponse = unlist(ydata),
                              fittedPredictedResponse = c(sim$sdList$yMu),
                              integerResponse = F, seed = 3333,method = "PIT")
  return(dharm_gj_un)
}

#apply the function
dharm_com <- dharm_obj(m_com, x_train, y_train)


#Plot residual diagnostics
plot(dharm_com)

#1. QQ-plot looks very concerning >> scaled residuals dont follow the expected uniform 
#distribution.
#2. KS-test is significant >> doesn't supports our assumption that the scaled residuals are 
#uniformly distributed >> indicates that our model is incorrectly specified.
#3. Residual vs. predicted quantile deviations is also significant >> another indication that our
#model has problems

#Plotting the residuals against all covariates to check whether we specified the
# functional relationships correctly.
plotResiduals(dharm_com, rep(train$Unit,3)) #looks okish, but unit 6 looks problematic
plotResiduals(dharm_com, rep(x_train$HIDRO_MEAN_500, 3))#significant problem with the .25 and 
#.75 quantiles. Combined quantile test also significant
plotResiduals(dharm_com, rep(x_train$NDVI_2000, 3))#looks bad >> many significant problems
plotResiduals(dharm_com, rep(x_train$DIST_ARROZ, 3))#looks bad >> many significant problems
plotResiduals(dharm_com, rep(x_train$DIST_urbano, 3)) #looks bad >> many significant problems
#the quadratic forms
plotResiduals(dharm_com, rep(x_train$`I(HIDRO_MEAN_500^2)`, 3)) #looks bad >> many significant problems
plotResiduals(dharm_com, rep(x_train$`I(NDVI_2000^2)`, 3)) #looks bad >> many significant problems
plotResiduals(dharm_com, rep(x_train$`I(DIST_ARROZ^2)`, 3)) #looks bad >> many significant problems
plotResiduals(dharm_com, rep(x_train$`I(DIST_urbano^2)`, 3)) #looks bad >> many significant problems
#the interactions
plotResiduals(dharm_com, rep(x_train$`HIDRO_MEAN_500:NDVI_2000`, 3))#looks bad >> many significant problems
plotResiduals(dharm_com, rep(x_train$`HIDRO_MEAN_500:DIST_ARROZ`, 3)) #looks bad >> many significant problems
plotResiduals(dharm_com, rep(x_train$`HIDRO_MEAN_500:DIST_urbano`, 3)) #looks bad >> many significant problems
plotResiduals(dharm_com, rep(x_train$`NDVI_2000:DIST_ARROZ`, 3)) #looks bad >> many significant problems
plotResiduals(dharm_com, rep(x_train$`NDVI_2000:DIST_urbano`, 3)) #looks bad >> many significant problems
plotResiduals(dharm_com, rep(x_train$`DIST_ARROZ:DIST_urbano`, 3))# #looks bad >> many significant problems
#Overall, there are many significant plots >> not for one covariate DHARMa doesn't flag no issues

hist(dharm_com)
#also, doesn't look uniform


# Test for spatial Autocorrelation
#Recalculate residuals per location (group variable)
dharm_com_spatial <- recalculateResiduals(dharm_com,
                                            group = rep(train$Norte, 3))
testSpatialAutocorrelation(dharm_com_spatial, 
                           x =  train$Oeste, 
                           y = train$Norte)
#Moran's I is insignificant >> No spatial autocorrelation
####4. Selecting the "best" model by step-wise selection based on DIC ####
#I proceed as follows. I drop each environmental covariate(Unit, NDVI, Hidromean...) 
#and the related terms (quadratic + interactions) separately and choose the model with
#the lowest DIC (the best reduced model). If this DIC is lower (meaning better)
#than the DIC of the complex model, I repeat the steps with the best reduced model as the new 
#complex model until the more complex model has a lower DIC (meaning better model).
#At the end of this first step,
#I have chosen the principal covariates of my model. Next, I want to choose the
#functional relationship of them with regard to the response. For this, I first
#check whether dropping the interaction terms reduces the DIC. So I drop them
#iteratevly and again choose the model with the lowest DIC. Afterwards I do 
#the same for the quadratic forms of the covariates.

#Note: The DIC should only be applied for model selection if the model is correctly specified
#>> This isn't the case here...
#Moreover, the DIC is not appropriate for all model averaging situations.
#("A guide to Bayesian model selection for ecologists", Hooten & Hobbs, 2015)

#Loading the step functions I created for dropping covariates.
#First argument: fitted model for which covariates should be dropped 
#Second argument: names of the covariates 