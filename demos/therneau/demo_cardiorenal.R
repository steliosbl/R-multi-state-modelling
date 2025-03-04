###########################################################
# Multi-state Cox PH demo for the cardio-renal AB model
# Script for model fitting and prediction
#
# Author: Stelios Boulitsakis Logothetis
# Email: sb2690@medschl.cam.ac.uk
###########################################################
#
# This demo script generates a sample dataset for a cardio-renal AB model.
# The state transition data follows this transition matrix:
# to
# from       CVD  CKD Death Both (censored)
# Healthy    1    1     1    0          1
# CVD        0    0     1    1          1 
# CKD        0    0     1    1          1
# Death      0    0     0    0          0
# Both       0    0     0    0          0

######################## 0. Setup #########################
# Load the required libraries
###########################################################
library(tidyverse)
library(survival)

######################## 1. Dataset #######################
# We load the fitting and validation data 
###########################################################
df_training <- read.csv('cardiorenal_training.csv')
df_validation <- read.csv('cardiorenal_validation.csv')

# Convert event columns (FROM and TO) to factors. Censoring must be the first level of the factor.
df_training <- df_training %>% 
    mutate(FROM = factor(FROM, levels=c("Healthy", "CVD", "CKD"))) %>% 
    mutate(TO = factor(TO, levels=c("Censor", "CVD", "CKD", "Death", "Both")))

##################### 2. Model Fitting ####################
# Fit the multi-state Cox model
###########################################################

m.formula <- Surv(TSTART, TSTOP, TO) ~ strata(cov_sex) + cov_age + cov_smoking + 
    cov_hba1c + cov_hdl + cov_sbp + cov_tchol + cov_lnegfr + 
    cov_diabetes + cov_vaccinated

m.cox <- coxph(
    m.formula,
    data=df_training,
    id=PERSON_ID,
    istate=FROM
)

###################### 3. Prediction ######################
# Draw predictions for new samples
###########################################################

# Perform one survfit.coxphms call per initial state sub-group
# Ordering of p0 indices is according to survcheck$states:
# > "Healthy" "CVD"     "CKD"     "Death"   "Both" 
pred.Healthy <- survfit(m.cox, newdata=filter(df_validation, istate=='Healthy'), p0=c(1,0,0,0,0))
pred.CVD <- survfit(m.cox, newdata=filter(df_validation, istate=='CVD'), p0=c(0,1,0,0,0))
pred.CKD <- survfit(m.cox, newdata=filter(df_validation, istate=='CKD'), p0=c(0,0,1,0,0))

########################### End ###########################
