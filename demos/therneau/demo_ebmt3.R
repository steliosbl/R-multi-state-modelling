###########################################################
# Multi-state Cox PH demo for the EBMT3 example dataset
# Script for model fitting and prediction
#
# Author: Stelios Boulitsakis Logothetis
# Email: sb2690@medschl.cam.ac.uk
###########################################################
# 
# This demo scripts uses the EBMT3 sample dataset. 
# The state transition data follows this transition matrix:
# to
# from       recovery relapse (censored)
# (s0)            1       1          1
# recovery        0       1          1
# relapse         0       0          0

######################## 0. Setup #########################
# Load the required libraries
###########################################################
library(tidyverse)
library(survival)

######################## 1. Dataset #######################
# We load the fitting and validation data 
###########################################################
df_training <- read.csv('ebmt3_training.csv')
df_validation <- read.csv('ebmt3_validation.csv')

# Convert event columns (event and istate) to factors. Censoring must be the first level of the factor.
df_training <- df_training %>%  
    mutate(event = factor(event, labels=c("censor", "recovery", "relapse"))) %>% 
    mutate(istate = factor(istate, labels=c("(s0)", "recovery")))
    

##################### 2. Model Fitting ####################
# Fit the multi-state Cox model
###########################################################
m.cox <- coxph(
    Surv(tstart, tstop, event, type='mstate') ~ age + drmatch + dissub + tcd,
    data=df_training,
    id=id,
    istate=istate
)

###################### 3. Prediction ######################
# Draw predictions for new samples
###########################################################

# Perform one survfit.coxphms call per initial state sub-group
# Ordering of p0 indices is according to survcheck$states:
# > "(s0)"     "recovery" "relapse" 
pred.s0 <- survfit(m.cox, newdata=filter(df_validation, istate=='(s0)'), p0=c(1,0,0))
pred.recovery <- survfit(m.cox, newdata=filter(df_validation, istate=='recovery'), p0=c(0,1,0))

########################### End ###########################
