###########################################################
# Multi-state Cox PH demo for the EBMT3 example dataset
# Script for generating sample data
#
# Author: Stelios Boulitsakis Logothetis
# Email: sb2690@medschl.cam.ac.uk
###########################################################
# 
# This demo script uses the EBMT3 sample dataset. 
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
library(mstate)

# Set seed for reproducibility
set.seed(123)

#################### 1. Data Generation ###################
# Generate some sample data
###########################################################
data(ebmt3)
df_ebmt3 <- ebmt3 %>% 
    mutate(prtime=prtime/365.25, rfstime=rfstime/365.25) %>% 
    tmerge(., ., id, relapse=event(rfstime, rfsstat), recovery=event(prtime, prstat)) %>% 
    tmerge(., ., id, enum=cumtdc(tstart)) %>% 
    mutate(event = ifelse(relapse==1, 2, recovery)) %>% 
    mutate(event = factor(event, 0:2, labels=c("censor", "recovery", "relapse"))) %>% 
    select(id, tstart, tstop, event, age, drmatch, dissub, tcd) 

# Introduce an initial state - which may be "healthy" but may also be "recovery"
df_ebmt3 <- df_ebmt3 %>% 
    mutate(istate = survcheck(Surv(tstart, tstop, event) ~ 1, df_ebmt3, id=id)$istate) %>% 
    mutate(
        modify_flag = event == "relapse" & runif(n()) < 0.3,  # Flag selected rows
        tstart = if_else(modify_flag, 0, tstart)  # Modify `tstart`
    ) %>% 
    # Keep only the flagged `relapse` rows and remove all other rows with matching IDs
    filter(modify_flag | !id %in% id[modify_flag]) %>%
    select(-modify_flag)  # Remove temporary flag column

# Ensure the transitions are valid using survcheck
check <- survcheck(Surv(tstart, tstop, event) ~ 1, data=df_ebmt3, istate=istate, id=id)
stopifnot(all(check$flag == 0))
print(check$transitions)

# Sample 80% of the IDs for the training set
train_ids <- sample(unique(df_ebmt3$id), size = 0.9 * length(unique(df_ebmt3$id)))

# Split dataframe by filtering on ID
df_training <- df_ebmt3[df_ebmt3$id %in% train_ids, ]
df_validation <- df_ebmt3[!df_ebmt3$id %in% train_ids, ] %>% distinct(id, .keep_all = TRUE)

write.csv(df_validation, "ebmt3_validation.csv")
write.csv(df_training, "ebmt3_training.csv")