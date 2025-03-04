###########################################################
# Multi-state Cox PH demo for the cardio-renal AB model
# Script for generating sample data
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

# Set seed for reproducibility
set.seed(123)

#################### 1. Data Generation ###################
# Generate some sample data
###########################################################

# Number of individuals
n_individuals <- 10000

# Create a dataframe with dummy covariate data
generate_predictors <- function(n_individuals) {
    data.frame(
        PERSON_ID = paste0("ID_", 1:n_individuals),
        cov_age = (runif(n_individuals, 30, 80)-60)/5,
        cov_sex = sample(c(0, 1), n_individuals, replace = TRUE), # 0 = Male, 1 = Female
        cov_smoking = sample(c(0, 1), n_individuals, replace = TRUE), # 0 = Non-smoker, 1 = Smoker
        cov_hba1c = (runif(n_individuals, 20, 150) - 31)/9.34, # Typical HbA1c range
        cov_hdl = (runif(n_individuals, 0.5, 2.5) - 1.3)/0.5, # HDL cholesterol range
        cov_sbp = (runif(n_individuals, 90, 180) - 120)/20, # Systolic blood pressure range
        cov_tchol = (runif(n_individuals, 3, 7) - 6), # Total cholesterol range
        cov_lnegfr = (log(runif(n_individuals, 30, 120)) - 4.5)/0.15,
        cov_diabetes = sample(c(0, 1), n_individuals, replace = TRUE), 
        cov_vaccinated = sample(c(0, 1), n_individuals, replace = TRUE)
    )
}
df_predictors <- generate_predictors(n_individuals)

# Generate transitions for each individual
generate_transitions <- function(id, time_max=300) {
    # Randomly assign a starting state
    states <- c(sample(c("Healthy", "CVD", "CKD"), 1, prob = c(0.6, 0.2, 0.2)))
    times <- c(0)
    
    # Define possible transitions
    while (!(tail(states, 1) %in% c("Death", "Both", "Censor"))) {
        last_state <- tail(states, 1)
        
        # Define next possible states based on the last state
        if (last_state == "Healthy") {
            next_state <- sample(c("CVD", "CKD", "Death", "Censor"), 1, prob = c(0.2, 0.2, 0.2, 0.4))
        } else if (last_state == "CVD") {
            next_state <- sample(c("Death", "Both", "Censor"), 1, prob = c(0.3, 0.2, 0.5))
        } else if (last_state == "CKD") {
            next_state <- sample(c("Death", "Both", "Censor"), 1, prob = c(0.1, 0.1, 0.8))
        } else {
            break
        }
        
        # Generate the next time point ensuring we do not exceed time_max
        next_time <- tail(times, 1) + sample(1:200, 1)
        
        if (next_state == "Censor" || next_time >= time_max) {
            next_time <- time_max
            next_state <- "Censor"
        }
        
        # Append state and time
        states <- c(states, next_state)
        times <- c(times, next_time)
        
        # Stop further transitions if a terminal state is reached
        if (next_state == "Censor" || next_state == "Death" || next_state == "Both") {
            break
        }
    }
    
    # Create dataframe for this individual
    data.frame(
        PERSON_ID = paste0("ID_",id),
        TSTART = head(times, -1),
        TSTOP = tail(times, -1),
        FROM = head(states, -1),
        TO = tail(states, -1)
    )
}

# Generate transitions for all individuals
df_states <- do.call(rbind, lapply(1:n_individuals, generate_transitions))

# Convert event columns (FROM and TO) to factors. Censoring must be the first level of the factor.
df_states <- df_states %>% 
    mutate(FROM = factor(FROM, levels=c("Healthy", "CVD", "CKD"))) %>% 
    mutate(TO = factor(TO, levels=c("Censor", "CVD", "CKD", "Death", "Both")))

# Check the validity of the state transitions with survcheck
check <- survcheck(Surv(TSTART, TSTOP, TO) ~ 1, data=df_states, istate=FROM, id=PERSON_ID)
stopifnot(all(check$flag == 0))
print(check$transitions)

# Create the final fitting dataset
df_training <- left_join(df_states, df_predictors, by="PERSON_ID")

# Generate some validation samples
df_validation <- generate_predictors(100) %>% 
    mutate(istate = c(sample(c("Healthy", "CVD", "CKD"), 100, replace=T, prob = c(0.6, 0.2, 0.2))))

write.csv(df_training, "cardiorenal_training.csv")
write.csv(df_validation, "cardiorenal_validation.csv")