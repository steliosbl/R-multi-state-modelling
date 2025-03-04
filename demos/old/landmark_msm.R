library(tidyverse)
library(magrittr)
library(survival)
library(mstate)


########### A) Study Cohort ###########
# We use the Second Primary Lung Cancer (SPLC) dataset from Fries et al., 2023 (https://doi.org/10.1093/ije/dyad122)
# It contains 3,332 ever-smoking patients from the Multiethnic Cohort Study who were diagnosed with initial primary lung cancer between 1993-2007, 
# and followed-up for SPLC through 2017. Patients without SPLC at 31 Dec 2017 were censored.
# The quantity of interest is time from IPLC diagnosis to SPLC incidence
# Competing events include lung cancer death pre-SPLC incidence and other-cause mortality. 
#######################################

devtools::install_github("thehanlab/dynamicLM")
library(dynamicLM)
data(splc)
splc %<>% rename(id=ID)


# For now, replace all competing events with censoring
# splc$event <- ifelse(splc$event %in% c(2, 3), 0, splc$event) 
# splc$event 

outcome <- list(time="Time", status="event")
fixed_covariates <- c("age.ix", "male", "bmi", "stage.ix",
                      "surgery.ix", "radiation.ix", "chemo.ix", "quityears")
time_varying_covariates <- c("smkstatus", "packyears")
all_covariates <- list(fixed = fixed_covariates, varying = time_varying_covariates)

W <- 5                       # 5-year prediction window
landmarks <- seq(0, 3, by=1) # 0, 1, 2, and 3 years

# We add the column "LM" to indicate the landmark for each row
landmark_data <- stack_data(
  data=splc, outcome=outcome, lms=landmarks, w=W, covs=all_covariates, 
  format="long", rtime="T.fup", id="id")
table(landmark_data$data$LM)

# We update the time-varying covariate of years since quitting smoking by incrementing it with each landmark year
former_smokers <- landmark_data$data$smkstatus == 2
landmark_data$data[former_smokers, "quityears"] <- landmark_data$data[former_smokers, "quityears"] + landmark_data$data[former_smokers, "LM"]

# We add linear and quadratic interaction terms with the landmark year for each of the below covariates
# The new columns are named `column_1`: column x landmark, `column_2`: column x landmark^2
# We additionally include transformations of the landmark year itself: LM_1 (unchanged from LM) and LM_2 (LM^2)
interaction_covariates <- c("packyears", "radiation.ix", "stage.ix")
landmark_data %<>% add_interactions(lm_covs=interaction_covariates, 
                                    func_covars = c("linear", "quadratic"),
                                    func_lms = c("linear", "quadratic"))
colnames(landmark_data$data)

# Expand status and time to prepare for msprep
landmark_data$data %<>%
  mutate(
    s_relapse = as.integer(event == 1), 
    s_death_cancer = as.integer(event == 2),
    s_death_other = as.integer(event == 3),
    t_relapse = Time,
    t_death_cancer = Time,
    t_death_other = Time
  )

tmat <- trans.comprisk(3, c("relapse", "death_cancer", "death_other"))

msbmt <- msprep(
  c(NA, "t_relapse", "t_death_cancer", "t_death_other"), 
  c(NA, "s_relapse", "s_death_cancer", "s_death_other"), 
  landmark_data$data, tmat, keep=c("male", "packyears", "quityears", "radiation.ix", "radiation.ix_1", "stage.ix", "stage.ix_1", "stage.ix_2", "LM_1", "LM_2"))



formula <- as.formula("Surv(Tstart, Tstop, status) ~ male + packyears + quityears + radiation.ix + 
            radiation.ix_1 + stage.ix + stage.ix_1 + stage.ix_2 + LM_1 + LM_2 + cluster(id) + strata(trans)")
c1 <- coxph(formula, msbmt, method="breslow")
print(c1)

formula <- as.formula("Surv(Time, event) ~ male + packyears + quityears + radiation.ix + 
            radiation.ix_1 + stage.ix + stage.ix_1 + stage.ix_2 + LM_1 + LM_2 + cluster(id)")
supermodel <- dynamic_lm(landmark_data, formula, "coxph", x=TRUE)
supermodel <- coxph(formula, landmark_data$data, method="breslow")



print(supermodel)


