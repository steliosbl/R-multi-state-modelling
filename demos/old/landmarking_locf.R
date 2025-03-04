library(joineR)
library(dplyr)
library(magrittr)
library(survival)
library(ggsurvfit)
library(nlme)
library(zeallot)

#
# Adapted from: https://github.com/survival-lumc/Landmarking2.0
#

standardise_years <- function(x) x * 365 / 365.25

data(liver)

# Dataframe of the survival outcome/censoring for each patient ID (one row per ID)
df_survival <- liver[, c("id", "survival", "cens", "treatment")] %>%
  subset(!duplicated(id))
names(df_survival) = c("id", "survival_years", "status", "treatment")
df_survival$survival_years %<>% standardise_years()
df_survival$treatment %<>% factor(levels=0:1, labels=c("Placebo", "Prednisone"))

# Dataframe of the longitudinal measurements for each patient ID (many rows per ID)
df_longitudinal <- liver[, c("id", "time", "prothrombin", "treatment")] %>%
  subset(time > 0)
names(df_longitudinal) <- c("id", "measurement_time", "prothrombin", "treatment")
df_longitudinal$measurement_time %<>% standardise_years()
df_longitudinal$treatment %<>% factor(levels=0:1, labels=c("Placebo", "Predinsone"))

makeLMdata <- function(df1, df2, col_id, col_time1, col_time2, LM) {
  # Left-truncate all records with event/censoring time prior to the landmark
  df1_lm <- df1[df1[[col_time1]] > LM, ]
  pat_lm <- sort(unique(df1_lm[[col_id]]))
  
  
  # Landmark df_2 (marker data) 
  df2_lm <- df2[df2[[col_id]] %in% pat_lm, ]
  df2_lm_before <- df2_lm[df2_lm[[col_time2]] <= LM, ]
  df2_lm_after <- df2_lm[df2_lm[[col_time2]] > LM, ]
  
  return (list(df1_lm=df1_lm, df2_lm_before=df2_lm_before, df2_lm_after=df2_lm_after))
}

c(df1_lm, df2_lm_before, df2_lm_after) %<-% makeLMdata(df_survival, df_longitudinal, "id", "survival_years", "measurement_time", 3)

makeLOCFdata <- function(df1_lm, df2_lm_before, col_id, col_time1, col_status, col_time2, col_marker, LM, width) {
  # Obtain the most recent ob per patient ID
  locf <- df2_lm_before %>% 
    group_by(id) %>%
    arrange(desc(measurement_time), .by_group=TRUE) %>%
    slice(1)
    
  
  # Join with the survival dataframe
  locf %<>% merge(df1_lm[, c(col_id, col_time1, col_status)], by=col_id, all.x=TRUE)
  
  # Administrative censoring at time LM+width
  cens_mask <- which(locf[[col_time1]] > LM+width)
  locf[cens_mask, col_time1] <- LM+width
  locf[cens_mask, col_status] <- 0
  
  return (locf)
}

df_locf <- makeLOCFdata(df1_lm, df2_lm_before, "id", "survival_years", "status", "measurement_time", "prothrombin", 3, 2)

fitLM <- function(df_locf) {
  cph <- coxph(Surv(survival_years, status) ~ prothrombin, data=df_locf)
  sf <- survfit(cph, newdata=data.frame(prothrombin=df_locf$prothrombin))
  y_pred <- as.vector(unlist(summary(sf, times=LM+width)$surv))
  
  return (list(cph=cph, sf=sf, y_pred=y_pred))
}

c(cph, sf, y_pred) %<-% fitLM(df_locf)


