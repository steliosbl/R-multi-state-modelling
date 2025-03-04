library(tidyverse)
library(survival)
library(mstate)
library(flexsurv)

library(glue)
library(ggplot2)
library(cowplot)

set.seed(123)  # For reproducibility

################################## 1. Dataset ##################################
data(ebmt3)
data.covs <- c("dissub", "age", "drmatch", "tcd")
data.tmat <- trans.illdeath(names = c("(s0)", "recovery", "relapse"))
data.mstate <- ebmt3 %>%  
    mutate(prtime=prtime/365.25, rfstime=rfstime/365.25) %>% 
    msprep(
        time=c(NA, "prtime", "rfstime"), status=c(NA, "prstat", "rfsstat"),
        keep=data.covs, id="id",
        trans=data.tmat
) %>% expand.covs(data.covs, append=TRUE, longnames=TRUE)
data.tgrid <- unique(pred.mstate$time) #TODO CHANGE THIS

data.survival <- ebmt3 %>%  
    mutate(prtime=prtime/365.25, rfstime=rfstime/365.25) %>% 
    tmerge(., ., id, relapse=event(rfstime, rfsstat), recovery=event(prtime, prstat)) %>% 
    tmerge(., ., id, enum=cumtdc(tstart)) %>% 
    mutate(event = ifelse(relapse==1, 2, recovery)) %>% 
    mutate(event = factor(event, 0:2, labels=c("censor", "recovery", "relapse"))) %>% 
    select(id, tstart, tstop, event, age, drmatch, dissub, tcd) %>% 
    rename(Tstart=tstart, Tstop=tstop)

data.spline <- data.mstate %>% 
    mutate(trans = factor(trans)) %>% 
    filter(Tstop >= 0.01 & Tstop <= 7.11)

list.datasets <- list(
    "mstate_cox" = data.mstate,
    "survival_cox" = data.survival,
    "flexsurv_spline" = data.spline
)

################################### 2. Models ##################################
formula.mstate <- as.formula(glue("Surv(Tstart, Tstop, status) ~ {paste(tail(colnames(data.mstate), 18), collapse=' + ')} + cluster(id) + strata(trans)"))
formula.survival <- Surv(Tstart, Tstop, event, type='mstate') ~ age + drmatch + dissub + tcd
formula.flexsurv <- Surv(Tstart, Tstop, status) ~ trans + age + drmatch + dissub + tcd

list.fitters <- list(
    "mstate_cox" = function(data) coxph(formula.mstate, data=data, method='breslow'),
    "survival_cox" = function(data) coxph(formula.survival, data=data, id=id, model=T),
    "flexsurv_spline" = function(data) flexsurvspline(formula.flexsurv, data=data, k=6)
)

list.fitted <- list() 
list.times <- list()
for (name in names(list.fitters)) {
    fit <- list.fitters[[name]]
    data <- list.datasets[[name]]
    t.start <- Sys.time()
    fitted <- fit(data)
    t.stop <- Sys.time()
    list.fitted[[name]] <- fitted 
    print(t.stop - t.start)
}

################################# 2. Prediction ################################
list.predictors <- list(
    "mstate_cox" = function(ndata) msfit(list.fitted[['mstate_cox']], newdata=ndata %>% rename(strata=trans), trans=data.tmat) %>% probtrans(predt=0),
    "survival_cox" = function(ndata) survfit(list.fitted[['survival_cox']], newdata=ndata),
    "flexsurv_spline" = function(ndata) msfit.flexsurvreg(list.fitted[['flexsurv_spline']], newdata=ndata, t=data.tgrid, trans=data.tmat)
)

list.predictions <- list()
list.times <- list()
for (name in names(list.predictors)) {
    predictor <- list.predictors[[name]]
    ndata <- list.datasets[[name]]
    t.start <- Sys.time()
    predictions <- predictor(ndata)
    t.stop <- Sys.time()
    list.predictions[[name]] <- predictions 
    print(t.stop - t.start)
}

