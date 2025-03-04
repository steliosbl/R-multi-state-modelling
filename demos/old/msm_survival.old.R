library(tidyverse)
library(survival)
library(mstate)
library(survivalROC)

data(cancer, package="survival")

df_mgus <- mgus2 %>% 
    # Solve ties in event times
    mutate(ptime := ifelse(ptime==futime & pstat==1, ptime-.1, ptime)) %>%
    # Rename and drop some columns
    rename(time_pcm=ptime, status_pcm=pstat, status_death=death, time_death=futime) %>% 
    select(id, age, sex, hgb, creat, mspike, time_pcm, status_pcm, time_death, status_death)
head(df_mgus)


# Data Format 1 - classic / counting process
# id, tstart, tstop, event, covariates
df_counting <- tmerge(
    df_mgus, df_mgus, id=id, 
    death=event(time_death, status_death), 
    pcm=event(time_pcm, status_pcm)
)
df_counting <- tmerge(df_counting, df_counting, id, enum=cumtdc(tstart)) %>% 
    mutate(event = ifelse(death==1, 2, pcm)) %>% 
    mutate(event = factor(event, 0:2, labels=c("censor", "pcm", "death"))) 
tail(df_counting)

mfit1 <- survfit(Surv(tstart, tstop, event) ~ sex, data=df_counting, id=id)
check1 <- survcheck(Surv(tstart, tstop, event)~ 1, id=id, df_counting)
stopifnot(all(check1$flag == 0))
print(mfit1, rmean=240, digits=2)

# Data Format 2 - new / timeline 
# id, time, event & a row at time=0 for covariates
add_indicator <- function(df, col_time, col_status, value) {
    df_ind <- df %>% 
        filter(.[,col_status] == 1) %>% 
        mutate(time=.[,col_time], event=value) %>%
        select(id, time, event)
    return(merge(df, df_ind, all=TRUE))
}

df_timeline <- df_mgus %>% 
    mutate(time=0, event=NA) %>%
    mutate(time_censor = time_death, status_censor = ifelse(status_death == 0, 1, 0)) %>% 
    add_indicator('time_censor', 'status_censor', 0) %>% 
    add_indicator('time_pcm', 'status_pcm', 1) %>% 
    add_indicator('time_death', 'status_death', 2) %>%
    mutate(event=factor(event, 0:2, c("censor", "pcm", "death"))) %>% 
    select(id, age, sex, time, event)
head(df_timeline)

mfit2 <- survfit(Surv2(time, event) ~ sex, data=df_timeline, id=id)
check2 <- survcheck(Surv2(time, event) ~ 1, id=id, df_timeline)
stopifnot(all(check2$flag == 0))
print(mfit2, rmean=240, digits=2)

# Data Format 1.1 classic but with msprep
tmat <- transMat(list(c(2,3,4), c(3,4), c(), c()), names=c("entry", "pcm", "death", "censor"))
tmat 

df_ms_tmp <- df_mgus %>%
    mutate(time_censor = time_death, status_censor = ifelse(status_death == 0, 1, 0))

msbmt <- msprep(
    data=df_ms_tmp,
    time=c(NA, 'time_pcm', 'time_death', 'time_censor'),
    status=c(NA, 'status_pcm', 'status_death', 'status_censor'),
    trans=tmat,
    id='id',
    keep=c('age', 'sex')
) %>% 
    rename(tstart=Tstart, tstop=Tstop) %>% 
    mutate(event = factor(to, c(4,2,3), labels=c("censor", "pcm", "death"))) %>% 
    mutate(istate = factor(from, 1:4, labels=c('s0', 'pcm', 'death', 'censor'))) %>%
    filter(status==1) %>% 
    select(id, age, sex, tstart, tstop, event, istate)
mfit3 <- survfit(Surv(tstart, tstop, event) ~ sex, data=msbmt, id=id)
check3 <- survcheck(Surv(tstart, tstop, event)~ 1, id=id, msbmt, istate=istate)
check3$flag
print(mfit3, rmean=240, digits=2)


# Ensure all datasets give the same stuff 
all.equal(check1$transitions, check2$transitions)
all.equal(check1$transitions, check3$transitions)
all.equal(mfit1$pstate, mfit2$pstate)
all.equal(mfit1$pstate, mfit3$pstate)

# Cox PH models
model.cox.1 <- coxph(Surv(tstart, tstop, event) ~ sex, data=df_counting, id=id)
model.cox.2 <- coxph(Surv2(time, event) ~ sex, data=df_timeline, id=id)
model.cox.3 <- coxph(Surv(tstart, tstop, event) ~ sex, data=msbmt, id=id, istate=istate)
# Verify results
all.equal(coef(model.cox.1), coef(model.cox.2))
all.equal(coef(model.cox.1), coef(model.cox.3))

# Inference
df_ndata <- expand.grid(sex=c("F", "M"), age=c(60, 80))

inf.csurv.1 <- survfit(model.cox.1, newdata=df_ndata)
inf.csurv.2 <- survfit(model.cox.2, newdata=df_ndata)
inf.csurv.3 <- survfit(model.cox.3, newdata=df_ndata)

plot(inf.csurv.3, col=1:2, xscale=365.25, xlab='Years', ylab='Survival')

# More inference

model.cox.simple <- coxph(Surv(time_death, status_death) ~ age+sex+age*mspike+hgb+age*hgb+age*age, id=id, data=df_mgus)
df_newdata <- df_mgus %>% filter(!is.na(mspike) & !is.na(hgb)) %>% sample_n(365)
inf.simple <- survfit(model.cox.simple, newdata=df_newdata)

y_pred_proba.s <- inf.simple$surv %>% .[nrow(.),]
plot(y_pred_proba.s)
c_stat <- Hmisc::rcorr.cens(y_pred_proba.s, Surv(df_newdata$time_death, df_newdata$status_death))
c_stat['C Index']

model.cox.multi <- coxph(Surv(tstart, tstop, event) ~ age+sex+age*mspike+hgb+age*hgb+age*age, id=id, data=df_counting)
inf.multi <- survfit(model.cox.multi, newdata=df_newdata)
y_pred_proba.m <- 1-inf.multi[,'death']$pstate[283,,1]
plot(y_pred_proba.m)
c_stat <- Hmisc::rcorr.cens(y_pred_proba.m, Surv(df_newdata$time_death, df_newdata$status_death))
c_stat['C Index']

