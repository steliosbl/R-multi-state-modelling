# This is a minimal, reproducible example of a bug impacting predict.coxphms in the R survival package version 3.8.3
# Relevant github issue: https://github.com/therneau/survival/issues/302
library(survival)

data(cancer, package="survival")

# Prepare the multi-state dataset. Code copied directly from the multi-state vignette section 2.2
ptemp <- with(mgus2, ifelse(ptime==futime & pstat==1, ptime-.1, ptime))
data3 <- tmerge(mgus2, mgus2,  id=id, death=event(futime, death),
                pcm = event(ptemp, pstat))
data3 <- tmerge(data3, data3, id, enum=cumtdc(tstart))
temp <- with(data3, ifelse(death==1, 2, pcm))
data3$event <- factor(temp, 0:2, labels=c("censor", "pcm", "death"))
data3 <- data3[complete.cases(data3), ]

# Fit Cox PH model and run survfit
msfit <- coxph(Surv(tstart, tstop, event) ~ age + strata(sex), data=data3, id=id, model=TRUE)

# Obtain a small new data frame
ndata <- data3[1:2,]

source('predict.coxphms.R')

predict_mine(msfit, newdata=ndata, type="survival")
