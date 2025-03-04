# This is a minimal, reproducible example of a bug impacting survfitAJ in the R survival package version 3.7.0
# Relevant github issue: https://github.com/therneau/survival/issues/278
library(survival)
library(ggsurvfit)

data(cancer, package="survival")

# First, we prepare the multi-state dataset.
# This code is copied directly from the survival package's multi-state vignette, section 2.2
ptemp <- with(mgus2, ifelse(ptime==futime & pstat==1, ptime-.1, ptime))
data3 <- tmerge(mgus2, mgus2,  id=id, death=event(futime, death),
                pcm = event(ptemp, pstat))
data3 <- tmerge(data3, data3, id, enum=cumtdc(tstart))
temp <- with(data3, ifelse(death==1, 2, pcm))
data3$event <- factor(temp, 0:2, labels=c("censor", "pcm", "death"))

# Finally, we fit the Cox PH model and run survfit2
msfit <- coxph(Surv(tstart, tstop, event) ~ sex + age, data=data3, id=id)
ndata <- expand.grid(sex=c("F", "M"), age=c(60, 80))
mssurv <- survfit2(msfit, newdata=ndata, p0=c(1,0,0))

# Run ggcuminc
ggcuminc(mssurv, outcome=c("censor", "pcm", "death"))
