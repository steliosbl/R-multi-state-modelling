# This is a minimal, reproducible example of a bug impacting survfitAJ in the R survival package version 3.7.0
# Relevant github issue: https://github.com/therneau/survival/issues/278
library(survival)

data(cancer, package="survival")

# Prepare the multi-state dataset. Code copied directly from the multi-state vignette section 2.2
ptemp <- with(mgus2, ifelse(ptime==futime & pstat==1, ptime-.1, ptime))
data3 <- tmerge(mgus2, mgus2,  id=id, death=event(futime, death),
                  pcm = event(ptemp, pstat))
data3 <- tmerge(data3, data3, id, enum=cumtdc(tstart))
temp <- with(data3, ifelse(death==1, 2, pcm))
data3$event <- factor(temp, 0:2, labels=c("censor", "pcm", "death"))

# Add an istate column
check <- survcheck(Surv(tstart, tstop, event) ~ 1, id=id, data3)
data3$istate <- check$istate

data3[1499, "istate"] <- "pcm"

# Fit Cox PH model and run survfit
msfit <- coxph(Surv(tstart, tstop, event) ~ sex + age, data=data3, id=id, istate=istate)
ndata <- expand.grid(sex=c("F", "M"), age=c(60, 80))
mssurv <- survfit(msfit, newdata=ndata)

# Fix by either specifying start.time (e.g., start.time=0.1) or p0 (e.g., p0=c(1,0,0))
# The latter indicates that everyone in `ndata` is assumed to start in the 1st state (in this case, (s0)). p0=c(0,1,0) would be the 2nd state (pcm), etc.