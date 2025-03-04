# This script implements a multi-state model using flexsurv, as demonstrated in the vignette: https://cran.r-project.org/web/packages/flexsurv/vignettes/multistate.pdf
library(tidyverse)
library(flexsurv)
library(survival)
library(mstate)

# Jointly fit a RP model:
cfrp <- flexsurvspline(Surv(Tstart, Tstop, status) ~ trans, data=bosms3)
summary(cfrp, type='cumhaz') %>%  str

# Or, with transition-specific models:
model.1 <- flexsurvspline(Surv(Tstart, Tstop, status) ~1, subset=(trans==1), data=bosms3)
model.2 <- flexsurvspline(Surv(Tstart, Tstop, status) ~1, subset=(trans==2), data=bosms3)
model.3 <- flexsurvspline(Surv(Tstart, Tstop, status) ~1, subset=(trans==3), data=bosms3)
tmat <- trans.illdeath(names = c("(s0)", "recovery", "relapse"))
crfs <- fmsm(model.1, model.2, model.3, trans=tmat)

pmatrix.fs(crfs, t=c(5,10), trans=tmat)

tgrid <- seq(0, 14, by = 0.1)
mf <- msfit.flexsurvreg(cfrp, t=tgrid, trans=tmat)

ptf <- probtrans(mf, predt=0, direction="forward")[[1]]

round(ptf[ptf$time %in% c(5, 10),], 3)
