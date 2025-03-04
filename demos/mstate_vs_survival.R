# This script tests out the function survival:::predict.coxphms and whether it produces valid state occupancy probabilities.
# Tests have been conducted using survival ver. 3.8-3
library(tidyverse)
library(survival)
library(mstate)
library(flexsurv)
library(glue)
library(ggplot2)
library(cowplot)

data(ebmt3)
covs <- c("dissub", "age", "drmatch", "tcd")
data.raw <- ebmt3 %>%
    mutate(prtime = prtime / 365.25, rfstime = rfstime / 365.25)
data.newd.ids <- c(1)

################ 1. Multi-State Cox (mstate) ##############
# Prepare data
model.mstate.tmat <- trans.illdeath(names = c("(s0)", "recovery", "relapse"))
data.mstate <- data.raw %>%
    msprep(
        time = c(NA, "prtime", "rfstime"),
        status = c(NA, "prstat", "rfsstat"),
        keep = covs,
        id = "id",
        trans = model.mstate.tmat
    ) %>%
    expand.covs(covs, append = TRUE, longnames = TRUE)

# Fit model
model.mstate.formula <- as.formula(glue("Surv(Tstart, Tstop, status) ~ {paste(tail(colnames(data.mstate), 18), collapse=' + ')} + cluster(id) + strata(trans)"))
model.mstate.cox <- coxph(model.mstate.formula, data.mstate, method = "breslow")

# Draw predictions
newd.mstate <- data.mstate %>%
    filter(id %in% data.newd.ids) %>%
    select(-c(from, to, Tstart, Tstop, time, status)) %>%
    rename(strata = trans)

pred.mstate <- msfit(
    model.mstate.cox,
    newdata = newd.mstate,
    trans = model.mstate.tmat
) %>%
    probtrans(predt = 0) %>%
    .[[1]] %>%
    head(-1) %>%
    select(time, "(s0)" = pstate1, recovery = pstate2, relapse = pstate3) %>%
    pivot_longer(cols = -time, names_to = "state", values_to = "pstate") %>%
    data.frame()

############## 2. Multi-State Cox (survival) ##############
# Prepare data
data.survival <- data.raw %>%
    tmerge(., ., id, relapse = event(rfstime, rfsstat), recovery = event(prtime, prstat)) %>%
    tmerge(., ., id, enum = cumtdc(tstart)) %>%
    mutate(event = ifelse(relapse == 1, 2, recovery)) %>%
    mutate(event = factor(event, 0:2, labels = c("censor", "recovery", "relapse"))) %>%
    select(id, tstart, tstop, event, age, drmatch, dissub, tcd) %>%
    rename(Tstart = tstart, Tstop = tstop)

# Fit model
model.survival.formula <- Surv(Tstart, Tstop, event, type = "mstate") ~ age + drmatch + dissub + tcd
model.survival.cox <- coxph(model.survival.formula, data.survival, id = id, model = T)

# Verify that the mstate and survival are identical
stopifnot(all.equal(sort(unname(coef(model.survival.cox))), sort(unname(coef(model.mstate.cox)))))

# Draw predictions
newd.survival <- data.raw %>%
    filter(id %in% 1:10) %>%
    select(id, all_of(covs))
pred.survival <- survfit(
    model.survival.cox,
    newdata = newd.survival,
) %>%
    summary(data.frame = T) %>%
    select(time, state, pstate) %>%
    distinct(., .keep_all = T) %>%
    rbind(data.frame(
        time = rep(0.0, 3),
        state = c("(s0)", "recovery", "relapse"),
        pstate = c(1.0, 0.0, 0.0)
    ))

############ 3. Joint Weibull Model (flexsurv) ############
# Prepare data (use the mstate one)
data.weibull <- data.mstate %>%
    mutate(trans = factor(trans))

# Fit model
model.weibull.j.formula <- as.formula(glue("Surv(Tstart, Tstop, status) ~ {paste(tail(colnames(data.weibull), 18), collapse=' + ')} + trans + shape(trans)"))
model.weibull.j <- flexsurvreg(model.weibull.j.formula, data = data.weibull, dist = "weibull")

# Draw predictions
model.weibull.tgrid <- unique(pred.mstate$time)
newd.weibull <- newd.mstate %>%
    distinct(id, .keep_all = T)

pred.weibull.j <- msfit.flexsurvreg(
    model.weibull.j,
    newdata = newd.weibull,
    t = model.weibull.tgrid,
    trans = model.mstate.tmat
) %>%
    probtrans(predt = 0, direction = "forward") %>%
    .[[1]] %>%
    head(-1) %>%
    select(time, "(s0)" = pstate1, recovery = pstate2, relapse = pstate3) %>%
    pivot_longer(cols = -time, names_to = "state", values_to = "pstate") %>%
    data.frame()

##### 4. Transition-Specific Weibull Model (flexsurv) #####
# Fit models and combine
model.weibull.t.formula <- Surv(Tstart, Tstop, status) ~ age + drmatch + dissub + tcd
model.weibull.t <- fmsm(
    flexsurvreg(model.weibull.t.formula, subset = (trans == 1), data = data.weibull, dist = "weibull"),
    flexsurvreg(model.weibull.t.formula, subset = (trans == 2), data = data.weibull, dist = "weibull"),
    flexsurvreg(model.weibull.t.formula, subset = (trans == 3), data = data.weibull, dist = "weibull"),
    trans = model.mstate.tmat
)

# Draw predictions
pred.weibull.t <- msfit.flexsurvreg(
    model.weibull.t,
    newdata = newd.weibull,
    t = model.weibull.tgrid,
    trans = model.mstate.tmat
) %>%
    probtrans(predt = 0, direction = "forward") %>%
    .[[1]] %>%
    head(-1) %>%
    select(time, "(s0)" = pstate1, recovery = pstate2, relapse = pstate3) %>%
    pivot_longer(cols = -time, names_to = "state", values_to = "pstate") %>%
    data.frame()

############# 5. Joint Spline Model (flexsurv) ############
# Prepare data (use the mstate one)
data.spline <- data.mstate %>%
    mutate(trans = factor(trans)) %>%
    filter(Tstop >= 0.01 & Tstop <= 7.11)

# Fit model
model.spline.j.formula <- as.formula(glue("Surv(Tstart, Tstop, status) ~ trans + {paste(tail(colnames(data.spline), 18), collapse=' + ')}"))
model.spline.j <- flexsurvspline(Surv(Tstart, Tstop, status) ~ trans + age + drmatch + dissub + tcd, data = data.spline, k = 6)

# Draw predictions
model.spline.tgrid <- unique(pred.mstate$time)
newd.spline <- newd.mstate %>%
    distinct(id, .keep_all = T)

pred.spline.j <- msfit.flexsurvreg(
    model.spline.j,
    newdata = newd.spline,
    t = model.spline.tgrid,
    trans = model.mstate.tmat
) %>%
    probtrans(predt = 0, direction = "forward") %>%
    .[[1]] %>%
    head(-1) %>%
    select(time, "(s0)" = pstate1, recovery = pstate2, relapse = pstate3) %>%
    pivot_longer(cols = -time, names_to = "state", values_to = "pstate") %>%
    data.frame()

###### 6. Transition-Specific Spline Model (flexsurv) #####
# Fit model
model.spline.t.formula <- Surv(Tstart, Tstop, status) ~ age + drmatch + dissub + tcd
model.spline.t <- fmsm(
    flexsurvspline(model.spline.t.formula, subset = (trans == 1), data = data.spline, k = 6),
    flexsurvspline(model.spline.t.formula, subset = (trans == 2), data = data.spline, k = 6),
    flexsurvspline(model.spline.t.formula, subset = (trans == 3), data = data.spline, k = 6),
    trans = model.mstate.tmat
)

# Draw predictions
pred.spline.t <- pmatrix.fs(
    model.spline.t,
    newdata = newd.spline,
    t = model.spline.tgrid,
    trans = model.mstate.tmat,
    tidy = T
) %>%
    filter(start == 1) %>%
    select(-start) %>%
    pivot_longer(cols = -time, names_to = "state", values_to = "pstate") %>%
    data.frame()


############### 7. Survival Curve Comparison ##############
pred.combined <- rbind(
    pred.mstate %>% mutate(source = "1.mstate"),
    pred.survival %>% mutate(source = "2.survival"),
    pred.weibull.j %>% mutate(source = "3.weibull.j"),
    pred.weibull.t %>% mutate(source = "4.weibull.t"),
    pred.spline.j %>% mutate(source = "5.spline.j"),
    pred.spline.t %>% mutate(source = "6.spline.t")
)

ggplot(pred.combined, aes(x = time, y = pstate, color = state)) +
    geom_line(linewidth = 1) +
    facet_wrap(~source, nrow = 3) +
    theme_cowplot()

################ 8. Instantaneous Prediction ##############
time <- 6.798084
pred.mstate[abs(pred.mstate$time - time) <= 0.01, ]
pred.survival[abs(pred.survival$time - time) <= 0.01, ]
pred.spline.t[abs(pred.spline.t$time - time) <= 0.01, ]
