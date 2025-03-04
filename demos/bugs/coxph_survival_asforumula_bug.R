# This is a minimal, reproducible example of a bug the definition of formulas for CoxPH using `as.formula(string)` when the call to `coxph` is inside a function (for example, when using lapply).
library(survival)

data(cancer, package="survival")
ndata <- expand.grid(sex=c("F", "M"), age=c(60, 80))

#### The following steps break: ####
my_formula <- as.formula("Surv(futime, death) ~ sex + age")

# Fit Cox PH model and run survfit
fit_cox <- function(df) {
    msfit <- coxph(my_formula, data=df, id=id)
    mssurv <- survfit(msfit, newdata=ndata)
}

fit_cox(mgus2)

#### The following change makes it work! ####
# Moving the my_formula definition into the fit_cox function

# Fit Cox PH model and run survfit
fit_cox <- function(df) {
    my_formula <- as.formula("Surv(futime, death) ~ sex + age")
    msfit <- coxph(my_formula, data=df, id=id)
    mssurv <- survfit(msfit, newdata=ndata)
}

fit_cox(mgus2)

