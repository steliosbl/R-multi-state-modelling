# This is a minimal, reproducible example of a bug impacting summary.survfit on multi-state data in the R survival package version 3.8.3
# Relevant github issue: https://github.com/therneau/survival/issues/305
library(survival)

# Prepare the mgus2 dataset.
mgus2$etime <- with(mgus2, ifelse(pstat == 0, futime, ptime))
event <- with(mgus2, ifelse(pstat == 0, 2 * death, 1))
mgus2$event <- factor(event, 0:2, labels = c("censor", "pcm", "death"))

# Fit Cox PH model
msfit <- coxph(
    Surv(etime, event) ~ age + sex,
    data = mgus2, id = id
)

# The following ndata causes an error
ndata <- mgus2[1:2, c("age", "sex")]
sf <- survfit(msfit, ndata)
summary(sf, data.frame = TRUE)
# > Error in (function (..., row.names = NULL, check.rows = FALSE, check.names = TRUE,  :
#  arguments imply differing number of rows: 1284, 1926

# We try all "ndata" lengths from 1-30 and find that multiples of 3 work
# But any other lengths cause the error
working_nrow <- unlist(lapply(1:30, function(i) {
    tryCatch(
        {
            ndata <- mgus2[1:i, c("age", "sex")]
            sf <- survfit(msfit, ndata)
            summary(sf, data.frame = TRUE)
            return(i)
        },
        error = function(e) NULL # Return NULL on error
    )
}))
working_nrow
