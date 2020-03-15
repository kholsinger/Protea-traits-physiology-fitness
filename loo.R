library(rstan)

rm(list = ls())

load("results-soft-centering-stom-phys.Rsave")
phys <- extract(fit$fit, "log_lik_perf")

load("results-soft-centering-stom-stct.RSave")
stct <- extract(fit$fit, "log_lik_perf")

STOMATES_STRUCT <- FALSE
source("traits.R")
phys_length <- length(Performance_Traits)

STOMATES_STRUCT <- TRUE
source("traits.R")
stct_length <- length(Performance_Traits)

Trait <- character(0)
elpd_diff <- numeric(0)
se <- numeric(0)
for (i in 1:min(phys_length, stct_length)) {
  Trait[i] <- Performance_Traits[i]
  a <- loo(phys$log_lik_perf[,i,])
  b <- loo(stct$log_lik_perf[,i,])
  comparison <- loo::compare(a, b)
  elpd_diff[i] <- comparison[1]
  se[i] = comparison[2]
}
sink("results-compare.txt")
comparisons <- data.frame(Trait = Trait,
                          elpd_diff = elpd_diff,
                          se = se)
print(comparisons)
sink()
