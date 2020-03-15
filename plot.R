library(tidyverse)
library(stringr)
library(rstan)
library(bayesplot)

rm(list = ls())

load("results-soft-centering-stom-phys.Rsave")
STOMATES_STRUCT <- FALSE
source("traits.R")

titles <- gsub("_.*", "", Performance_Traits)
columns <- gsub("_.*", "", Morph_Traits)
ordered_pars <- names(sort(get_posterior_mean(fit$fit, pars = "beta_1_fit_trait")[, 5],
                           decreasing = TRUE))
ordered_pars <- as.numeric(gsub(".*\\[(.*)\\]", "\\1", ordered_pars))
for (i in 1:length(Performance_Traits)) {
  df <- as.data.frame(fit$fit) %>%
    select(matches(str_c("beta_1_perf_trait\\[[0-9],[", i, "]\\]")))
  colnames(df) <- columns
  p <- mcmc_intervals(df, pars = columns[ordered_pars],
                      prob = 0.8, prob_outer = 0.95) + ggtitle(Performance_Traits[i])
  print(p)
  ggsave(str_c(Performance_Traits[i], "-stom-phys.pdf"))
}

load("results-soft-centering-stom-stct.RSave")
STOMATES_STRUCT <- TRUE
source("traits.R")

titles <- gsub("_.*", "", Performance_Traits)
columns <- gsub("_.*", "", Morph_Traits)
ordered_pars <- names(sort(get_posterior_mean(fit$fit, pars = "beta_1_fit_trait")[, 5],
                           decreasing = TRUE))
ordered_pars <- as.numeric(gsub(".*\\[(.*)\\]", "\\1", ordered_pars))
for (i in 1:length(Performance_Traits)) {
  df <- as.data.frame(fit$fit) %>%
    select(matches(str_c("beta_1_perf_trait\\[[0-9],[", i, "]\\]")))
  colnames(df) <- columns
  p <- mcmc_intervals(df, pars = columns[ordered_pars],
                      prob = 0.8, prob_outer = 0.95) + ggtitle(Performance_Traits[i])
  print(p)
  ggsave(str_c(Performance_Traits[i], "-stom-stct.pdf"))
}




