library(tidyverse)
library(reshape2)
library(rstan)
library(rstanarm)
library(ggplot2)
library(bayesplot)

rm(list=ls())

## functions for multiple regression of performance on traits
##
source("mr.R")
## functions for full path analysis
##
source("path.R")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

perf_scale <- 1.0
beta_scale <- 1.0
sigma_scale <- 5.0

DEBUG <- FALSE

## set STOMATES_STRUCT <- TRUE to include stomates in structural traits,
## otherwise include them in physiological traits
##
STOMATES_STRUCT <- TRUE

## set HARD_CENTER <- TRUE for hard centering of variable intercepts
## set HARD_CENTER <- FALSE for soft centering of variable intercepts
##
HARD_CENTER <- FALSE
if (HARD_CENTER) {
  stan_file="path-hard-center.stan"
} else {
  stan_file="path.stan"
}

## set REGRESSION <- TRUE for regression of performance traits on morphological
## traits
## set REGRESSION <- FALSE for full path analysis
##
REGRESSION <- FALSE

## N.B.: PLOT_RESIDUAL must be TRUE in order to get comparison of mean squared
## residuals from Stan and stan_lmer() code
##
PLOT_RESIDUAL <- TRUE

## N.B: code in path.R must be adjusted manually to compare path coefficients
## from Stan code with regression coefficients from stan_lmer()
##
CHECK_LMER <- FALSE
if (CHECK_LMER) {
  library(rstanarm)
}

if (REGRESSION) {
  ## N.B.: This is a horrible hack. If the number of traits included is changed,
  ## this has to be changed to match. Execution will stop if this doesn't match
  ## the number of columns specified
  ##
  n_traits <- 9

  PATH <- FALSE
} else {
  PATH <- TRUE
}

## total number of posterior samples desired
if (DEBUG) {
  n_sample <- 50
  n_chains <- 1
  n_thin <- 1
} else {
  n_sample <- 5000
  n_chains <- 4
  n_thin <- 1
}
## calculate required no. of iterations per chain
## half of iterations are burnin
n_iter <- (n_sample*2/n_chains)*n_thin

scale.by.column <- function(dat) {
  for (i in 1:ncol(dat)) {
    dat[,i] <- as.numeric(scale(dat[,i]))
  }
  return(dat)
}

complete.row <- function(row) {
  return(sum(is.na(row)) == 0)
}

complete <- function(dat) {
  return(apply(dat, 1, complete.row))
}

dat <- read.csv("Individual_Data_for_Analysis.csv", na.strings=c("NA", "."))

## Pare down the data to only have Protea, and to exclude the "bad" Ks data
##
dat <- subset(dat, Species != "LDSA")
dat <- subset(dat, Site_ID != "McGregor")
dat <- subset(dat, Site_ID != "Jonaskop_2")

dat <- dat %>% select(Site, Species, Ind_No, Age, Basal_cm, Height_in,
                      Growth_Rate, Canopy_area_in2, Total_Fruit_No,
                      HV_m2_cm2, WD, Ks, KL_Units, BW_dry_avg, Photo, Trmmol,
                      Cond, Thick_cm, Area_cm2, LMA_g.cm2, LWR, LD_g.cm3,
                      Density_mm2_Top, Length_mm_Top, SPI_Top,
                      Height_in, Canopy_area_in2, Total_Fruit_No) %>%
  mutate(wue = Photo/Cond,
         TotAssim = Photo/HV_m2_cm2) %>%
  rename(KL = KL_Units)

dat <- droplevels(dat)

source("traits.R")

fitness <- dat$Total_Fruit_No

if (REGRESSION) {
  ## parameters for default reporting from Stan
  ##
  report_pars <- c("beta_1_perf_trait",
                   "sigma_perf",
                   "sigma_perf_species",
                   "sigma_perf_species_site")

  fit <- run_mr_analysis(dat)

  old_opts <- options(max.print=10000)

  sink("results-mr.txt", append=FALSE)
  print(fit$fit, pars=report_pars, digits=3)
  run_post_mr(fit)
  sink()

  options(old_opts)
} else if (PATH) {
  ## parameters for default reporting from Stan
  ##
  report_pars <- c("beta_1_fit_trait",
                   "beta_1_fit_size",
                   "beta_1_fit_perf",
                   "beta_1_perf_trait",
                   "beta_1_size_trait",
                   "beta_1_size_perf",
                   "sigma_perf",
                   "sigma_perf_species",
                   "sigma_perf_species_site",
                   "sigma_size",
                   "sigma_size_species",
                   "sigma_size_species_site",
                   "sigma_fit_species",
                   "sigma_fit_species_site")

  fit <- run_path_analysis(dat, stan_file)

  old_opts <- options(max.print=10000)

  if (HARD_CENTER) {
    base_name <- "results-hard-centering"
  } else {
    base_name <- "results-soft-centering"
  }
  if (STOMATES_STRUCT) {
    base_name <- str_c(base_name, "-stom-stct")
  } else {
    base_name <- str_c(base_name, "-stom-phys")
  }
  if (!DEBUG) {
    sink(paste(base_name, ".txt", sep=""), append=FALSE)
  }
  if (HARD_CENTER) {
    cat("Hard centering of variable intercepts...\n\n")
  } else {
    cat("Soft centering of variable intercepts...\n\n")
  }
  cat("Traits:      ", Morph_Traits, "\n")
  cat("Performance: ", Performance_Traits, "\n")
  cat("Size:        ", Size_Traits, "\n")
  cat("Fitness:     Total_Fruit_No\n\n")
  cat("Matrix indexing...\n",
      "beta_1_perf_trait[trait,performance]\n",
      "beta_1_size_trait[trait,size]\n",
      "beta_1_size_perf[performance,size]\n\n\n")
  print(fit$fit, pars=report_pars, digits=3)
  run_post_path(fit, dat)
  if (!DEBUG) {
    sink()
    save(fit, STOMATES_STRUCT, file=paste(base_name, ".Rsave", sep=""))
  }

  options(old_opts)
}
