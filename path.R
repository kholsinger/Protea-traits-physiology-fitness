r2 <- function(mu, obs) {
  n_obs <- length(obs)
  n_samp <- nrow(mu)
  r_2 <- numeric(n_samp)
  for (i in 1:n_samp) {
    ## mu[i, ] is the vector of predictions for MCMC sample i
    ## eps is the residual vector
    ##
    eps <- mu[i,] - obs
    r_2[i] <- var(mu[i,])/(var(mu[i,]) + var(eps))
  }
  ## The Gelman et al. approach
  ##
  ## return(median(r_2))
  ##
  ## My approach to use with partitioning
  ##
  return(mean(r_2))
}

report_r2 <- function(mu, obs, prefix="  ") {
  n_cov <- ncol(obs)
  for (i in 1:n_cov) {
    cat(prefix, colnames(obs)[i], ": ", sep="")
    cat(r2(mu[,,i], obs[,i]), "\n")
  }
}

plot_residuals <- function(mu, obs, species, traits, label) {
  Predicted <- numeric(0)
  Observed <- numeric(0)
  Species <- character(0)
  Trait <- character(0)
  n_obs <- nrow(obs)
  n_traits <- length(traits)
  for (i in 1:n_traits) {
    tmp <- mu[,,i]
    Predicted <- c(Predicted, apply(tmp, 2, mean))
    Observed <- c(Observed, obs[,traits[i]])
    Species <- c(Species, as.character(species))
    Trait <- c(Trait, rep(traits[i], n_obs))
  }
  for.plot <- data.frame(Predicted=Predicted,
                         Observed=Observed,
                         Species=Species,
                         Trait=Trait)
  p <- ggplot(for.plot, aes(x=Predicted, y=Observed, color=Species)) +
    geom_point(size=0.5) +
    geom_abline(intercept=0.0, slope=1.0) +
    facet_wrap(~ Trait)
  print(p)
  if (!DEBUG) {
    if (HARD_CENTER) {
      path <- "hard-centering-plots/"
    } else {
      path <- "soft-centering-plots/"
    }
    ggsave(filename=paste(path, label, "-observed-vs-predicted.pdf", sep=""))
  }
}

run_path_analysis <- function(dat, stan_file) {
  ## size
  ##
  size <- dat[, Size_Traits]
  size <- scale.by.column(size)

  ## performance
  ##
  performance <- dat[, Performance_Traits]
  performance <- scale.by.column(performance)

  ## traits
  ##
  trait_obs <- dat[, Morph_Traits]
  trait_obs <- scale.by.column(trait_obs)

  ## exclude individuals without performance data
  ## N.B.: will need different approach for final analysis
  ##
  include <- !is.na(fitness) & complete(size) &
    complete(performance) & complete(trait_obs)
  site <- as.numeric(dat$Site)[include]
  species <- as.numeric(dat$Species)[include]
  fitness <- fitness[include]
  i_size <- size[include,]
  performance <- performance[include,]
  trait_obs <- trait_obs[include,]

  n_species <- max(species)
  n_site <- max(site)
  n_obs <- nrow(trait_obs)
  n_trait <- ncol(trait_obs)
  n_perf <- ncol(performance)
  n_size <- ncol(size)

  cat("n_site: ", n_site, "\n",
      "n_species: ", n_species, "\n",
      "n_trait: ", n_trait, "\n", sep="")

  ## path model
  ##
  stan_data <- list(site=site,
                    species=species,
                    i_size=i_size,
                    fitness=fitness,
                    performance=performance,
                    trait_obs=trait_obs,
                    n_site=n_site,
                    n_species=n_species,
                    n_obs=n_obs,
                    n_trait=n_trait,
                    n_perf=n_perf,
                    n_size=n_size,
                    perf_scale=perf_scale,
                    beta_scale=beta_scale,
                    sigma_scale=sigma_scale)
  stan_pars <- c(report_pars,
                 "mu_perf_obs",
                 "mu_size_obs",
                 "mu_fit_obs",
                 "mu_size_species_site",
                 "log_lik_perf",
                 "log_lik_size",
                 "log_lik_fit",
                 "log_lik")
  fit <- stan(file=stan_file,
              data=stan_data,
              pars=stan_pars,
              chains=n_chains,
              iter=n_iter,
              thin=n_thin,
              control=list(adapt_delta=0.99999, max_treedepth=20))

  if (CHECK_LMER) {
    lmer_dat <- data.frame(fit=fitness,
                           bw=trait_obs[,1],
                           wd=trait_obs[,2],
                           lt=trait_obs[,3],
                           la=trait_obs[,4],
                           lma=trait_obs[,5],
                           hv=trait_obs[,6],
                           sl=trait_obs[,7],
                           sd=trait_obs[,8],
                           spi=trait_obs[,9],
                           ht=i_size[,1],
                           can=i_size[,2],
                           ks=performance[,1],
                           kl=performance[,2],
                           amax=performance[,3],
                           atot=performance[,4],
                           tran=performance[,5],
                           cond=performance[,6],
                           wue=performance[,7],
                           species=species,
                           site=site)
    fit_lmer <- stan_glmer(fitness ~ bw + wd + lt + la + lma + hv + sl + sd +
                             spi + ht + can + ks + kl + amax + atot + tran +
                             cond + wue + (1|species/site), data=lmer_dat,
                           family=poisson,
                           adapt_delta=0.999)
  }

  if (CHECK_LMER) {
    ret_list <- list(fit=fit,
                     fit_lmer=fit_lmer,
                     fitness=fitness,
                     performance=performance,
                     size=i_size,
                     trait_obs=trait_obs,
                     species=species,
                     site=site,
                     site_table=
                       unique(data.frame(Site=as.character(dat$Site[include]),
                                            number=site)),
                     species_table=
                       unique(data.frame(Species=as.character(dat$Species[include]),
                                         number=species)),
                     include=incldue)
  } else {
    ret_list <- list(fit=fit,
                     fitness=fitness,
                     size=i_size,
                     performance=performance,
                     trait_obs=trait_obs,
                     species=species,
                     site=site,
                     site_table=
                       unique(data.frame(Site=as.character(dat$Site[include]),
                                            number=site)),
                     species_table=
                       unique(data.frame(Species=as.character(dat$Species[include]),
                                         number=species)),
                     include=include)
  }
  return(ret_list)
}

run_post_path <- function(fit, dat) {
  mu <- extract(fit$fit, pars=c("mu_perf_obs",
                                "mu_size_obs",
                                "mu_size_species_site",
                                "beta_1_size_trait",
                                "beta_1_size_perf"))
  cat("\n\n")
  cat("R^2(perf):\n")
  report_r2(mu$mu_perf_obs, fit$performance)
  cat("R^2(size):\n")
  report_r2(mu$mu_size_obs, fit$size)
  if (PLOT_RESIDUAL) {
    plot_residuals(mu$mu_perf_obs, fit$performance, dat$Species[fit$include],
                   Performance_Traits, "Performance");
    plot_residuals(mu$mu_size_obs, fit$size, dat$Species[fit$include],
                   Size_Traits, "Size");
    obs <- fit$fitness
    pred <- get_posterior_mean(fit$fit, pars="mu_fit_obs")
    n_obs <- length(obs)
    Observed <- obs
    Predicted <- numeric(0)
    Method <- character(0)
    ## Stan
    ##
    for (j in 1:n_obs) {
      Predicted <- c(Predicted, exp(pred[j]))
    }
    Method <- c(Method, rep("Stan", length(obs)))
    if (CHECK_LMER) {
      ## stan_lmer()
      ##
      Observed <- c(Observed, obs)
      Predicted <- c(Predicted, fitted(fit$fit_lmer))
      Method <- c(Method, rep("stan_lmer()", length(obs)))
      for.plot <- data.frame(Observed=Observed,
                             Predicted=Predicted,
                             Method=Method)
      for.plot$Residual <- for.plot$Predicted - for.plot$Observed
      p <- ggplot(for.plot, aes(x=Predicted, y=Observed, color=Method)) +
        geom_point(size=0.5) +
        geom_abline(intercept=0.0, slope=1.0)
      print(p)
      ggsave(filename="path-observed-vs-predicted.pdf")
      p <- ggplot(for.plot, aes(x=Predicted, y=Residual, color=Method)) +
        geom_point(size=0.5) +
        geom_abline(intercept=0.0, slope=0.0)
      print(p)
      ggsave(filename="path-residual-vs-predicted.pdf")
      x <- subset(for.plot, Method=="Stan")
      y <- subset(for.plot, Method=="stan_lmer()")
      tmp <- data.frame(cbind(x[,c("Predicted")], y[,c("Predicted")]))
      colnames(tmp) <- c("Stan", "stan_lmer")
      p <- ggplot(tmp, aes(x=Stan, y=stan_lmer)) +
        geom_point(size=0.5) +
        geom_abline(intercept=0.0, slope=1.0)
      print(p)
      ggsave(filename="path-stan_lmer-vs-Stan.pdf")
    } else {
      if (HARD_CENTER) {
        path <- "hard-centering-plots/"
      } else {
        if (STOMATES_STRUCT) {
          path <- "soft-centering-plots/Stomata-structural"
        } else {
          path <- "soft-centering-plots/Stomata-physiology"
        }
      }
      Species <- dat$Species[fit$include]
      for.plot <- data.frame(Observed=Observed,
                             Predicted=Predicted,
                             Species=Species)
      for.plot$Residual <- for.plot$Predicted - for.plot$Observed
      p <- ggplot(for.plot, aes(x=Predicted, y=Observed, color=Species)) +
        geom_point(size=0.5) +
          geom_abline(intercept=0.0, slope=1.0)
      print(p)
      if (!DEBUG) {
        ggsave(filename=paste(path, "Fitness-observed-vs-predicted.pdf", sep=""))
      }
      p <- ggplot(for.plot, aes(x=Predicted, y=Residual, color=Species)) +
        geom_point(size=0.5) +
          geom_abline(intercept=0.0, slope=0.0)
      print(p)
      if (!DEBUG) {
        ggsave(filename=paste(path, "Fitness-residuals.pdf", sep=""))
      }
    }
  }
}

