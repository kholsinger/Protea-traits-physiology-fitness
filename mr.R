run_mr_analysis <- function(dat) {
  ## performance
  ##
  performance <- dat[, Performance_Traits]
  performance <- scale.by.column(performance)

  ## traits
  ##
  bw <- as.numeric(scale(dat$BW_dry_avg))
  wd <- as.numeric(scale(dat$WD))
  lt <- as.numeric(scale(dat$Thick_cm))
  la <- as.numeric(scale(dat$Area_cm2))
  lma <- as.numeric(scale(dat$LMA_g.cm2))
  sl <- as.numeric(scale(dat$Length_mm_Top))
  sd <- as.numeric(scale(dat$Density_mm2_Top))
  spi <- as.numeric(scale(dat$SPI_Top))
  trait_obs <- cbind(bw, wd, lt, la, lma, sl, sd, spi)
  stopifnot(ncol(trait_obs) == n_traits)

  ## exclude individuals without performance data
  ## N.B.: will need different approach for final analysis
  ##
  include <- complete(performance) & complete(trait_obs)
  site <- as.numeric(dat$Site)[include]
  species <- as.numeric(dat$Species)[include]
  performance <- performance[include,]
  trait_obs <- trait_obs[include,]

  n_species <- max(species)
  n_site <- max(site)
  n_obs <- nrow(trait_obs)
  n_trait <- ncol(trait_obs)
  n_perf <- ncol(performance)

  cat("n_site: ", n_site, "\n",
      "n_species: ", n_species, "\n",
      "n_trait: ", n_trait, "\n", sep="")

  ## multiple regression model
  ##
  stan_data <- list(site=site,
                    species=species,
                    performance=performance,
                    trait_obs=trait_obs,
                    n_site=n_site,
                    n_species=n_species,
                    n_obs=n_obs,
                    n_trait=n_trait,
                    n_perf=n_perf,
                    perf_scale=perf_scale,
                    beta_scale=beta_scale,
                    sigma_scale=sigma_scale)
  stan_pars <- c(report_pars,
                 "mu_obs")
  fit <- stan(file="mr.stan",
              data=stan_data,
              pars=stan_pars,
              chains=n_chains,
              iter=n_iter,
              thin=n_thin,
              control=list(adapt_delta=0.99))
  if (CHECK_LMER) {
    lmer_list <- vector("list", ncol(performance))
    if (DEBUG) {
      n_chains = 1
      n_iter = 200
    } else {
      n_chains = 4
      n_iter = 2000
    }
    for (i in 1:ncol(performance)) {
      lmer_list[[i]] <- stan_lmer(as.matrix(performance)[,i] ~ trait_obs + (1|species/site),
                                  adapt_delta=0.99,
                                  prior=normal(0.0,1.0),
                                  prior_intercept=normal(0.0,1.0),
                                  prior_aux=cauchy(0.0,5.0),
                                  prior_covariance=decov(regularization=1.0,
                                                         concentration=0.01,
                                                         shape=1.0,
                                                         scale=1.0),
                                  chains=n_chains,
                                  iter=n_iter)
    }
    ret_list <- list(fit=fit,
                     performance=performance,
                     trait_obs=trait_obs,
                     site=unique(data.frame(Site=as.character(dat$Site[include]),
                                            number=site)),
                     species=unique(data.frame(Species=as.character(dat$Species[include]),
                                               number=species)),
                     lmer_list=lmer_list)
  } else {
    ret_list <- list(fit=fit,
                     performance=performance,
                     trait_obs=trait_obs,
                     site=unique(data.frame(Site=as.character(dat$Site[include]),
                                            number=site)),
                     species=unique(data.frame(Species=as.character(dat$Species[include]),
                                               number=species)))
  }
  return(ret_list)
}

run_post_mr <- function(fit) {
  if (CHECK_LMER) {
    lmer_matrix <- matrix(nrow=n_traits, ncol=length(Performance_Traits))
    for (i in 1:length(Performance_Traits)) {
      cat(paste("\n\nstan_lmer() results for", Performance_Traits[i], "\n"))
      print(summary(fit$lmer_list[[i]], digits=3,
                    pars=c("trait_obsbw",
                           "trait_obswd",
                           "trait_obslt",
                           "trait_obsla",
                           "trait_obslma",
                           "trait_obssl",
                           "trait_obssd",
                           "trait_obsspi")))
      lmer_matrix[1,i] <- fit$lmer_list[[i]]$stan_summary["trait_obsbw", "mean"]
      lmer_matrix[2,i] <- fit$lmer_list[[i]]$stan_summary["trait_obswd", "mean"]
      lmer_matrix[3,i] <- fit$lmer_list[[i]]$stan_summary["trait_obslt", "mean"]
      lmer_matrix[4,i] <- fit$lmer_list[[i]]$stan_summary["trait_obsla", "mean"]
      lmer_matrix[5,i] <- fit$lmer_list[[i]]$stan_summary["trait_obslma", "mean"]
      lmer_matrix[6,i] <- fit$lmer_list[[i]]$stan_summary["trait_obssl", "mean"]
      lmer_matrix[7,i] <- fit$lmer_list[[i]]$stan_summary["trait_obssd", "mean"]
      lmer_matrix[8,i] <- fit$lmer_list[[i]]$stan_summary["trait_obsspi", "mean"]
    }
    cat("\n\n\nComparison of Stan and stan_lmer() results...\n")
    cat("        Stan   stan_lmer()   Difference\n")
    for (i in 1:length(Performance_Traits)) {
      cat(Performance_Traits[i], "\n")
      betas <- extract(fit$fit, pars="beta_1_perf_trait")$beta_1_perf_trait
      for (j in 1:n_traits) {
        line <- sprintf("      % 5.3f        % 5.3f       % 5.3f\n",
                        mean(betas[,j,i]), lmer_matrix[j,i], mean(betas[,j,i])-lmer_matrix[j,i])
        cat(line)
      }
    }
  }

  if (PLOT_RESIDUAL) {
    obs <- fit$performance
    pred <- extract(fit$fit, pars="mu_obs")$mu_obs
    n_obs <- nrow(obs)
    n_col <- ncol(obs)
    perf_traits <- colnames(obs)
    Observed <- numeric(0)
    Predicted <- numeric(0)
    Trait <- character(0)
    Method <- character(0)
    for (i in 1:n_col) {
      ## Stan
      ##
      Observed <- c(Observed, obs[,i])
      for (j in 1:n_obs) {
        Predicted <- c(Predicted, mean(pred[,j,i]))
      }
      Trait <- c(Trait, rep(perf_traits[i], length(obs[,i])))
      Method <- c(Method, rep("Stan", length(obs[,i])))
      ## stan_lmer()
      ##
      if (CHECK_LMER) {
        Observed <- c(Observed, obs[,i])
        Predicted <- c(Predicted, fitted(fit$lmer_list[[i]]))
        Trait <- c(Trait, rep(perf_traits[i], length(obs[,i])))
        Method <- c(Method, rep("stan_lmer()", length(obs[,i])))
      }
    }
    for.plot <- data.frame(Observed=Observed,
                           Predicted=Predicted,
                           Trait=Trait,
                           Method=Method)
    for.plot$Residual <- for.plot$Predicted - for.plot$Observed
    if (!CHECK_LMER) {
      ## Stan
      ##
      for.plot.tmp <- subset(for.plot, Method="Stan")
      p <- ggplot(for.plot.tmp, aes(x=Predicted, y=Observed)) +
        geom_point(size=0.5) +
        geom_abline(intercept=0.0, slope=1.0) +
        ggtitle("Stan") +
        facet_wrap(~ Trait)
      print(p)
      p <- ggplot(for.plot.tmp, aes(x=Predicted, y=Residual)) +
        geom_point(size=0.5) +
        geom_abline(intercept=0.0, slope=0.0) +
        ggtitle("Stan") +
        facet_wrap(~ Trait)
      print(p)
    } else if (CHECK_LMER) {
      if (0) {
        ## stan_lmer()
        ##
        for.plot <- subset(for.plot.tmp, Method="stan_lmer")
        p <- ggplot(for.plot.tmp, aes(x=Predicted, y=Observed)) +
          geom_point(size=0.5) +
          geom_abline(intercept=0.0, slope=1.0) +
          ggtitle("stan_lmer") +
          facet_wrap(~ Trait)
        print(p)
        p <- ggplot(for.plot.tmp, aes(x=Predicted, y=Residual)) +
          geom_point(size=0.5) +
          geom_abline(intercept=0.0, slope=0.0) +
          ggtitle("stan_lmer") +
          facet_wrap(~ Trait)
        print(p)
      }
      ## Combined
      ##
      p <- ggplot(for.plot, aes(x=Predicted, y=Observed, color=Method)) +
        geom_point(size=0.5) +
        geom_abline(intercept=0.0, slope=1.0) +
        facet_wrap(~ Trait)
      print(p)
      ggsave(filename="mr-observed-vs-predicted.pdf")
      p <- ggplot(for.plot, aes(x=Predicted, y=Residual, color=Method)) +
        geom_point(size=0.5) +
        geom_abline(intercept=0.0, slope=0.0) +
        facet_wrap(~ Trait)
      print(p)
      ggsave(filename="mr-residual-vs-predicted.pdf")
      x <- subset(for.plot, Method=="Stan")
      y <- subset(for.plot, Method=="stan_lmer()")
      tmp <- data.frame(cbind(x[,c("Predicted")], y[,c("Predicted", "Trait")]))
      colnames(tmp) <- c("Stan", "stan_lmer", "Trait")
      p <- ggplot(tmp, aes(x=Stan, y=stan_lmer)) +
        geom_point(size=0.5) +
        geom_abline(intercept=0.0, slope=1.0) +
        facet_wrap(~ Trait)
      print(p)
      ggsave(filename="mr-stan_lmer-vs-Stan.pdf")
      ## Compare squared residuals
      ##
      for.plot$resid.2 <- for.plot$Residual^2
      resid <- ddply(for.plot,
                     c("Trait", "Method"),
                     summarize,
                     Residual.2=mean(resid.2))
      resid <- dcast(resid, Trait ~ Method, value.var="Residual.2")
      resid$Difference <- resid$`stan_lmer()` - resid$Stan
      sink("results-mr.txt", append=TRUE)
      cat("\n\nComparison of mean squared residuals\n",
          "  A positive value of difference meands Stan code fit better\n",
          "  than corresponding stan_lmer() code.\n",
          "  N.B.: Results based on using posterior mean as prediction.\n\n")
      print(resid)
      sink()
    }
  }
}

