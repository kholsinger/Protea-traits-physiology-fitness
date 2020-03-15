data {
  int<lower=0> n_site;
  int<lower=0> n_species;
  int<lower=0> n_obs;
  int<lower=0> n_trait;
  int<lower=0> n_perf;
  int<lower=0> n_size;
  int<lower=0> site[n_obs];
  int<lower=0> species[n_obs];
  real<lower=0> perf_scale;
  real<lower=0> beta_scale;
  real<lower=0> sigma_scale;
  int fitness[n_obs];
  vector[n_size] i_size[n_obs];
  vector[n_perf] performance[n_obs];
  vector[n_trait] trait_obs[n_obs];
}

transformed data {
  vector[n_species] species_1_0;
  vector[n_site] site_1_0;

  for (i in 1:n_species) {
    species_1_0[i] = 1.0;
  }
  for (i in 1:n_site) {
    site_1_0[i] = 1.0;
  }
}

parameters {
  real beta_0_fit;
  vector[n_perf] beta_0_perf;
  vector[n_size] beta_0_size;
  vector<lower=0>[n_perf] sigma_perf;
  vector<lower=0>[n_perf] sigma_perf_species;
  vector<lower=0>[n_perf] sigma_perf_species_site;
  vector<lower=0>[n_size] sigma_size;
  vector<lower=0>[n_size] sigma_size_species;
  vector<lower=0>[n_size] sigma_size_species_site;
  real<lower=0> sigma_fit_species;
  real<lower=0> sigma_fit_species_site;
  real beta_1_fit_trait[n_trait];
  real beta_1_fit_size[n_size];
  real beta_1_fit_perf[n_perf];
  vector[n_size] beta_1_size_perf[n_perf];
  vector[n_size] beta_1_size_trait[n_trait];
  vector[n_perf] beta_1_perf_trait[n_trait];
  // hard centering variable intercepts
  //
  // performance
  //
  vector<lower=0>[n_perf] gamma_perf_species;
  vector<lower=0>[n_perf] gamma_perf_species_site;
  simplex[n_species] alpha_perf_species[n_perf];
  simplex[n_site] alpha_perf_species_site[n_perf, n_species];
  // size
  //
  vector<lower=0>[n_size] gamma_size_species;
  vector<lower=0>[n_size] gamma_size_species_site;
  simplex[n_species] alpha_size_species[n_size];
  simplex[n_site] alpha_size_species_site[n_size, n_species];
  // species
  //
  real<lower=0> gamma_fit_species;
  real<lower=0> gamma_fit_species_site;
  simplex[n_species] alpha_fit_species;
  simplex[n_site] alpha_fit_species_site[n_species];
}

transformed parameters {
  vector[n_obs] mu_fit_obs;
  vector[n_perf] mu_perf_obs[n_obs];
  vector[n_size] mu_size_obs[n_obs];
  vector[n_species] mu_fit_species;
  matrix[n_perf,n_species] mu_perf_species;
  matrix[n_size,n_species] mu_size_species;
  real mu_fit_species_site[n_species, n_site];
  vector[n_perf] mu_perf_species_site[n_species, n_site];
  vector[n_size] mu_size_species_site[n_species, n_site];
  // variable intercepts
  //
  matrix[n_perf, n_species] beta_perf_species;
  vector[n_perf] beta_perf_species_site[n_species, n_site];
  matrix[n_size, n_species] beta_size_species;
  vector[n_size] beta_size_species_site[n_species, n_site];
  vector[n_species] beta_fit_species;
  matrix[n_species, n_site] beta_fit_species_site;

  real t_sum;

  // variable intercepts
  //
  // performance
  //
  for (i in 1:n_perf) {
    for (j in 1:n_species) {
      beta_perf_species[i,j] = gamma_perf_species[i]
        * (alpha_perf_species[i,j] - 1.0/n_species);
      for (k in 1:n_site) {
        // N.B.: indexing reversed because alpha_* is an n_perf x n_species
        // array of simplexes (length n_site)
        //
        beta_perf_species_site[j,k][i] =
          gamma_perf_species_site[i] * (alpha_perf_species_site[i,j][k]
                                        - 1.0/n_site);
      }
    }
  }
  // size
  //
  for (i in 1:n_size) {
    for (j in 1:n_species) {
      beta_size_species[i,j] = gamma_size_species[i]
        * (alpha_size_species[i,j] - 1.0/n_species);
      for (k in 1:n_site) {
        // N.B.: indexing reversed because alpha_* is an n_perf x n_species
        // array of simplexes (length n_site)
        //
        beta_size_species_site[j,k][i] =
          gamma_size_species_site[i] * (alpha_size_species_site[i,j][k]
                                        - 1.0/n_site);
      }
    }
  }
  // species
  //
  beta_fit_species = gamma_fit_species*alpha_fit_species;
  for (j in 1:n_species) {
    beta_fit_species_site[j]
      = (gamma_fit_species_site*alpha_fit_species_site[j])';
  }

  // site within species effect for performance
  //
  for (i in 1:n_perf) {
    for (j in 1:n_species) {
      mu_perf_species[i,j] = beta_0_perf[i] + beta_perf_species[i,j];
      for (k in 1:n_site) {
        mu_perf_species_site[j,k][i] = mu_perf_species[i,j] +
                                       beta_perf_species_site[j,k][i];
      }
    }
  }

  // site within species effect for size
  //
  for (i in 1:n_size) {
    for (j in 1:n_species) {
      mu_size_species[i,j] = beta_0_size[i] + beta_size_species[i,j];
      for (k in 1:n_site) {
        mu_size_species_site[j,k][i] = mu_size_species[i,j] +
                                       beta_size_species_site[j,k][i];
      }
    }
  }

  // site withn species effect for fitness
  //
  for (j in 1:n_species) {
    mu_fit_species[j] = beta_0_fit + beta_fit_species[j];
    for (k in 1:n_site) {
      mu_fit_species_site[j,k] = mu_fit_species[j] + beta_fit_species_site[j,k];
    }
  }

  // expectation for individual performance observation
  //    site and species determine intercept
  //    shared regression coefficients for traits across species
  //
  for (i in 1:n_perf) {
    for (j in 1:n_obs) {
      mu_perf_obs[j,i] = mu_perf_species_site[species[j], site[j]][i];
      for (k in 1:n_trait) {
        mu_perf_obs[j,i] = mu_perf_obs[j,i] + beta_1_perf_trait[k,i]*trait_obs[j,k];
      }
    }
  }

  // expectation for individual size observation
  //    site and species determine intercept
  //    shared regression coefficients for performance across species
  //
  for (i in 1:n_size) {
    for (j in 1:n_obs) {
      mu_size_obs[j,i] = mu_size_species_site[species[j], site[j]][i];
      for (k in 1:n_trait) {
        mu_size_obs[j,i] = mu_size_obs[j,i] + beta_1_size_trait[k,i]*trait_obs[j,k];
      }
      for (k in 1:n_perf) {
        mu_size_obs[j,i] = mu_size_obs[j,i] + beta_1_size_perf[k,i]*performance[j,k];
      }
    }
  }

  // expectaton for individual fitness observation
  //    site and species determine intercept
  //    shared regression coefficients for traits and size across species
  for (j in 1:n_obs) {
    mu_fit_obs[j] = mu_fit_species_site[species[j], site[j]];
    for (k in 1:n_trait) {
      mu_fit_obs[j] = mu_fit_obs[j] + beta_1_fit_trait[k]*trait_obs[j,k];
    }
    for (k in 1:n_size) {
      mu_fit_obs[j] = mu_fit_obs[j] + beta_1_fit_size[k]*i_size[j,k];
    }
    for (k in 1:n_perf) {
      mu_fit_obs[j] = mu_fit_obs[j] + beta_1_fit_perf[k]*performance[j,k];
    }
  }

}

model {
  // observational model
  //
  // fitness
  //
  for (j in 1:n_obs) {
    fitness[j] ~ poisson(exp(mu_fit_obs[j]));
  }
  // performance
  //
  for (i in 1:n_perf) {
    for (j in 1:n_obs) {
      performance[j,i] ~ normal(mu_perf_obs[j,i], sigma_perf[i]);
    }
  }
  // size
  //
  for (i in 1:n_size) {
    for (j in 1:n_obs) {
      i_size[j,i] ~ normal(mu_size_obs[j,i], sigma_size[i]);
    }
  }

  // priors
  //
  // performance
  //
  for (i in 1:n_perf) {
    beta_0_perf[i] ~ normal(0.0, perf_scale);
  }
  // size
  //
  for (i in 1:n_size) {
    beta_0_size[i] ~ normal(0.0, perf_scale);
  }
  // fitness
  //
  beta_0_fit ~ normal(0.0, perf_scale);

  // variable intercepts
  //
  // performance
  //
  gamma_perf_species ~ normal(0.0, sigma_perf_species);
  gamma_perf_species_site ~ normal(0.0, sigma_perf_species_site);
  for (i in 1:n_perf) {
    alpha_perf_species[i] ~ dirichlet(species_1_0);
    for (j in 1:n_species) {
      alpha_perf_species_site[i,j] ~ dirichlet(site_1_0);
    }
  }
  // size
  //
  gamma_size_species ~ normal(0.0, sigma_size_species);
  gamma_size_species_site ~ normal(0.0, sigma_size_species_site);
  for (i in 1:n_size) {
    alpha_size_species[i] ~ dirichlet(species_1_0);
    for (j in 1:n_species) {
      alpha_size_species_site[i,j] ~ dirichlet(site_1_0);
    }
  }
  // fitness
  //
  gamma_fit_species ~ normal(0.0, sigma_fit_species);
  alpha_fit_species ~ dirichlet(species_1_0);
  gamma_fit_species_site ~ normal(0.0, sigma_fit_species_site);
  for (j in 1:n_species) {
    alpha_fit_species_site[j] ~ dirichlet(site_1_0);
  }

  // regression coefficients
  //
  // performance
  //
  for (i in 1:n_trait) {
    beta_1_perf_trait[i] ~ normal(0.0, beta_scale);
  }
  // size
  //
  for (i in 1:n_perf) {
    beta_1_size_perf[i] ~ normal(0.0, beta_scale);
  }
  for (i in 1:n_trait) {
    beta_1_size_trait[i] ~ normal(0.0, beta_scale);
  }
  // fitness
  //
  for (i in 1:n_trait) {
      beta_1_fit_trait[i] ~ normal(0.0, beta_scale);
  }
  for (i in 1:n_perf) {
      beta_1_fit_perf[i] ~ normal(0.0, beta_scale);
  }
  for (i in 1:n_size) {
      beta_1_fit_size[i] ~ normal(0.0, beta_scale);
  }

  // variable intercept standard deviations
  //
  // performance
  //
  for (i in 1:n_perf) {
    sigma_perf[i] ~ cauchy(0, sigma_scale);
    sigma_perf_species[i] ~ cauchy(0, sigma_scale);
    sigma_perf_species_site[i] ~ cauchy(0, sigma_scale);
  }
  // size
  //
  for (i in 1:n_size) {
    sigma_size[i] ~ cauchy(0, sigma_scale);
    sigma_size_species[i] ~ cauchy(0, sigma_scale);
    sigma_size_species_site[i] ~ cauchy(0, sigma_scale);
  }
  // fitness
  //
  sigma_fit_species ~ cauchy(0, sigma_scale);
  sigma_fit_species_site ~ cauchy(0, sigma_scale);
}
