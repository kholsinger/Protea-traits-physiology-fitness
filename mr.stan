data {
  int<lower=0> n_site;
  int<lower=0> n_species;
  int<lower=0> n_obs;
  int<lower=0> n_trait;
  int<lower=0> n_perf;
  int<lower=0> site[n_obs];
  int<lower=0> species[n_obs];
  real<lower=0> perf_scale;
  real<lower=0> beta_scale;
  real<lower=0> sigma_scale;
  vector[n_perf] performance[n_obs];
  vector[n_trait] trait_obs[n_obs];
}

parameters {
  vector[n_perf] beta_0_perf;
  matrix[n_perf,n_species] beta_perf_species;
  vector[n_perf] beta_perf_species_site[n_species, n_site];
  vector<lower=0>[n_perf] sigma_perf;
  vector<lower=0>[n_perf] sigma_perf_species;
  vector<lower=0>[n_perf] sigma_perf_species_site;
  vector[n_perf] beta_1_perf_trait[n_trait];
}

transformed parameters {
  vector[n_perf] mu_obs[n_obs];
  matrix[n_perf,n_species] mu_perf_species;
  vector[n_perf] mu_perf_species_site[n_species, n_site];

  for (i in 1:n_perf) {
    for (j in 1:n_species) {
      mu_perf_species[i,j] = beta_0_perf[i] + beta_perf_species[i,j];
      for (k in 1:n_site) {
        mu_perf_species_site[j,k][i] = mu_perf_species[i,j] +
                                       beta_perf_species_site[j,k][i];
      }
    }
  }

  // expectation for individual observation
  //    site and species determine intercept
  //    shared regression coefficients for traits across species
  //
  for (i in 1:n_perf) {
    for (j in 1:n_obs) {
      mu_obs[j,i] = mu_perf_species_site[species[j], site[j]][i];
      for (k in 1:n_trait) {
        mu_obs[j,i] = mu_obs[j,i] + beta_1_perf_trait[k,i]*trait_obs[j,k];
      }
    }
  }

}

model {
  // observational model
  //
  for (i in 1:n_perf) {
    for (j in 1:n_obs) {
      performance[j,i] ~ normal(mu_obs[j,i], sigma_perf[i]);
    }
  }

  // priors
  //
  for (i in 1:n_perf) {
    beta_0_perf[i] ~ normal(0.0, perf_scale);
  }

  for (i in 1:n_perf) {
    for (j in 1:n_species) {
//      beta_perf_species[i,j] ~ normal(0.0, beta_scale);
      beta_perf_species[i,j] ~ normal(0.0, sigma_perf_species[i]);
      for (k in 1:n_site) {
//        beta_perf_species_site[j,k][i] ~ normal(0.0, beta_scale);
        beta_perf_species_site[j,k][i] ~ normal(0.0, sigma_perf_species_site[i]);
      }
    }
  }

  for (i in 1:n_perf) {
    beta_1_perf_trait[i] ~ normal(0.0, beta_scale);
    sigma_perf[i] ~ cauchy(0, sigma_scale);
    sigma_perf_species[i] ~ cauchy(0, sigma_scale);
    sigma_perf_species_site[i] ~ cauchy(0, sigma_scale);
  }
}
