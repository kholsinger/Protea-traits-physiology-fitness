# Protea-traits-physiology-fitness

Data and code for analysis of associations among traits, physiology,
and fitness in five species of _Protea_.

The data and code provided here should allow anyone who is interested
to reproduce results reported in the _Annals of Botany_ paper that
links to this repository. Note: Some of the scripts were used during
development and debugging. The final structure of the running model
(represented by the configuration variables in analysis.R) may not
allow one or more of the development/debugging scripts to run as is.

## Data

* Individual_Data_for_Analysis.csv - The complete set of data used in
the analyses consisting of

    - 450 individual lines of data, each corresponding to a different
      individual plant
    - 4 sites (Observations from McGregor were not included in analyses)
    - 6 species (Observsation on _Leucadendron salignum_ and _Protea
      neriifolia_ were not included in analyses)
    - 51 variables
* Protea_Trait_Data_Full.csv - Data from Mitchell et al. (_The American
  Naturalist_ 185:525-537; 2015)

## R scripts

* analysis.R - Driver script for path analyses. Variables that affect the
  analysis that is run are in ALL CAPS near the head of the script. The
  values in the script as posted are those used to produce results
  presented in the paper. Other options invoke various configurations
  used during development and debugging.
  
* Bivariate_models_R2.R - Script to estimate R2 for bivariate
  associations. 

* loo.R - Script for leave-one-out cross validation of various models
  explored during development and debugging.

* mr.R - Multiple regression analysis of structural and performance
  traits. Used to validate path model during development.
  
* path.R - Path analysis.

* plot.R - Plotting functions to visualize certain aspects of
  posterior distributions.
  
* traits.R - A utility module to define the trait sets used and ensure
  consistency across models.
  
* Variance_Partitioning.R - Script to estimate variance partitioning.
  
## Stan scripts

* mr.stan - Multiple regression model coded in Stan to match the path
  model.
  
* path-hard-center.stan - Path model with hard centering of random
  effects. See Stan documentation for details on the difference
  between hard- and soft-centering.
  
* path.stan - Path model as used in analysis presented in the paper
  with soft-centering.
  
