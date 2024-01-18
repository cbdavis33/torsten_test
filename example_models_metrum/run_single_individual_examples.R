rm(list = ls())
cat("\014")

library(mrgsolve)
library(tidybayes)
library(posterior)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("Torsten/cmdstan")

# Analytical
model_analytical <- cmdstan_model("example_models_metrum/pk2cpt/pk2cpt.stan")

fit_analytical <- model_analytical$sample(
  data = "example_models_metrum/pk2cpt/pk2cpt.data.R", 
  init = "example_models_metrum/pk2cpt/pk2cpt.init.R", 
  seed = 3191951,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  adapt_delta = 0.80)

# Matrix exponential
model_mat_exp <- cmdstan_model(
  "example_models_metrum/pk2cpt_linode/pk2cpt_linode.stan")

fit_mat_exp <- model_mat_exp$sample(
  data = "example_models_metrum/pk2cpt_linode/pk2cpt_linode.data.R", 
  init = "example_models_metrum/pk2cpt_linode/pk2cpt_linode.init.R", 
  seed = 3191951,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  adapt_delta = 0.80)

# Fit pk2cpt data with pk2cpt_linode model just to check
fit_mat_exp_2 <- model_mat_exp$sample(
  data = "example_models_metrum/pk2cpt/pk2cpt.data.R", 
  init = "example_models_metrum/pk2cpt/pk2cpt.init.R", 
  seed = 3191951,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  adapt_delta = 0.80)

# Fit pk2cpt_ode data with pk2cpt_linode model just to check
fit_mat_exp_3 <- model_mat_exp$sample(
  data = "example_models_metrum/pk2cpt_ode/pk2cpt_ode.data.R", 
  init = "example_models_metrum/pk2cpt_ode/pk2cpt_ode.init.R", 
  seed = 3191951,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  adapt_delta = 0.80)


# rk45
model_rk45 <- cmdstan_model(
  "example_models_metrum/pk2cpt_ode/pk2cpt_ode.stan")

stan_data_rk45 <- 
  rstan::read_rdump("example_models_metrum/pk2cpt_ode/pk2cpt_ode.data.R")

stan_data_rk45$prior_only <- 1

fit_rk45 <- model_rk45$sample(
  data = stan_data_rk45, 
  init = "example_models_metrum/pk2cpt_ode/pk2cpt_ode.init.R", 
  seed = 3191951,
  chains = 1,
  parallel_chains = 1,
  iter_warmup = 1,
  iter_sampling = 1,
  adapt_delta = 0.80)

# Fit pk2cpt data with pk2cpt_ode model just to check
stan_data_rk45 <- 
  rstan::read_rdump("example_models_metrum/pk2cpt/pk2cpt.data.R")

stan_data_rk45$prior_only <- 0

fit_rk45_2 <- model_rk45$sample(
  data = stan_data_rk45,
  init = "example_models_metrum/pk2cpt/pk2cpt.init.R",
  seed = 3191951,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  adapt_delta = 0.80)

# Fit pk2cpt_linode data with pk2cpt_ode model just to check
stan_data_rk45 <- 
  rstan::read_rdump("example_models_metrum/pk2cpt_linode/pk2cpt_linode.data.R")

stan_data_rk45$prior_only <- 0

fit_rk45 <- model_rk45$sample(
  data = stan_data_rk45,
  init = "example_models_metrum/pk2cpt_linode/pk2cpt_linode.init.R",
  seed = 3191951,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000,
  adapt_delta = 0.80)
