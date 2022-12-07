# A collection of benchmarks for logistic SuSiE

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed. # nolint


library(future)
library(future.callr)
plan(callr)

# Set target options:
tar_option_set(
  packages = c("tibble", "logisticsusie", "dplyr", "tidyr", "purrr", "ggplot2", "susieR", "stringr"), # packages that your targets need to run
  format = "rds", # default storage format
  workspace_on_error = TRUE # Save a workspace file for a target that errors out.
  # Set other options as needed.
)


tar_source()


# Generate y --------
simulate_null <- function(X){
  # simulate across multiple settings
  beta0 <- c(-2, -1, -.5, 0)
  re_var <- c(0)
  reps <- 1:1

  # generate simulations
  sims <- tidyr::crossing(beta0=beta0, re_var=re_var, rep=reps) %>%
    dplyr::mutate(X_name = attributes(X)$name) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sim = list(logisticsusie:::sim_y_null(X, beta0, re_var))) %>%
    ungroup()
  return(sims)
}

simulate_half_normal <- function(X){
  # simulate across multiple settings
  beta0 <- c(-2, -1, -.5, 0)
  beta_sigma <- c(0.2, 0.4, 0.6, 0.8, 1.0, 2.0) / sqrt(2/pi)
  L <- c(1, 3, 5)
  reps <- 1:1

  # generate simulations
  sims <- tidyr::crossing(beta0=beta0, beta_sigma=beta_sigma, L = L, rep=reps) %>%
    dplyr::mutate(X_name = attributes(X)$name) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(beta = list(abs(rnorm(L) * beta_sigma))) %>%
    dplyr::mutate(sim = list(logisticsusie:::sim_y_susie(X, beta0, beta))) %>%
    ungroup()
  return(sims)
}

simulate_half_normal_re <- function(X){
  # simulate across multiple settings
  beta0 <- c(-2, -1, -.5, 0)
  beta_sigma <- c(0.2, 0.4, 0.6, 0.8, 1.0, 2.0) / sqrt(2/pi)
  L <- c(1, 3, 5)
  reps <- 1:1

  # generate simulations
  sims <- tidyr::crossing(beta0=beta0, beta_sigma=beta_sigma, L = L, rep=reps) %>%
    dplyr::mutate(X_name = attributes(X)$name) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(beta = list(abs(rnorm(L) * beta_sigma))) %>%
    dplyr::mutate(sim = list(logisticsusie:::sim_y_susie(X, beta0, beta, re_var=1.0))) %>%
    ungroup()
  return(sims)
}

simulate_half_normal_L1 <- function(X){
  # simulate across multiple settings
  beta0 <- c(-2)
  beta_sigma <- c(0.2, 0.4, 0.6, 0.8, 1.0, 2.0) / sqrt(2/pi)
  L <- c(1)
  reps <- 1:20

  # generate simulations
  sims <- tidyr::crossing(beta0=beta0, beta_sigma=beta_sigma, L = L, rep=reps) %>%
    dplyr::mutate(X_name = attributes(X)$name) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(beta = list(abs(rnorm(L) * beta_sigma))) %>%
    dplyr::mutate(sim = list(logisticsusie:::sim_y_susie(X, beta0, beta))) %>%
    ungroup()
  return(sims)
}

single_ser_sim <- function(X, beta0, beta, reps=20){
  # simulate across multiple settings
  reps <- 1:reps

  # generate simulations
  sims <- tidyr::crossing(beta0=beta0, beta=beta, rep=reps) %>%
    dplyr::mutate(X_name = attributes(X)$name) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sim = list(logisticsusie:::sim_y_susie(X, beta0, beta))) %>%
    ungroup()
  return(sims)
}


# Target definitions --------

# Yusha example -------
# analysis function
# take yushas example, take the n_top marginal enrichments and fit logistic ibss
yusha_example_ibss <- function(L, n_top, maxit=50){
  example <- readRDS('data/yusha_sc_tumor/pdac_example.rds')
  bindata <- with(example, gseasusie::prep_binary_data(genesets, data, thresh = 0.01))
  gs_names <- colnames(bindata$X)
  n_gene_sets <- dim(bindata$X)[2]

  ora <- with(bindata, gseasusie::fit_ora(X, y))
  top_idx <- head(order(ora$pFishersExact), n_top)
  bindata_sub <- bindata
  bindata_sub$X <- bindata_sub$X[, top_idx]

  fit <- with(bindata_sub, logisticsusie::binsusie(X, y, L=L, estimate_prior_variance = F, prior_variance = 1))

  ibss_vb <- with(bindata_sub, logisticsusie::ibss_from_ser(
    X, y, L=L, prior_variance = 1, ser_function = logisticsusie::fit_bin_ser, maxit = maxit))

  ibss_vbc <- with(bindata_sub, logisticsusie::ibss_from_ser(
    X, y, L=L, prior_variance = 1, ser_function = logisticsusie::fit_bin_ser_corrected, maxit = maxit))

  ibss_uvb <- with(bindata_sub, logisticsusie::ibss_from_ser(
    X, y, L=L, prior_variance = 1, ser_function = logisticsusie::fit_uvb_ser, maxit = maxit))

  res <- list(
    fit = fit,
    ibss_vb = ibss_vb,
    ibss_vbc = ibss_vbc,
    ibss_uvb = ibss_uvb
  )
  return(res)
}

yusha_example_targets <- list(
  tar_target(yusha_ibss_l3_n100, yusha_example_ibss(L=3, n_top=100)),
  tar_target(yusha_ibss_l5_n100, yusha_example_ibss(L=5, n_top=100)),
  tar_target(yusha_ibss_l3_n500, yusha_example_ibss(L=3, n_top=500)),
  tar_target(yusha_ibss_l5_n500, yusha_example_ibss(L=5, n_top=500))
)

# All targets -------
list(
  make_X_targets(),
  half_normal_target,
  constant_sim_target,
  website
)

