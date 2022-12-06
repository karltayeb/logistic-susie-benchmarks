# A collection of benchmarks for logistic SuSiE

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed. # nolint


library(future)
library(future.callr)
plan(callr)

# Set target options:
tar_option_set(
  packages = c("tibble", "logisticsusie", "dplyr", "tidyr", "purrr", "ggplot2", "susieR"), # packages that your targets need to run
  format = "rds", # default storage format
  workspace_on_error = TRUE # Save a workspace file for a target that errors out.
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
# options(clustermq.scheduler = "multiprocess")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# METHODS ---------
# all the fit functions take a simplified interface
# accepting ONLY two arguments: X and y

ser_fit_functions <- tidyr::tribble(
  ~fit_method, ~fit_fun, ~fit_args,
  'vb_ser', 'fit_bin_ser', list(),
  'veb_ser', 'fit_veb_ser', list(),
  'uvb_ser', 'fit_uvb_ser', list(),
  'quad_ser', 'fit_quad_ser', list(),
  'glm_ser', 'fit_glm_ser', list(),
  'uvb_ser_re', 'fit_uvb_ser_re', list(),
  'vb_ser_corrected', 'fit_bin_ser_corrected', list(),
  'vb_ser2', 'fit_bin_ser', list(estimate_prior_variance=T),
  'uvb_ser2', 'fit_uvb_ser', list(estimate_prior_variance=T),
  'glm_ser2', 'fit_glm_ser', list(estimate_prior_variance=T),
  'jj_abf_ser', 'fit_jjabf_ser', list(estimate_prior_variance=F)
)

# all the fit functions take a simplified interface
# accepting three arguments: X, y, L
# NOTE: by quoting the function names we avoid rerunning every time the functions change
# this is helpful while developing the `logisticsusie` package, but it means we need to be
# sure to MANUALLY re-trigger runs when there is a non-trivial change to the method
logistic_ibss_functions <- tidyr::tribble(
  ~fit_method, ~fit_fun, ~fit_args,
  'ibss_vb_L5', 'fit_ibss_vb', list(L=5),
  'ibss_vb2_L5', 'fit_ibss_vb2', list(L=5),
  'ibss_vbc_L5', 'fit_ibss_vbc', list(L=5),
  'ibss_uvb_L5', 'fit_ibss_uvb', list(L=5),
  'ibss_uvb2_L5', 'fit_ibss_uvb2', list(L=5),
  'ibss_glm_L5', 'fit_ibss_glm', list(L=5),
  'binsusie_L5', 'fit_binsusie_wrapped', list(L=5, estimate_prior_variance=F, prior_variance=1),
  'binsusie2_L5', 'fit_binsusie_wrapped', list(L=5, estimate_prior_variance=T, prior_variance=1),
  'ibss2m_L5', 'ibss2m', list(L=5, maxit=50, track_elbo=T)
  #'ibss2m2_L5', 'ibss2m', list(L=5, maxit=50, estimate_prior_variance=T, track_elbo=T)
)


# Generating X ------
sim_X_sparse <- logisticsusie:::sim_X_sparse

sim_X_dense <- function(...){logisticsusie:::sim_X(...) %>% scale}


#' Simulate dense X with varying correlation structure
make_X_dense_spec <- function(){
  X_spec <- tidyr::tribble(
    ~X_name, ~X_fun, ~X_args, ~X_seed,
    'X_dense_l=1', 'sim_X_dense', list(n=500, p=100, length_scale=1), 3,
    'X_dense_l=10', 'sim_X_dense', list(n=500, p=100, length_scale=10), 4,
    'X_dense_l=50', 'sim_X_dense', list(n=500, p=100, length_scale=50), 5
  )
  return(X_spec)
}

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

# Scoring -------

# compute CSs for a single simulation group
score_cs_row <- function(sim, fit){
  res <- bind_cols(sim, fit) %>%
    rowwise() %>%
    mutate(cs = list(logisticsusie:::get_all_cs2(fit$alpha))) %>%
    mutate(cs = list(logisticsusie:::get_all_coverage(cs, sim$idx))) %>%
    select(cs) %>% ungroup()
}

# compute CSs for multiple simulation groups in a tibble
score_cs <- function(fits){
  scored_cs <- fits %>%
    rowwise() %>%
    mutate(fits = fit$fit) %>%
    select(fit_method, X_name, y_name, sims, fits) %>%
    mutate(score_cs = list(score_cs_row(sims, fits))) %>%
    select(-c(fits))
  return(scored_cs)
}

# compute empirical coverage of each credible set
score_cs_coverage <- function(score_cs){
  # 1 row, 1 simulation
  score_cs2 <- score_cs %>%
    hoist(score_cs, cs='cs') %>%
    unnest_wider(sims) %>%
    unnest_longer(everything())

  # 1 row, 1 credible set
  score_cs3 <- score_cs2 %>%
    unnest_longer(cs) %>%
    hoist(cs, 'covered') %>%
    rowwise() %>% mutate(covered = any(covered))

  # coverage of all credible sets
  # score_cs3 %>%
  #   group_by(fit_method, X_name, y_name) %>%
  #   summarise(empirical_coverage = mean(covered))

  # coverage of non-trivial credible sets
  # coverage_tbl <- score_cs3 %>%
  #   hoist(cs, 'size') %>%
  #   select(fit_method, X_name, y_name, L, covered, size) %>%
  #   filter(size < 10) %>%
  #   group_by(fit_method, X_name, y_name) %>%
  #   summarise(empirical_coverage = mean(covered))

  return(score_cs3)
}

#' extract pips and causal status of each feature
#' return one row per simulation
score_pips <- function(fits){
  score_pip_tbl <- fits %>%
    rowwise() %>%
    mutate(fit = list(fit$fit[[1]])) %>%
    select(X_name, y_name, fit_method, sims, fit) %>%
    unnest_wider(c(sims, fit)) %>%
    unnest_longer(everything()) %>%
    rowwise() %>%
    mutate(
      pip = list(logisticsusie:::get_pip(fit$alpha)),
      causal = list(as.integer(1:length(pip) %in% sim$idx))
    ) %>%
    select(-c(sim, method, fit))
  return(score_pip_tbl)
}

make_pip_roc <- function(pip_scores){
  roc_plot_data <- pip_scores %>%
    group_by(fit_method) %>%
    summarise(
      pip = list(unlist(pip)),
      causal = list(unlist(causal))
    ) %>%
    unnest_longer(everything())
  return(roc_plot_data)
}


# Target definitions --------

# To specify a new simulation
#   1. get the combinations of X, y, and fit specification you want
#   2. call the appropriate scoring functiosn for each target
exec2 <- function(seed, fun, args, ...){
  set.seed(seed)
  args <- c(args, ...)
  if(length(args) < 1){
    rlang::exec(fun)
  } else{
    rlang::exec(fun, !!!args)
  }
}

#' fit model for multiple simulations with a shared X
#' @param X the n x p design matrix
#' @param sims a tibble with m rows, each row is a simulated data set
#' @param method a label for the fit
#' @param fit_fun the fit function name (accepts X and y as first two arguments)
#' @param fit_args optional arguments for the fit function
#' @return and m x 2 tibble with columns 'method' and 'fit'
fit_model <- function(X, sims, method, fit_fun, fit_args = list()){
  fits <- sims %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      method=method,
      fit = list(rlang::exec(fit_fun, X=X, y = sim$y, !!!fit_args))) %>%
    select(method, fit) %>%
    ungroup()
  return(fits)
}

#' Driver function
#' simulates X and y (with random seed) and then fits the model
#' For all combinations of `X_spec`, `y_spec`, and `fit_spec`
#' @param spec a tibble specifying how to simulate X, y, and fit model
fit_from_spec <- function(spec){
  fits <- spec %>%
    rowwise() %>%
    mutate(
      X = list(exec2(X_seed , X_fun, X_args)),
      sims = list(exec2(y_seed, y_fun, y_args, list(X=X))),
      fit = list(fit_model(X, sims, fit_method, fit_fun, fit_args))
    ) %>%
    ungroup() %>%
    nest(  # compact
      y = starts_with('y'),
      X = starts_with('X'),
      fit=starts_with('fit')
    ) %>%
    mutate(  # pull out some identifying info
      fit_method = map_chr(fit, ~pluck(.x, 'fit_method')),
      X_name = map_chr(X, ~pluck(.x, 'X_name')),
      y_name = map_chr(y, ~pluck(.x, 'y_name'))
    )
  return(fits)
}


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
  half_normal_target,
  website
)

