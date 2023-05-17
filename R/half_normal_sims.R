# Half normal simulations----

single_half_normal_sim <- function(X, beta0, beta_sigma, L){
  # simulate across multiple settings
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

# define simulation parameters
make_half_normal_sim_spec <- function(){
  beta0 <- c(-2)
  beta_sigma <- c(0.2, 0.4, 0.6, 0.8, 1.0, 2.0) / sqrt(2/pi)
  L <- c(1, 3, 5)
  half_normal_spec <- tidyr::crossing(beta0=beta0, beta_sigma=beta_sigma, L=L) %>%
    rowwise() %>%
    mutate(y_args = list(c_across())) %>%
    ungroup() %>%
    select(y_args) %>%
    mutate(
      y_name = 'half_normal',
      y_fun = 'single_half_normal_sim',
      y_seed = row_number()
    ) %>%
    mutate(y_sym = rlang::syms(y_fun))
  return(half_normal_spec)
}

# define models to fit
source('R/susie_methods.R')
make_half_normal_fit_spec <- make_fit_spec

# build simulation spec
make_half_normal_spec <- function(){
  spec <- tidyr::crossing(
    make_X_dense_spec(),
    make_half_normal_sim_spec(),
    make_half_normal_fit_spec()
  )
  return(spec)
}

half_normal_target <- list(
  tar_target(half_normal_spec, make_half_normal_spec()),
  tar_target(
    half_normal_fits, fit_from_spec(half_normal_spec),
    pattern = map(half_normal_spec)
  ),
  tar_target(half_normal_pips, score_pips(half_normal_fits), pattern=map(half_normal_fits)),
  tar_target(half_normal_cs, score_cs(half_normal_fits), pattern=map(half_normal_fits))
)
