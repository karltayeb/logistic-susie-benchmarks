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
make_half_normal_fit_spec <- function(){
  spec <- tidyr::tribble(
    ~fit_method, ~fit_fun, ~fit_args,
    # SERs
    'vb_ser', 'fit_bin_ser', list(),
    'uvb_ser', 'fit_uvb_ser', list(),
    'glm_ser', 'fit_glm_ser', list(),
    'bayes_ser', 'fit_bayes_ser', list(),
    'linear_ser', 'fit_linear_susie', list(L=1),
    'quad_ser', 'fit_quad_ser', list(n=2^8),
    # SuSiE
    'linear_susie_L5', 'fit_linear_susie', list(L=5),
    'binsusie_L5', 'fit_binsusie_wrapped', list(L=5, prior_variance=1, estimate_prior_variance=F),
    'binsusie2_L5', 'fit_binsusie_wrapped', list(L=5, estimate_prior_variance=T, prior_variance=1),
    'ibss2m_L5', 'ibss2m_jax', list(L=5),
    'ibss_uvb_L5', 'fit_ibss_uvb', list(L=5),
    'ibss_vb_L5', 'fit_ibss_vb', list(L=5),
    'ibss_glm_L5', 'fit_ibss_glm', list(L=5),
    # GLM augmented,
    'glm_ser_aug', 'fit_glm_ser_aug', list(),
    'ibss_glm_aug_L5', 'fit_ibss_glm_aug', list(L=5),
    # Polynomial susie
    'poly_susie_M2_L5', 'fit_logistic_polysusie', list(L=5, left=-5, right=5, M=2),
    'poly_susie_M6_L5', 'fit_logistic_polysusie', list(L=5, left=-5, right=5, M=6),
    'poly_susie_M10_L5', 'fit_logistic_polysusie', list(L=5, left=-5, right=5, M=10)
  ) %>%
  mutate(fit_sym = rlang::syms(fit_fun))
  return(spec)
}

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
