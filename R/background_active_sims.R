sim_y_ba <- function(X, background, active, L, N=1){
  n <- dim(X)[1]
  p <- dim(X)[2]
  if (length(beta) >= p) {
    stop("length(beta) must be less than number of columns of X")
  }
  idx <- sample(p, L)
  active_idx = matrixStats::rowSums2(X[, idx, drop=F]) > 0
  logits <- rep(background, n)
  logits[which(active_idx)] = active
  p <- sigmoid(logits)
  y <- rbinom(n, N, p)
  data <- list(y = y, logits = logits, N = N, beta = beta,
               background = background, active=active, idx = idx)
  return(data)
}

sim_ba <- function(X, background, active, L){
  # simulate across multiple settings
  reps <- 1:5

  # generate simulations
  sims <- tidyr::crossing(bacgkround=background, active=active, L=L, rep=reps) %>%
    dplyr::mutate(X_name = attributes(X)$name) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sim = list(sim_y_ba(X, background, active, L))) %>%
    ungroup()

  return(sims)
}

make_sim_ba_spec <- function(){
  background = c(-10, -6, -4, -2)
  active = c(0)
  L = c(1)

  ba_spec <- tidyr::crossing(background=background, active=active, L=L) %>%
    rowwise() %>%
    mutate(y_args = list(c_across())) %>%
    ungroup() %>%
    select(y_args) %>%
    mutate(
      y_name = 'background_active',
      y_fun = 'sim_ba',
      y_seed = row_number()
    ) %>%
    mutate(y_sym = rlang::syms(y_fun))
  return(ba_spec)
}

# define models to fit
make_ba_fit_spec <- function(){
  spec <- tidyr::tribble(
    ~fit_method, ~fit_fun, ~fit_args,
    'vb_ser', 'fit_bin_ser', list(),
    'uvb_ser', 'fit_uvb_ser', list(),
    'linear_ser', 'fit_linear_susie', list(L=1),
    'linear_susie_L5', 'fit_linear_susie', list(L=5),
    'binsusie2_L5', 'fit_binsusie_wrapped', list(L=5, estimate_prior_variance=T, prior_variance=1),
    #'ibss2m_L5', 'ibss2m', list(L=5, maxit=50, track_elbo=T)
  ) %>%
    mutate(fit_sym = rlang::syms(fit_fun))
  return(spec)
}

# build simulation spec
make_ba_spec <- function(){
  spec <- tidyr::crossing(
    make_X_binary_spec(),
    make_sim_ba_spec(),
    make_ba_fit_spec()
  )
  return(spec)
}

background_active_target <- list(
  tar_target(ba_spec, make_ba_spec()),
  tar_target(
    ba_fits, fit_from_spec(ba_spec),
    pattern = map(ba_spec)
  ),
  tar_target(ba_pips, score_pips(ba_fits)),
  tar_target(ba_cs, score_cs(ba_fits))
)


