# Simulations where there is a background rate and and active rate
# Genes in the gene set are observed at the active rate
# all other genes are observed at the background rate
# for multiple active gene sets, just use the union of these gene sets

sim_y_constant <- function(X, background, active, L, N=1){
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

sim_constant <- function(X, background, active, L){
  # simulate across multiple settings
  reps <- 1:20

  # generate simulations
  sims <- tidyr::crossing(constantcgkround=background, active=active, L=L, rep=reps) %>%
    dplyr::mutate(X_name = attributes(X)$name) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sim = list(sim_y_constant(X, background, active, L))) %>%
    ungroup()

  return(sims)
}

make_constant_sim_spec <- function(){
  background = c(-10)
  active = c(0)
  L = c(1)

  constant_spec <- tidyr::crossing(background=background, active=active, L=L) %>%
    rowwise() %>%
    mutate(y_args = list(c_across())) %>%
    ungroup() %>%
    select(y_args) %>%
    mutate(
      y_name = 'background_active',
      y_fun = 'sim_constant',
      y_seed = row_number()
    ) %>%
    mutate(y_sym = rlang::syms(y_fun))
  return(constant_spec)
}

# define models to fit
make_constant_fit_spec <- function(){
  spec <- tidyr::tribble(
    ~fit_method, ~fit_fun, ~fit_args,
    'vb_ser', 'fit_bin_ser', list(),
    'uvb_ser', 'fit_uvb_ser', list(),
    'linear_ser', 'fit_linear_susie', list(L=1),
    'linear_susie_L5', 'fit_linear_susie', list(L=5),
    'binsusie2_L5', 'fit_binsusie_wrapped', list(L=5, estimate_prior_variance=T, prior_variance=1),
    'ibss2m_L5', 'ibss2m', list(L=5, maxit=50, track_elbo=T)
  ) %>%
    mutate(fit_sym = rlang::syms(fit_fun))
  return(spec)
}

# build simulation spec
make_constant_spec <- function(){
  spec <- tidyr::crossing(
    make_X_binary_spec(),
    make_constant_sim_spec(),
    make_constant_fit_spec()
  )
  return(spec)
}

constant_sim_target <- list(
  tar_target(constant_spec, make_constant_spec()),
  tar_target(
    constant_fits, fit_from_spec(constant_spec),
    pattern = map(constant_spec)
  ),
  tar_target(constant_pips, score_pips(constant_fits)),
  tar_target(constant_cs, score_cs(constant_fits))
)