# Generating X ------
sim_X_sparse <- function(...){as.matrix(logisticsusie:::sim_X_sparse(...))}

#' Simulate dense X with varying correlation structure
make_X_binary_spec <- function(){
  X_spec <- tidyr::tribble(
    ~X_name, ~X_fun, ~X_args, ~X_seed,
    'X_bin_strong', 'sim_X_sparse', list(n=500, p=100, pi1=0.1, p_flip=0.02), 3,
    'X_bin_medium', 'sim_X_sparse', list(n=500, p=100, pi1=0.1, p_flip=0.1), 4,
    'X_bin_weak', 'sim_X_sparse', list(n=500, p=100, pi1=0.1, p_flip=0.5), 5
  ) %>%
    mutate(X_sym = rlang::syms(X_fun))
  return(X_spec)
}

sim_X_dense <- function(...){logisticsusie:::sim_X(...) %>% scale}

#' Simulate dense X with varying correlation structure
make_X_dense_spec <- function(){
  X_spec <- tidyr::tribble(
    ~X_name, ~X_fun, ~X_args, ~X_seed,
    'X_dense_l=1', 'sim_X_dense', list(n=500, p=100, length_scale=1), 3,
    'X_dense_l=10', 'sim_X_dense', list(n=500, p=100, length_scale=10), 4,
    'X_dense_l=50', 'sim_X_dense', list(n=500, p=100, length_scale=50), 5
  ) %>%
  mutate(X_sym = rlang::syms(X_fun))
  return(X_spec)
}

make_X_c2_random_1k <- function(min_genes = 10, dense = F){
  set.seed(1)
  X <- gseasusie::load_msigdb_geneset_x('C2')$X
  X <- X[, BiocGenerics::colSums(X) > min_genes]
  X <- X[, sample(dim(X)[2], 1000)]
  X <- X[BiocGenerics::rowSums(X) > 0,]

  if(dense){
    X <- as.matrix(X)
  }
  return(X)
}

make_X_yusha <- function(){
  data <- yusha_example_data()
  X <- as.matrix(data$bindata$X)
  return(X)
}

#' generate X from real gene sets
#' includes Yusha's example, a random 1k sample of Msigdb C2
make_realX_spec <- function(){
  X_spec <- tidyr::tribble(
    ~X_name, ~X_fun, ~X_args, ~X_seed,
    'X_yusha', 'make_X_yusha', list(), 3
    # 'X_c2_1k', 'make_X_c2_random_1k', list(min_genes=10, dense=T), 4
  ) %>%
  mutate(X_sym = rlang::syms(X_fun))
  return(X_spec)
}

make_X_targets <- function(){
  list(
    tar_target(X_c2_random_1k, make_X_c2_random_1k())
  )
}