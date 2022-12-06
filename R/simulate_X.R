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
  ) %>%
  mutate(X_sym = rlang::syms(X_fun))
  return(X_spec)
}
