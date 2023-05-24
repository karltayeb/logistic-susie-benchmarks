
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
fit_model_to_sims <- function(X, sims, fit_method, fit_fun, fit_args = list()){
  fits <- sims %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      method=fit_method,
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
      fit = list(fit_model_to_sims(X, sims, fit_method, fit_fun, fit_args)) 
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
