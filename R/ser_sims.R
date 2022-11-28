sigmoid <- function(x){1/(1 + exp(-x))}

simulate_y_ser <- function(X, b0, b, idx=NULL){
  X <- scale(X)
  p <- ncol(X)
  # select a random index if not provied
  if(is.null(idx)){
    idx <- sample(1:p, 1)
  }
  logits <- b0 + X[, idx] * b
  y <- rbinom(length(logits), 1, sigmoid(logits))
  if(sum(y) == 0){
    y[which.max(y)] <- 1  # at least one case
  }
  return(list(y=y, logits=logits, b0=b0, b=b, idx=idx))
}

# simulate X
get_dense_Xs <- function(){
  lengthscales <- c(1, 10, 100, 200)
  Xs <- map(lengthscales, ~(logisticsusie:::sim_X(n=200, p=50, length_scale=.x) %>% scale))
  names(Xs) <- paste0('lengthscale=', lengthscales)
  return(Xs)
}

simulate_ser <- function(Xs){
  # simulation parameters
  b0 <- c(-2, -1, -0.5, 0)
  b <- c(0.2, 0.4, 0.6, 0.8, 1.)
  reps <- 1:10
  spec <- tidyr::crossing(b0, b, rep=reps, X=names(Xs))

  # simulate y
  sim <- spec %>%
    rowwise() %>%
    mutate(
      sim = list(simulate_y(Xs[[X]], b0, b))
    )
  return(sim)
}

fit_sers <- function(sim, Xs){
  # fit SERs
  fit <- sim %>%
    rowwise() %>%
    mutate(
      vb_ser = list(fit_bin_ser(Xs[[X]], sim$y, prior_variance=1)),
      uvb_ser = list(fit_uvb_ser(Xs[[X]], sim$y, prior_variance = 1)),
      veb_ser = list(fit_veb_ser(Xs[[X]], sim$y, prior_variance = 1)),
      quad_ser = list(fit_quad_ser(Xs[[X]], sim$y, prior_variance = 1))
    ) %>%
    mutate(
      vb_ser_corrected = list(correct_bin_ser(Xs[[X]], sim$y, 0, fit = vb_ser)),
    )
  return(fit)
}


