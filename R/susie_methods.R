#' Make fit spec
#' 
#' This function returns a tibble with information for fitting different models
#' All methods are found below in this script
#' All methods take a uniform set of arguments (X, y, ...)
#' fit_args get passed to ...
make_fit_spec <- function(){
  spec <- tidyr::tribble(
    ~fit_method, ~fit_fun, ~fit_args,
    # SERs
    'vb_ser', 'fit_bin_ser', list(),
    'uvb_ser', 'fit_uvb_ser', list(),
    'glm_ser', 'fit_glm_ser', list(),
    'glm_ser2', 'fit_glm_ser2', list(),
    'linear_ser', 'fit_linear_susie', list(L=1),
    # SuSiE
    'linear_susie_L5', 'fit_linear_susie', list(L=5),
    'binsusie_L5', 'fit_binsusie_wrapped', list(L=5, prior_variance=1., estimate_prior_variance=F),
    'binsusie2_L5', 'fit_binsusie_wrapped', list(L=5, prior_variance=1., estimate_prior_variance=T),
    # IBSS2m?
    # 'ibss2m_L5', 'ibss2m_jax', list(L=5),
    'ibss_uvb_L5', 'fit_ibss_uvb', list(L=5),
    'ibss_vb_L5', 'fit_ibss_vb', list(L=5),
    'ibss_glm_L5', 'fit_ibss_glm', list(L=5),
    'ibss_glm2_L5', 'fit_ibss_glm2', list(L=5),
    # GLM augmented,
    # 'glm_ser_aug', 'fit_glm_ser_aug', list(),
    # 'ibss_glm_aug_L5', 'fit_ibss_glm_aug', list(L=5),
    # 'ibss_bayes_L5', 'fit_ibss_bayes', list(L=5, maxit=20),
    # Polynomial susie
    'poly_susie_M2_L5', 'fit_logistic_polysusie', list(L=5, left=-5, right=5, M=2, max_iter=100),
    'poly_susie_M6_L5', 'fit_logistic_polysusie', list(L=5, left=-5, right=5, M=6, max_iter=100),
    'poly_susie_M10_L5', 'fit_logistic_polysusie', list(L=5, left=-5, right=5, M=10, max_iter=100)
    ) %>%
    mutate(fit_sym = rlang::syms(fit_fun))
  return(spec)
}

# SERS ----------
fit_bin_ser <- logisticsusie::binser
fit_uvb_ser <- logisticsusie::uvbser
fit_glm_ser <- logisticsusie::fit_glm_ser
fit_glm_ser2 <- logisticsusie::fit_glm_ser2


# IBSS ---------
# Define IBSS algorithms with different base SER approximations

fit_ibss_vb <- function(X, y, L, ...){
  logisticsusie::ibss_from_ser(X, y, L=L, ser_function = logisticsusie::binser, ...)
}

fit_ibss_glm <- function(X, y, L, ...){
  logisticsusie::ibss_from_ser(X, y, L=L, ser_function = logisticsusie::fit_glm_ser, ...)
}

# uses corrected ABF
fit_ibss_glm2 <- function(X, y, L, ...){
  logisticsusie::ibss_from_ser(X, y, L=L, ser_function = logisticsusie::fit_glm_ser2, ...)
}

fit_ibss_uvb <- function(X, y, L, ...){
  logisticsusie::ibss_from_ser(X, y, L=L, ser_function = logisticsusie::uvbser, ...)
}

# Polynomial susie---------
fit_logistic_polysusie <- function(X, y, L, ...){
  polysusie::logistic_polysusie(X, y, L, ...)
}


# binsusie -----------

#' remove the lth effect from a susie fit
#' useful for computing the ELBO of the alternative nodel where b_l=0
remove_effect <- function(fit, l=1){
  fit2 <- fit
  fit2$params$alpha <- fit2$params$alpha[-l,]
  fit2$params$mu <- fit2$params$mu[-l,]
  fit2$params$var <- fit2$params$var[-l,]
  fit2$params$delta <- matrix(fit2$params$delta[-l,], ncol = 1)
  
  fit2$hypers$L <- fit2$hypers$L - 1
  fit2$hypers$pi <- fit2$hypers$pi[-l,]
  fit2$hypers$prior_mean <- fit2$hypers$prior_mean[-l]
  fit2$hypers$prior_variance <- fit2$hypers$prior_variance[-l]
  return(fit2)
}

#' compute ElBO - ELBO_{-l}
#' where ELBO_{-l} is the model with b_l=0
compute_component_lbfs <- function(fit){
  lbfs <- rep(0, fit$L)
  for(l in 1:fit$L){
    if(fit$hypers$prior_variance[l] > 0.01){
      elbo_null <- fit %>%
        remove_effect(l) %>% 
        logisticsusie:::fit.binsusie(fit_prior_variance = F, fit_alpha = F) %>%
        logisticsusie:::compute_elbo.binsusie()
      lbfs[l] <- tail(fit$elbo, 1) - elbo_null
    }
  }
  return(lbfs)
}

fit_binsusie_wrapped <- function(X, y, L, ...){
  tictoc::tic()
  fit <- logisticsusie::binsusie(X, y, L=L, ...)
  timer <- tictoc::toc()
  
  fit$elapsed_time <- timer$toc - timer$tic

  data <- logisticsusie::binsusie_prep_data(X, y, N=rep(1, nrow(X)), Z=NULL)
  fit$psi <- list(
    mu = logisticsusie:::compute_Xb.binsusie(fit, data),
    mu2 = logisticsusie:::compute_Xb2.binsusie(fit, data)
  )

  fit$lbf <- rep(0, fit$L) #compute_component_lbfs(fit)
  fit$data <- NULL
  fit$xi <- fit$params$xi
  fit$params <- NULL
  return(fit)
}

# Linear SuSiE-----

fit_linear_susie <- function(X, y, L, ...){
  fit <- susieR::susie(X, y, L, scaled_prior_variance = 1.0, standardize = F)
  fit$mu <- drop(fit$mu)
  fit$var <- drop(fit$mu2) - fit$mu^2
  fit$alpha <- drop(fit$alpha)
  return(fit)
}

# Augmented GLM ------------
# add data to stabilize glm estimation

#' add two data points, ensure
augment_data <- function(X, y){
  p <- dim(X)[2]

  newX <- rbind(rep(0, p), rep(0.1, p))
  Xaug <- rbind(X, newX)
  yaug <- c(y, c(0, 1))
  return(list(y=yaug, X=Xaug))
}

fit_glm_ser_aug <- function(X, y, ...){
  aug <- augment_data(X, y)
  return(logisticsusie::fit_glm_ser(aug$X, aug$y, ...))
}

fit_ibss_glm_aug <- function(X, y, L, ...){
  aug <- augment_data(X, y)
  logisticsusie::ibss_from_ser(aug$X, aug$y, L=L, ser_function = logisticsusie::fit_glm_ser, ...)
}