# IBSS ---------
# Define IBSS algorithms with different base SER approximations

fit_ibss_vb <- function(X, y, L){
  logisticsusie::ibss_from_ser(X, y, L=L, ser_function = logisticsusie::fit_bin_ser)
}

fit_bin_ser2 <- purrr::partial(logisticsusie::fit_bin_ser, estimate_prior_variance=T)
fit_ibss_vb2 <- function(X, y, L){
  logisticsusie::ibss_from_ser(X, y, L=L, ser_function = fit_bin_ser2)
}

fit_ibss_vbc <- function(X, y, L){
  logisticsusie::ibss_from_ser(X, y, L=L, ser_function = logisticsusie::fit_bin_ser_corrected)
}

fit_ibss_glm <- function(X, y, L){
  logisticsusie::ibss_from_ser(X, y, L=L, ser_function = logisticsusie::fit_glm_ser)
}

fit_ibss_uvb <- function(X, y, L){
  logisticsusie::ibss_from_ser(X, y, L=L, ser_function = logisticsusie::fit_uvb_ser)
}

fit_uvb2 <- purrr::partial(logisticsusie::fit_uvb_ser, estimate_prior_variance=T)
fit_ibss_uvb2 <- function(X, y, L){
  logisticsusie::ibss_from_ser(X, y, L=L, ser_function = fit_uvb2)
}

fit_ibss_veb <- function(X, y, L){
  logisticsusie::ibss_from_ser(X, y, L=L, ser_function = logisticsusie::fit_veb_ser)
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
  lbfs <- rep(0, fit$hypers$L)
  for(l in 1:fit$hypers$L){
    if(fit$hypers$prior_variance[l] > 0.01){
      elbo_null <- remove_effect(fit, l) %>% 
      {
        fit <- .
        fit$params$xi <- logisticsusie:::update_xi.binsusie(fit)
        fit$params$tau <- logisticsusie:::compute_tau(fit)
        fit
      } %>%
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

  fit$psi <- list(
    mu = logisticsusie:::compute_Xb.binsusie(fit),
    mu2 = logisticsusie:::compute_Xb2.binsusie(fit)
  )

  fit$lbf <- compute_component_lbfs(fit)
  fit$data <- NULL
  fit$xi <- fit$params$xi
  fit$params <- NULL
  fit$elapsed_time <- timer$toc - timer$tic
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