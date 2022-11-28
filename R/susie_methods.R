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

# IBSS2 ---------
# Define IBSS algorithms with different base SER approximations


