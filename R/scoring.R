# Scoring -------

# Score CSs-------

#' Check if a CS covers a particular index
get_coverage <- function(cs, idx) {
  cs$covered <- idx %in% cs$cs
  cs$which_covered <- intersect(idx, cs$cs)
  names(cs$covered) <- idx
  return(cs)
}

#' Check if set of indices `idx` in a list of credible sets
get_all_coverage <- function(css, idx) {
  purrr::map(css, ~ get_coverage(.x, idx))
}

#' compute CS and annotate with coverage info
make_scored_cs <- function(alpha, idx){
  cs <- logisticsusie::compute_cs(alpha) %>%
    { 
      # special case for SER
      if(is.null(nrow(alpha))){
        list(L1 = .)
      } else{
        .
      }
    } %>%
    get_all_coverage(idx)
  return(cs)
}

#' extract log Bayes factors from each model fit
get_lbfs <- function(fit_method, fit){
  if(stringr::str_detect(fit_method, 'ser')){
    lbf <- fit$lbf_model
  } else{
    lbf <- fit$lbf
  }
  return(lbf)
}

score_cs <- function(fits){
  cs1 <- fits %>%
    select(X_name, y_name, fit_method, sims, fit) %>%
    unnest_wider(sims) %>%
    hoist(fit,  fitted = list('fit')) %>%
    unnest_longer(-c(X_name, y_name, fit_method)) %>%
    rowwise() %>%
    mutate(
      alpha = list(fitted$alpha),
      idx = list(sim$idx)
    ) %>%
    ungroup() %>%
    mutate(hash = purrr::map_chr(sim, rlang::hash)) %>%
    select(-c(sim, fitted, fit))

  cs2 <- cs1 %>%
    rowwise() %>%
    mutate(
      cs = list(make_scored_cs(alpha, idx))
    ) %>%
    select(-c(alpha))
  return(cs2)
}

# Score PIP ------------

compute_pip <- function(alpha){
  if(is.null(nrow(alpha))){
    return(alpha)
  } else{
    p <- ncol(alpha)
    pip <- purrr::map_dbl(1:p, ~1 - prod(1 - alpha[, .x]))
    return(pip)
  }
}

#' Score PIPs
#' 
#' Function computes pips and extracts true causal status of each variable from simulaiton
#' to be used e.g. for power/FDR plots
#' 
#' @return a tibble each row is a simulation, 
#'   columns included identitying information and list-columns 
#'   `pip` and `causal`
score_pips <- function(fits){
  score_pip_tbl <- fits %>%
    select(X_name, y_name, fit_method, sims, fit) %>%
    unnest_wider(sims) %>%
    hoist(fit,  fitted = list('fit')) %>%
    unnest_longer(-c(X_name, y_name, fit_method)) %>%
    rowwise() %>%
    mutate(
      pip = list(compute_pip(fitted$alpha)),
      causal = list(as.integer(1:length(pip) %in% sim$idx))
    ) %>%
    select(-c(sim, fitted, fit))
  return(score_pip_tbl)
}

make_pip_roc <- function(pip_scores){
  roc_plot_data <- pip_scores %>%
    group_by(fit_method) %>%
    summarise(
      pip = list(unlist(pip)),
      causal = list(unlist(causal))
    ) %>%
    unnest_longer(everything())
  return(roc_plot_data)
}
