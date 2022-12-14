# Scoring -------

make_scored_cs <- function(alpha, idx){
  cs <- logisticsusie:::get_all_cs2(alpha) %>%
    logisticsusie:::get_all_coverage(idx)
  return(cs)
}

get_lbfs <- function(fit_method, fit){
  if(stringr::str_detect(fit_method, 'ser')){
    lbf <- fit$lbf_model
  } else{
    lbf <- fit$lbf
  }
  return(lbf)
}

score_cs <- function(fits){
  # 1 row 1 simulation
  tmp <- fits %>%
    rowwise() %>%
    mutate(fits = fit$fit) %>%
    select(fit_method, X_name, y_name, sims, fits) %>%
    unnest_wider(c(sims, fits)) %>%
    unnest_longer(c(sim, fit))

  # score each simulation
  tmp2 <- tmp %>%
    rowwise() %>%
    mutate(
      cs = list(make_scored_cs(fit$alpha, sim$idx)),
      lbf = list(get_lbfs(fit_method, fit))
    )

  return(tmp2)
}

#' extract pips and causal status of each feature
#' return one row per simulation
score_pips <- function(fits){
  score_pip_tbl <- fits %>%
    rowwise() %>%
    mutate(fit = list(fit$fit[[1]])) %>%
    select(X_name, y_name, fit_method, sims, fit) %>%
    unnest_wider(c(sims, fit)) %>%
    unnest_longer(everything()) %>%
    rowwise() %>%
    mutate(
      pip = list(logisticsusie:::get_pip(fit$alpha)),
      causal = list(as.integer(1:length(pip) %in% sim$idx))
    ) %>%
    select(-c(sim, method, fit))
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
