# Scoring -------

# compute CSs for a single simulation group
score_cs_row <- function(sim, fit){
  res <- bind_cols(sim, fit) %>%
    rowwise() %>%
    mutate(cs = list(logisticsusie:::get_all_cs2(fit$alpha))) %>%
    mutate(cs = list(logisticsusie:::get_all_coverage(cs, sim$idx))) %>%
    select(cs) %>% ungroup()
}

# compute CSs for multiple simulation groups in a tibble
score_cs <- function(fits){
  scored_cs <- fits %>%
    rowwise() %>%
    mutate(fits = fit$fit) %>%
    select(fit_method, X_name, y_name, sims, fits) %>%
    mutate(score_cs = list(score_cs_row(sims, fits))) %>%
    select(-c(fits))
  return(scored_cs)
}

# compute empirical coverage of each credible set
score_cs_coverage <- function(score_cs){
  # 1 row, 1 simulation
  score_cs2 <- score_cs %>%
    hoist(score_cs, cs='cs') %>%
    unnest_wider(sims) %>%
    unnest_longer(everything())

  # 1 row, 1 credible set
  score_cs3 <- score_cs2 %>%
    unnest_longer(cs) %>%
    hoist(cs, 'covered') %>%
    rowwise() %>% mutate(covered = any(covered))

  # coverage of all credible sets
  # score_cs3 %>%
  #   group_by(fit_method, X_name, y_name) %>%
  #   summarise(empirical_coverage = mean(covered))

  # coverage of non-trivial credible sets
  # coverage_tbl <- score_cs3 %>%
  #   hoist(cs, 'size') %>%
  #   select(fit_method, X_name, y_name, L, covered, size) %>%
  #   filter(size < 10) %>%
  #   group_by(fit_method, X_name, y_name) %>%
  #   summarise(empirical_coverage = mean(covered))

  return(score_cs3)
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
