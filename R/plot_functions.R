
get_y_range <- function(y){
  res <- quantile(purrr::map_dbl(1:50, ~mean(sample(y, replace = T))), c(0.05, 0.5, 0.95))
  names(res) <- c('ymin', 'y', 'ymax')
  return(as.list(res))
}

make_pip_calibration_plot <- function(pips, facet = 'fit_method'){
  pip_bins <- seq(0, 1.0, by=0.05)

  pip_calibration_data <- pips %>%
    mutate(pip_bin = cut(pip, pip_bins, labels=F)) %>%
    mutate(pip_bin = seq(0.05, 1.0, by=0.05)[pip_bin]) %>%
    select(-pip) %>%
    group_by(across(c(facet, 'pip_bin'))) %>%
    summarise(
      size = n(),
      y_range = list(get_y_range(causal))
    ) %>% 
    unnest_wider(y_range) %>%
    ungroup()
  
  plot <- pip_calibration_data %>%
    ggplot(aes(x=pip_bin, y=y, ymin=ymin, ymax=ymax)) + 
    geom_pointrange(size=0.1) +
    geom_abline(intercept = 0, slope = 1, color='red') +
    facet_wrap(as.formula(paste("~", facet))) + 
    theme_bw()
    
  return(plot)
}

make_coverage_plot <- function(cs, max_size=1e10){
  coverage_plot_data <- cs %>% 
    filter(cs_size < max_size) %>%
    group_by(fit_method, X_name, L) %>%
    summarise(
      size = n(),
      y_range = list(get_y_range(covered))
    ) %>% 
    unnest_wider(y_range) %>%
    ungroup()

  plot <- coverage_plot_data %>%
    ggplot(aes(x=as.factor(L), y=y, ymin=ymin, ymax=ymax, color=as.factor(X_name))) +
    geom_pointrange(position=position_dodge(0.3)) + 
    geom_hline(yintercept = 0.95, size=0.1) + 
    facet_wrap(vars(fit_method)) +
    theme_bw()
  
  return(plot)
}