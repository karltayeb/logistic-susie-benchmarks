
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


# PIPs Power vs FDR
compute_power_fdr_thresh <- function(pip, causal){
  n_causal <- sum(causal)

  sorter <- order(pip, decreasing = T)
  sorted_pips <- pip[sorter]
  sorted_causal <- causal[sorter]
  n <- 1:length(sorted_causal)

  power <- cumsum(sorted_causal) / n_causal
  fdr <- (n - cumsum(sorted_causal)) / n
  thresh <- sorted_pips

  # prune
  include <- logical(length(thresh))
  include[1] <- T
  last <- 1
  for(i in 2:length(thresh)){
    if(((fdr[i] - fdr[last]) > 0.005) | ((power[i] - power[last]) > 0.005)){
      include[i] <- T
      last <- i
    }
  }

  power <- power[include]
  fdr <- fdr[include]
  thresh <- thresh[include]
  return(list(power=power, fdr=fdr, thresh=thresh))
}

#' Plot power vs FDP for each fit_metho
#' @params pips a data frame with columns fit_method, pip, causal
#' @params colors a list matching fit_methdods to colors
#' @params max_fdr the maximum FDR level to show in the plot (ie xlim)
make_power_fdr_plot <- function(pips, colors, max_fdr=0.25){
  methods <- unique(pips$fit_method)
  colors <- colors[methods]

  plot_data <- pips %>%
    group_by(fit_method) %>%
    summarise(pft = list(compute_power_fdr_thresh(pip, causal))) %>%
    unnest_wider(pft) %>%
    unnest_longer(c(power, fdr, thresh))

  max_power <- plot_data %>%
    filter(fdr <= max_fdr) %>%
    {max(.$power)}

  plot <- plot_data %>%
    ggplot(aes(x=fdr, y=power, color=fit_method)) + 
    geom_path() + 
    theme_bw() +
    scale_color_manual(values=colors) +
    xlim(0, max_fdr) + ylim(0, max_power)
  return(plot)
}

