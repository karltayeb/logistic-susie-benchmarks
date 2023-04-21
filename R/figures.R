# code and targets for making figures

library(ggplot2)
library(dplyr)
library(tidyr)
#---- PIP Calibration


#pips <- tar_read(half_normal_pips)

pips <- tar_read(constant_pips)

pips %>% 
  filter_ser() %>%
  filter(L==1) %>%
  make_pip_calibration_plot() + facet_wrap(vars(fit_method))


pips %>% 
  filter_ser() %>%
  filter(L>1) %>%
  make_pip_calibration_plot() + facet_wrap(vars(fit_method))


susie_methods <- c(
    'binsusie_L5',
    'ibss_vb_L5',
    'ibss_glm_aug_L5',
    'ibss_uvb_L5',
    'ibss2m_L5',
    'linear_susie_L5'
)

pips %>% 
  filter_susie() %>%
  filter(fit_method %in% susie_methods) %>%
  filter(L==1) %>%
  make_pip_calibration_plot() + facet_wrap(vars(fit_method))


pips %>% 
  filter_susie() %>%
  filter(fit_method %in% susie_methods) %>%
  filter(L>1) %>%
  make_pip_calibration_plot() + facet_wrap(vars(fit_method))



pips %>%
  filter_ser() %>%
  filter(L==1) %>%
  make_fdp_plot(max_fdp=0.25)

pips %>%
  filter_ser() %>%
  filter(L>1) %>%
  make_fdp_plot(max_fdp=0.25)

pips %>%
  filter_susie() %>%
  filter(L==1) %>%
  filter(fit_method %in% susie_methods) %>%
  make_fdp_plot(max_fdp=0.25)

pips %>%
  filter_susie() %>%
  filter(L>1) %>%
  filter(fit_method %in% susie_methods) %>%
  make_fdp_plot(max_fdp=0.25)


cs <- tar_read(constant_cs)
X_order <- unique(pips$X_name)[c(3,1,2)]
cs$X_name <- factor(cs$X_name, levels = X_order) 

cs %>%
  filter_ser() %>%
  make_cs_coverage_plot() +
  facet_wrap(vars(X_name))

cs %>%
  filter_susie() %>%
  filter(fit_method %in% susie_methods) %>%
  make_cs_coverage_plot() +
  facet_wrap(vars(X_name))

cs %>%
  filter_ser() %>%
  make_coverage_by_cs_size() +
  facet_wrap(vars(X_name))

cs %>%
  filter_susie() %>%
  filter(fit_method %in% susie_methods) %>%
  make_coverage_by_cs_size() +
  facet_wrap(vars(X_name))

# Polysusie comparisons------
methods <- c(
  'binsusie_L5',
  #'ibss_vb_L5',
  'ibss_glm_aug_L5',
  'ibss_uvb_L5',
  #'ibss2m_l5',
  'linear_susie_L5',
  'poly_susie_M2_L5',
  'poly_susie_M6_L5',
  'poly_susie_M10_L5'
)

pips %>%
  filter_susie() %>%
  filter(fit_method %in% methods) %>%
  filter(L>1) %>%
  make_fdp_plot(max_fdp=0.3)
