# Yusha example -------
# this is an example where binsusie essentially fits the null model.
# in the below analysis we fit various SERs and SuSiEs
# to see if our new strategies improve the situation

# analysis function
# take yushas example, take the n_top marginal enrichments and fit logistic ibss

#' Function returns list with X and y
yusha_example_data <- function(n_top=1e10){
  example <- readRDS('data/yusha_sc_tumor/pdac_example.rds')
  bindata <- with(example, gseasusie::prep_binary_data(genesets, data, thresh = 0.01))
  gs_names <- colnames(bindata$X)
  n_gene_sets <- dim(bindata$X)[2]

  ora <- with(bindata, gseasusie::fit_ora(X, y))
  top_idx <- head(order(ora$pFishersExact), n_top)
  bindata_sub <- bindata
  bindata_sub$X <- bindata_sub$X[, top_idx]

  return(list(bindata = bindata_sub, ora=ora))
}

ibss2m_jax_nosparse <- function(X, y, ...){
  logisticsusie::ibss2m_jax(as.matrix(X), y, ...)
}

#' Specify which models to fit
make_yusha_fit_spec <- function(){
  spec <- tidyr::tribble(
    ~fit_method, ~fit_fun, ~fit_args,
    # SERs
    'vb_ser', 'fit_bin_ser', list(),
    'uvb_ser', 'fit_uvb_ser', list(),
    'glm_ser', 'fit_glm_ser', list(),
    'linear_ser', 'fit_linear_susie', list(L=1),
    'bayes_ser', 'fit_bayes_ser', list(),
    # SuSiE
    'linear_susie_L10', 'fit_linear_susie', list(L=10),
    'binsusie_L10', 'fit_binsusie_wrapped', list(L=10, prior_variance=1, estimate_prior_variance=F),
    #'binsusie2_L10', 'fit_binsusie_wrapped', list(L=10, estimate_prior_variance=T, prior_variance=1),
    #'ibss2m_L10', 'ibss2m_jax_nosparse', list(L=10),
    'ibss_uvb_L10', 'fit_ibss_uvb', list(L=10),
    #'ibss_vb_L10', 'fit_ibss_vb', list(L=10)
  ) %>%
  mutate(fit_sym = rlang::syms(fit_fun))
  return(spec)
}

#' fit model
fit_model2 <- function(X, y, method, fit_fun, fit_args = list()){
    fit <- rlang::exec(fit_fun, X=X, y=y, !!!fit_args)
    return(fit)
}

yusha_example_targets <- list(
    tar_target(yusha_data, yusha_example_data()$bindata),
    tar_target(yusha_fit_spec, make_yusha_fit_spec()),
    tar_target(yusha_fits,
      yusha_fit_spec %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            fit = list(fit_model2(
              yusha_data$X, yusha_data$y, fit_method, fit_fun, fit_args))
        ),
        pattern = map(yusha_fit_spec)
    ),
    tar_target(yusha_data500, yusha_example_data(500)$bindata),
    tar_target(yusha_fits500,
      yusha_fit_spec %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            fit = list(with(yusha_data500,
              fit_model2(X, y, fit_method, fit_fun, fit_args)))
        ),
        pattern = map(yusha_fit_spec)
    )
)
