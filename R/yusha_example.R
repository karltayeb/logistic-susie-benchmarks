# Yusha example -------
# this is an example where binsusie essentially fits the null model.
# in the below analysis we fit various SERs and SuSiEs
# to see if our new strategies improve the situation

#' Load background genes for pdac example
load_pdac_example <- function(){

  # load background genes
  genes_in_study <- read.csv2('data/yusha_sc_tumor/pdac/gene_list.txt', header = F)$V1
  hs <- org.Hs.eg.db::org.Hs.eg.db


  # map to entrez IDs
  idmap <- AnnotationDbi::select(hs, keys = genes_in_study,
                                columns = c('SYMBOL', 'ENTREZID'),
                                keytype = 'SYMBOL')
  genes_in_study_entrez <- idmap$ENTREZID[!is.na(idmap$ENTREZID)]


  # get file names of gene lists
  files <- list.files('data/yusha_sc_tumor/pdac/')
  keep <- purrr::map_lgl(files, ~stringr::str_detect(.x, 'factor'))
  files <- files[keep]

  # load gene lists, convert to ENTREZ
  gene_lists <- list()
  gene_lists_symbol <- list()
  for(file in files){
    name <- strsplit(file, '\\.')[[1]][1]
    path <- paste0('data/yusha_sc_tumor/pdac/', file)
    gene_list <- read.csv2(path, header=F)$V1
    gene_list_entrez <- idmap %>%
      dplyr::filter(!is.na(ENTREZID)) %>%
      dplyr::filter(SYMBOL %in% gene_list) %>% 
      {.$ENTREZID}
    gene_lists[[name]] <- gene_list_entrez
    gene_lists_symbol[[name]] <- gene_list
  }

  data <- list(background = genes_in_study_entrez, gene_lists = gene_lists)
  return(data)
}

# for static branching over factors/gene lists
make_pdac_tbl <- function(){
  # get file names of gene lists
  files <- list.files('data/yusha_sc_tumor/pdac/')
  keep <- purrr::map_lgl(files, ~stringr::str_detect(.x, 'factor'))
  files <- files[keep]
  names <- purrr::map_chr(files, ~strsplit(.x, '\\.')[[1]][1])

  pdac <- dplyr::as_tibble(data.frame(study='pdac', factor = names))
  return(pdac)
}

geneSet2X <- function(geneSet){
  # construct sparse matrix
  genes <- unique(geneSet$gene)
  sets <- unique(geneSet$geneSet)
  
  geneSet$row <- match(geneSet$gene, genes)
  geneSet$col <- match(geneSet$geneSet, sets)
  
  X <- Matrix::sparseMatrix(
    i = geneSet$row, 
    j = geneSet$col,
    x = rep(1, nrow(geneSet)),
    dimnames = list(genes, sets)
  )
}

load_msigdb_X <- function(background){
  msigdb <- gseasusie::load_gene_sets(c('c2', 'h')) # load mSigDB "C2" and "H"
  geneSet <- 
    purrr::map_dfr(msigdb, ~purrr::pluck(.x, "geneSet")) %>% # concat "C2" and "H"
    dplyr::select(geneSet, gene) %>%
    dplyr::distinct() %>% 
    dplyr::filter(gene %in% background) %>% # subset to genes in study
    dplyr::group_by(geneSet) %>%
    dplyr::mutate(n = length(gene)) %>%
    dplyr::filter(n >= 10) %>% # filter to gene sets of a minimum size
    dplyr::select(-n)
  X <- geneSet2X(geneSet) # make gene set matrix
  return(X)
}

fit_lasso <- function(X, gene_list){
  y <- as.integer(rownames(X) %in% gene_list)
  lasso_fit <- glmnet::cv.glmnet(X, y, family = "binomial")
  return(lasso_fit)
}

fit_lasso_all_gene_lists <- function(X, gene_lists){
  lasso_fits <- purrr::map(gene_lists, ~fit_lasso(X, .x))
  names(lasso_fits) <- names(gene_lists)
  return(lasso_fits)
}

simulate_from_lasso <- function(X, lasso_fit, reps){
  coef_1se <- coef(lasso_fit, s = 'lambda.1se')[,1]
  logit_1se <- (X %*% tail(coef_1se, -1))[,1] + coef_1se[1]
  ysim_1se <- purrr::map(1:reps, ~rbinom(length(logit_1se), 1, 1/(1 + exp(-logit_1se))))
}

simulate_from_lasso_map <- function(X, lasso_fits, reps){
  simulations <- purrr::map(lasso_fits, ~simulate_from_lasso(X, .x, reps))
  names(simulations) <- names(lasso_fits)
  return(simulations)
}

fit_ibss <- function(gene_list, X, L, maxit=20){
  # make response variable y_i = 1[gene i in gene list]
  y <- as.integer(rownames(X) %in% gene_list)
  
  # define ser function
  ser_fun <- purrr::partial(
    logisticsusie::fit_glm_ser2,
    laplace=T,
    estimate_prior_variance = T)
  
  # fit ibss
  ibss_fit <- logisticsusie::ibss_from_ser(X, y, 
                                           L=L, 
                                           ser_function = ser_fun,
                                           maxit = maxit)
  
  # compute CS, add purity
  ibss_fit$cs <- logisticsusie::compute_cs(ibss_fit$alpha)
  for(l in names(ibss_fit$cs)){
    ibss_fit$cs[[l]]$purity <- logisticsusie::compute_purity(ibss_fit$cs[[l]]$cs, X)
  }

  # collect 
  ibss_fit$prior_variance <- purrr::map_dbl(ibss_fit$fits, ~.x$prior_variance)
  ibss_fit$lbf_ser <- purrr::map_dbl(ibss_fit$fits, ~.x$lbf_model)
  ibss_fit$cs_size <- purrr::map_dbl(ibss_fit$cs, ~length(.x$cs))

  return(ibss_fit)
}

fit_ibss_select_L <- function(gene_list, X){
  L <- 8
  L_is_big_enough <- F
  tictoc::tic()
  while(!L_is_big_enough){
    ibss_fit <- fit_ibss(gene_list, X, L=L, maxit=20)
    if(min(ibss_fit$lbf_ser) < 0){
      L_is_big_enough <- T
    } else if(max(ibss_fit$cs_size) > (ncol(X) / 2)){
      L_is_big_enough <- T
    } else if(min(ibss_fit$prior_variance) < 0.1){
      L_is_big_enough <- T
    } else{
      message('L was too small-- doubling')
      L <- L * 2
    }
  }
  tictoc::toc()
  return(ibss_fit)
}

#' Fit IBSS marginal only
#' 
#' fits ibss using only gene sets that are marginally enriched
#' sets L to be the smallest multiple of 10 thats at least as large
#' as the number of effects returned by lasso
fit_ibss_marginal_only <- function(gene_list, X, maxit=50){
  y <- as.integer(rownames(X) %in% gene_list)
  ora <- gseasusie::fit_ora(X, y)
  marginally_enriched <- ora %>%
    dplyr::mutate(pAdj = p.adjust(pFishersExact, method='BH')) %>%
    dplyr::filter(pAdj < 0.05) %>%
    { .$geneSet }
  print(length(marginally_enriched))

  lasso_fit <- fit_lasso(X, gene_list)
  lasso_L <- sum(coef(lasso_fit, s='lambda.1se') !=0) - 1
  L <- max(10, ceiling(lasso_L / 10) * 10)
  print(L)
  ibss_fit <- fit_ibss(gene_list, X[, marginally_enriched], L = L, maxit=maxit)
  ibss_fit$var_names <- marginally_enriched
  return(ibss_fit)
}

fit_ibss_lasso_init <- function(gene_list, X, L=20, maxit=100){
  y <- as.integer(rownames(X) %in% gene_list)
  fit <- logisticsusie::gibss_lasso_init(X, y, L=L, family='binomial', maxit=maxit)
  fit$var_names <- colnames(X)
}

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


pdac_tbl <- make_pdac_tbl()
pdac_targets <- list(
    tar_target(pdac_data, load_pdac_example()),
    tar_target(pdac_msigdb_X, load_msigdb_X(pdac_data$background)),
    tar_map(pdac_tbl, 
      tar_target(ibss_fit_marginal_only, fit_ibss_marginal_only(pdac_data$gene_lists[[factor]], pdac_msigdb_X)),
      tar_target(ibss_fit, fit_ibss_select_L(pdac_data$gene_lists[[factor]], pdac_msigdb_X)),
      tar_target(ibss_fit_lasso_init, fit_ibss_lasso_init(pdac_data$gene_lists[[factor]], pdac_msigdb_X))
    ),
    tar_target(pdac_lasso_fits, fit_lasso_all_gene_lists(pdac_msigdb_X, pdac_data$gene_lists)),
    tar_target(pdac_lasso_sims, simulate_from_lasso_map(pdac_msigdb_X, pdac_lasso_fits, 20))
)

#     tar_target(yusha_data, yusha_example_data()$bindata),
#     tar_target(yusha_fit_spec, make_yusha_fit_spec()),
#     tar_target(yusha_fits,
#       yusha_fit_spec %>%
#         dplyr::rowwise() %>%
#         dplyr::mutate(
#             fit = list(fit_model2(
#               yusha_data$X, yusha_data$y, fit_method, fit_fun, fit_args))
#         ),
#         pattern = map(yusha_fit_spec)
#     ),
#     tar_target(yusha_data500, yusha_example_data(500)$bindata),
#     tar_target(yusha_fits500,
#       yusha_fit_spec %>%
#         dplyr::rowwise() %>%
#         dplyr::mutate(
#             fit = list(with(yusha_data500,
#               fit_model2(X, y, fit_method, fit_fun, fit_args)))
#         ),
#         pattern = map(yusha_fit_spec)
#     )
# )
