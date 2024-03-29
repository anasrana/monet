## All class definitions

#' monetDataSuper superclass
#'
#' @slot path to gene expression
#' @slot path to ATACseq data
#'
#' @export
setClass("monet",
         contains = "VIRTUAL",
         slots = list(gene_exp_path = "character",
                      atac_seq_path = "character"))

# =============================================================================
# MONET main data classes
# =============================================================================

#' MONET data input class
#'
#' @slot gene_exp data.table. Gene expression time-course data.
#' @slot atac_seq data.table. ATAC-seq data for each time-point.
#' @slot no_tpts number of time-points.
#' @slot no_genes number of genes in the data, derived from `gene_exp`.
#' @slot no_tf_try number of trascription factors to try.
#' @slot gene_names vector of gene.names.
#'
#' @export
setClass("monetDataRaw",
         contains = "monet",
         slots = list(gene_exp = "data.table",
                      atac_seq = "data.table",
                      gene_exp_init = "data.table",
                      no_tpts = "numeric",
                      no_genes = "numeric",
                      gene_names = "vector",
                      data_test = "vector"
                      ))

#' MONET filtered data input class
#'
#'
#' @slot exp_filt filtered list of gene
#'
#' @export
setClass("monetData",
         contains = "monetDataRaw",
         slots = list(gene_filt_list = "vector"))

# ============================================================================
# MONET fit class
# ============================================================================

#' MONET optimization output class
#'
#' @slot x_est matrix
#' @slot b_est matrix
#' @slot w vector.
#'
#' @export
setClass("monetOptim",
         slots = list(x_est = "matrix",
                      b_est = "matrix",
                      w = "vector",
                      gene_names = "vector",
                      rstan_optim = "list"))

#' MONET sampling output class
#'
#' This class contains a stanfit object and some useful variables.
#' @slot x_est array of the mcmc chains.
#' @slot b_est array of the mcmc cahins.
#' @slot w
setClass("monetSample",
         slots = list(x_est = "array",
                      b_est = "array",
                      w = "matrix",
                      gene_names = "vector",
                      stanfit = "list",
                      n_iter = "numeric",
                      n_warmup = "numeric",
                      n_chains = "numeric",
                      w_up_saved = "logical"
                      ))
