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
setClass("monetData",
         contains = "monet",
         slots = list(gene_exp = "data.table",
                      atac_seq = "data.table",
                      no_tpts = "numeric",
                      no_genes = "numeric",
                      no_tf_try = "numeric",
                      gene_names = "vector"
                      ))

#' MONET filtered data input class
#'
#'
#' @slot exp_filt filtered list of gene
#'
#' @export
setClass("monetDataFilt",
         contains = "monetData",
         slots = list(exp_filt = "vector",
                      atac_filt = "vector"))