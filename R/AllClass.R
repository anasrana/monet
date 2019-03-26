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
#' @slot
#'
#' @export
setClass("monetInput",
         contains = "monet",
         slots = list(gene_exp = "data.table",
                      atac_seq = "data.table",
                      no_tpts = "numeric",
                      no_genes = "numeric",
                      no_tf_try = "numeric",
                      gene_names = "vector"
                      ))
