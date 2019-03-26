#' initilize monetInput
#'
#' @param monetInput
#'
#' @return  monetInput object
#'
setMethod("initialize",
          signature = "monetInput",
          function(.Object,
                   gene_dt,
                   atac_dt,
                   gene_nm,
                   gene_path = "NA",
                   atac_path = "NA",
                   no_tf = 30) {
    # assigning values
    .Object@gene_exp_path <- gene_path
    .Object@atac_seq_path <- atac_path
    .Object@atac_seq <- atac_dt
    .Object@gene_exp <- gene_dt
    # calcualte based on supplied varaible
    tmp_dim <- dim(gene_dt)
    .Object@no_genes <- tmp_dim[1]
    .Object@no_tpts <- tmp_dim[2] - 1
    .Object@gene_names <- gene_nm

    return(.Object)
    }
)

#' Print method for monetInput
#'
#' @param monetInput
#'
#' @export
#'
setMethod("show", "monetInput",
  function(object) {
    g_names <- object@gene_names
    cat("----------------------------------\n",
        "----- MONET input data class -----\n",
        object@no_genes, " genes and   ", object@no_tpts, " time-points\n",
        "----------------------------------\n")
    cat("Genes names: ", g_names[1], ", ", g_names[2], "...",
        g_names[length(g_names) - 1], ", ", g_names[length(g_names)], "\n")
    if (nrow(object@gene_exp) != nrow(object@atac_seq)) {
        warning("No of rows different between gene exp and atac seq\n",
                "MONET won't run until the data has been filtered",
                call. = FALSE)
    }
  })

# =============================================================================
# GENERICS
# =============================================================================

setGeneric(name = "getGeneExp",
           def = function(.Object) {
            standardGeneric("getGeneExp")
           })

setGeneric(name = "getAtac",
           def = function(.Object) {
            standardGeneric("getAtac")
           })
