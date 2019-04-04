# =============================================================================
# Initialising classes
# =============================================================================

#' initilize monetData
#'
#' @param monetData
#'
#' @return  monetData object
#'
setMethod("initialize",
          signature = "monetData",
          function(.Object,
                   gene_dt,
                   atac_dt,
                   gene_nm,
                   gene_t0_dt,
                   gene_path = "NA",
                   atac_path = "NA") {
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
    .Object@gene_exp_init <- gene_t0_dt

    .Object@data_test <-
        all.equal(gene_dt[, "gene"], atac_dt[, "gene"],
                  ignore.row.order = TRUE)

    return(.Object)
    }
)

# =============================================================================
# Printing methods
# =============================================================================

#' Print method for monetData
#'
#' @param monetData
#'
#' @export
#'
setMethod("show", "monetData",
  function(object) {
    g_names <- object@gene_names
    cat("----------------------------------\n",
        "----- MONET input data class -----\n",
        object@no_genes, " genes and   ", object@no_tpts, " time-points\n",
        "----------------------------------\n")
    if (object@gene_exp_path != "NA") {
        cat("Data imported from: ", basename(object@gene_exp_path), "\n")
        cat("Data imported from: ", basename(object@atac_seq_path), "\n")
    }
    cat("Genes names:\n", g_names[1], ", ", g_names[2], "...",
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

setGeneric(name = "getAtacSeq",
           def = function(.Object) {
            standardGeneric("getAtacSeq")
           })

setGeneric(name = "setDataTest",
           def = function(.Object, data_test) {
            standardGeneric("setDataTest")
           })

# =============================================================================
# GENERIC implemented
# =============================================================================

#' Extract gene expression slot
#'
#' @param monetData monetData class object.
#'
#' @export
#'
setMethod(f = "getGeneExp",
          signature = "monetData",
          definition = function(.Object) {
            gene_exp <- .Object@gene_exp
            return(gene_exp)
          })
#' Extract ATACseq slot
#'
#' @param monetData monetData class object.
#'
#' @export
setMethod(f = "getAtacSeq",
          signature = "monetData",
          definition = function(.Object) {
            return(.Object@atac_seq)
          })

#' Update data_test slot
#'
#' @param monetData monetData class object.
#' @param data_test new data_test entry.
#'
#' @export
setMethod(f = "setDataTest",
          signature = "monetData",
          definition = function(.Object, data_test) {
            .Object@data_test <- data_test
            return(.Object)
          })
