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

#' initialize monetDataFilt
#'
setMethod("initialize",
          signature = "monetDataFilt",
          function(.Object,
                   monetData,
                   gene_filt) {
    # import from existing class
    .Object@gene_exp_path <- monetData@gene_exp_path
    .Object@atac_seq_path <- monetData@atac_seq_path
    .Object@atac_seq <- monetData@atac_seq
    .Object@gene_exp <- monetData@gene_exp
    .Object@no_genes <- monetData@no_genes
    .Object@no_tpts <- monetData@no_tpts
    .Object@gene_names <- monetData@gene_names
    .Object@gene_exp_init <- monetData@gene_exp_init
    .Object@data_test <- TRUE
    .Object@gene_overlap_list <- gene_filt
    return(.Object)
          })

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
    cat("------------------------------------\n",
        "----- MONET input data class -----\n",
        object@no_genes, " genes and   ", object@no_tpts, " time-points\n",
        "----------------------------------\n")
    if (object@gene_exp_path != "NA") {
        cat("Data imported from: ", basename(object@gene_exp_path), "\n")
        cat("Data imported from: ", basename(object@atac_seq_path), "\n")
    }
    cat("Genes names:\n", g_names[1], ", ", g_names[2], "...",
        g_names[length(g_names) - 1], ", ", g_names[length(g_names)], "\n")
    if (object@data_test != TRUE) {
        warning("No of rows different between gene exp and atac seq\n",
                "MONET won't run until the data has been filtered",
                call. = FALSE)
    }
  })

#' Print method for monetDataFilt
#'
#' @param monetDataFilt
#'
#' @export
#'
setMethod("show", "monetDataFilt",
  function(object) {
    g_names <- object@gene_names
    cat("------------------------------------\n",
        "------------ FILTERED -------------\n",
        "----- MONET input data class -----\n",
        length(object@gene_overlap_list), "of", object@no_genes,
        "genes and   ", object@no_tpts, " time-points\n",
        "----------------------------------\n")
    if (object@gene_exp_path != "NA") {
        cat("Data imported from: ", basename(object@gene_exp_path), "\n")
        cat("Data imported from: ", basename(object@atac_seq_path), "\n")
    }
    cat("Genes names:\n", g_names[1], ", ", g_names[2], "...",
        g_names[length(g_names) - 1], ", ", g_names[length(g_names)], "\n")
    if (object@data_test != TRUE) {
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
