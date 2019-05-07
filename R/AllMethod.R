# =============================================================================
# Initialising classes
# =============================================================================

#' initilize monetDataRaw
#'
#' @param monetDataRaw
#'
#' @return  monetDataRaw object
#'
#' @importFrom data.table melt dcast setkeyv
#'
setMethod("initialize",
          signature = "monetDataRaw",
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
    # calcualte based on supplied varaible

    tmp_noG <- nrow(gene_dt)

    g0_name <- colnames(gene_t0_dt)[2]
    g_cname <- colnames(gene_dt)
    sd_test <-
        gene_dt[gene_t0_dt, nomatch = NA] %>%
        melt(id.vars = "gene") %>%
        setkeyv("gene")

    gene_dt <-
    sd_test[, sd_v := sd(value), by = "gene"
          ][sd_v > 0.1
          ][, -"sd_v"] %>%
          dcast(gene ~ variable)

    if(nrow(gene_dt) != tmp_noG) {
        message("\n", tmp_noG - nrow(gene_dt),
                " genes with small sd removed.\n")
        gene_nm <- gene_dt[, get("gene")]
    }

    init_nm <- c("gene", g0_name)
    gene_t0_dt <- gene_dt[, ..init_nm]
    gene_dt <- gene_dt[, ..g_cname]

    .Object@data_test <-
        all.equal(gene_dt[, "gene"], atac_dt[, "gene"],
                  ignore.row.order = TRUE)

    tmp_dim <- dim(gene_dt)
    .Object@gene_exp <- gene_dt
    .Object@gene_exp_init <- gene_t0_dt
    .Object@no_genes <- tmp_dim[1]
    .Object@no_tpts <- tmp_dim[2] - 1
    .Object@gene_names <- gene_nm

    return(.Object)
    }
)

#' initialize monetData
#'
setMethod("initialize",
          signature = "monetData",
          function(.Object,
                   monetDataRaw,
                   gene_filt) {
      #   [, value := (value - mean(value)) / sd_v, by = "gene"] %>%
      # dcast(gene ~ variable) %>%
      # setkeyv("gene")
    # import from existing class
    .Object@gene_exp_path <- monetDataRaw@gene_exp_path
    .Object@atac_seq_path <- monetDataRaw@atac_seq_path
    .Object@atac_seq <- monetDataRaw@atac_seq
    .Object@gene_exp <- monetDataRaw@gene_exp
    .Object@no_genes <- monetDataRaw@no_genes
    .Object@no_tpts <- monetDataRaw@no_tpts
    .Object@gene_names <- monetDataRaw@gene_names
    .Object@gene_exp_init <- monetDataRaw@gene_exp_init
    .Object@data_test <- TRUE
    .Object@gene_filt_list <- gene_filt
    return(.Object)
          })

#' initialize monetOptim class
#'
#' @param monetOptim
#'
#' @return monetOptim object
#'
#' @importFrom magrittr set_colnames set_rownames set_names
#' @importFrom stringr str_c
#'
setMethod("initialize",
          signature = "monetOptim",
          function(.Object,
                   optim_out,
                   monet_dat) {
    # extract parameters
    x_mat <- optim_out$par[grepl("x", names(optim_out$par))]
    b_mat <- optim_out$par[grepl("b", names(optim_out$par))]
    w_vec <- optim_out$par[grepl("w", names(optim_out$par))]

    if (class(monet_dat) == "monetData") {
        g_names <- monet_dat@gene_names
    } else if (class(monet_dat) == "monetDataFilt") {
        g_names <- monet_dat@gene_filt_list
    }
    no_tf <- length(w_vec)
    no_tpts <- monet_dat@no_tpts
    no_gns <- length(g_names)

    x_mat <- matrix(x_mat, nrow = no_tpts, ncol = no_tf, byrow = F) %>%
        set_colnames(stringr::str_c("TF.", 1:no_tf)) %>%
        set_rownames(stringr::str_c("t = ", 1:no_tpts))

    b_mat <- matrix(b_mat, nrow = no_gns, ncol = no_tf, byrow = F) %>%
        set_colnames(stringr::str_c("TF.", 1:no_tf)) %>%
        set_rownames(g_names)

    w_vec <- w_vec %>%
        set_names(stringr::str_c("TF.", 1:no_tf))

    .Object@x_est <- x_mat
    .Object@b_est <- b_mat
    .Object@w <- w_vec
    .Object@gene_names <- g_names
    .Object@rstan_optim <- optim_out

    return(.Object)
          })

# =============================================================================
# Printing methods
# =============================================================================

#' Print method for monetDataRaw
#'
#' @param monetDataRaw
#'
#' @export
#'
setMethod("show", "monetDataRaw",
  function(object) {
    g_names <- object@gene_names
    cat("------------------------------------\n",
        "---- MONET Raw input data class ----\n",
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
        "------------ FILTERED -------------\n",
        "----- MONET input data class -----\n",
        length(object@gene_filt_list), "of", object@no_genes,
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
          signature = "monetDataRaw",
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
          signature = "monetDataRaw",
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
          signature = "monetDataRaw",
          definition = function(.Object, data_test) {
            .Object@data_test <- data_test
            return(.Object)
          })
