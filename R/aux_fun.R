#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


#' data.table operator
#'
#' @name :=
#' @rdname assignment
#' @keywords internal
#' @export
#' @importFrom data.table ":="
`:=` <- function(...) NULL


#' Wrapper around file.exists
#'
#' Need NA return for monet if file does not exist.
#'
#' @param file_path File path to be tested.
#'
#' @return file_path or "NA"
filePath <- function(file_path) {
    if ("character" %in% class(file_path)) {
        if (file.exists(file_path)) {
            file_path
        }
    } else {
        "NA"
    }
}

#' Check if file exists
#'
#' Check if the specified file exists and give better error message.
#'
#' @param file_path character. String specifying location of file to be read.
checkFile <- function(file_path) {
    if (!file.exists(file_path)) {
        stop("(MONET data) - file specified:\n", file_path,
             "\ndoes not exist, check and provide the full path.",
             call. = FALSE)
    } else {
        message("Reading file:\t", file_path)
    }
}

#' testing gene_col var
#'
#' Testing if the gene_col variable provided is contained within the data.
#' Also converts numerical column indicator into characters.
#'
#' @param gene_col character or integer. The column which contains gene names.
#' @param file_path character. Location of file to be read.
#'
#' @importFrom data.table fread
#' @importFrom dplyr nth
#'
testGeneCol <- function(gene_col, file_path) {
    gene_num <- suppressWarnings(as.integer(gene_col))

    if (!is.na(gene_num)) {
        gene_col  <- gene_num
    }

    if ("character" %in% class(file_path) & !is.na(gene_num)) {
        gene_col <- fread(file_path, nrows = 1) %>%
            colnames() %>%
            nth(gene_col)
    } else if (!is.na(gene_num) & is.data.frame(file_path)) {
        gene_col <- colnames(file_path)[gene_col]
    }

    if (is.na(gene_num)) {
        if (!(gene_col %in% colnames(fread(file_path, nrows = 1)))) {
            stop("in monetInput-\nNo column '",
                 gene_col,
                 "' found in input data:- ",
                 basename(file_path),
                 call. = FALSE)
        }
    }

    return(gene_col)
}

#' geneColCheck
#'
#' if no gene col passed check if character column exists otherwise use index
#'
#' @param data_dt data.table. Data in `data.table` format
#'
geneColCheck <- function(data_dt) {
    class_v <- sapply(data_dt, class)
    col_nm_v <- colnames(data_dt)
    if ("character" %in% class_v) {
        gene_col <- col_nm_v["character" == class_v]
    } else {
         gene_col <- "index"
    }

    return(gene_col)
}

#' Set key in data.table
#'
#' helper function for [`setkeyv()`][data.table::setkeyv()] function for more
#'  useful error messages in the context of monet.
#'
#' @param data_dt data.table. A `data.table` object.
#' @param col_name the column containing the key variable.
#'
#' @importFrom data.table setkeyv
#'
set_key_dt <- function(data_dt, col_name) {
    tryCatch(setkeyv(data_dt, col_name),
             error = function(e) {
                message("The column name for genes '", col_name,
                        "' not found in the data.")
             })
}

#' Transforming data function
#'
#' helper function
#'
#' @param data_dt data.table. A `data.table` object of gene expression.
#' @param trnsfrm_ge character.
#'
#' @importFrom data.table key
trnsfrmGe <- function(data_dt, trnsfrm_ge) {
    data_key <- key(data_dt)
    if (trnsfrm_ge == "none") {
        return(data_dt)
    } else if (trnsfrm_ge == "asinh") {
        cbind(data_dt[, ..data_key], asinh(data_dt[, -..data_key])) %>%
            return()
    } else if (trnsfrm_ge == "log") {
        cbind(data_dt[, ..data_key], log2(data_dt[, -..data_key]) + 1) %>%
            return()
    } else {
        stop("Transformation (", trnsfrm_ge, ") not implemented in monet.\n",
             "Contact the developer (a.a.rana@bham.ac.uk) if you think it",
             " should be implemented")
    }
}

#' Melt gene_dt
#'
#' Reshape the gene_dt object and prepare for filtering
#'
#' @param gene_dt data.table. FPKM transformed gene expression data.table.
#' @param mu_th numeric. Threshold for mean values.
#' @param sd_th numeric. Threshold for sd values.
#'
#' @importFrom data.table melt setkeyv
#'
geneDtMelt <- function(gene_dt, mu_th = 0.1, sd_th = 0.1) {
    no_gene <- nrow(gene_dt)
    gene_melt <-
        gene_dt %>%
        melt(id.vars = "gene") %>%
        setkeyv("gene")

    gene_melt <-
    gene_melt[, c("sd_v", "mu_v") :=
                .(sd(value, na.rm = TRUE), mean(value, na.rm = TRUE)),
                by = "gene"
          ][, sn := mu_v / sd_v
          ][mu_v > mu_th & sd_v > sd_th]

    no_gn <- length(unique(gene_melt$gene))

    if (no_gn < no_gene) {
        message(no_gene - no_gn, " genes removed due to lack of variation or ",
         "very low expression")
    }

    return(gene_melt)
}


#' standardise gene expression
#'
#' standardise gene expression in long format.
#'
#' @param gene_long data.table.
#'
geneStanderdise <- function(gene_long) {
gene_long[, value := (value - mean(value, na.rm = TRUE)) /
                        sd(value, na.rm = TRUE), by = "gene"]

}
