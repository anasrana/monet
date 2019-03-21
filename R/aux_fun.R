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
    if (!is.integer(gene_col) | !is.character(gene_col)) {
        stop("The gene column needs to be an integer location or the name.",
             call. = FALSE)
    }

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
