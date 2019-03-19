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


checkFile <- function(file_path) {
    if (!file.exists(file_path)) {
        stop("(MONET data) - file specified:\n", file_path,
             "\ndoes not exist, check and provide the full path.",
             call. = FALSE)
    } else {
        message("Reading file:\t", file_path)
    }
}

testGeneCol <- function(gene_col, file_path) {
    gene_num <- suppressWarnings(as.numeric(gene_col))

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

    return(gene_col)
}
