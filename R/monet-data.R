#' monetInput
#'
#' Function to initialise monetInput class
#'
#' @param gene_exp gene expression. Either path to gene expression file or
#'  data.frame of gene expression. Rows represent genes and columns represent
#' @param atac_seq ATAC-Seq.
#' @param gene_col column of gene names in gene_exp either in the file or
#'  data.frame
#' @param gene_atac column of gene names in atac_seq if different from
#'  gene_col.
#'
#' @return monetInput class based on input to function.
#' @export
#'
#' @examples
#'
monetInput <- function(gene_exp = NULL,
                       atac_seq = NULL,
                       gene_col = NULL,
                       gene_atac = NULL) {

    gene_exp_dt <- readGeneExp(gene_exp, gene_col)

    return(gene_exp_dt)
}

#' readGeneExp
#'
#' Function to set format of gene expression
#'
#' @param gene_exp data or path to data.file.
#' @param gene_col gene name column.
#'
#' @importFrom data.table fread setkeyv as.data.table setnames .N
#' @importFrom dplyr nth
readGeneExp <- function(gene_exp = NULL, gene_col = NULL) {
    # TODO: if gene_col doesn't exist give warning
    if (!is.null(gene_col)) {
        gene_col <- testGeneCol(gene_col = gene_col, file_path = gene_exp)
    }

    if ("character" %in% class(gene_exp)) {
        checkFile(gene_exp)

        if (file.exists(gene_exp) & !is.null(gene_col)) {
            gene_exp <- fread(gene_exp, key = gene_col)
        } else if (file.exists(gene_exp)) {
            gene_exp <- fread(gene_exp)
            message("No gene names provided using row number...\n")
            gene_exp <- gene_exp[, "gene" := 1:.N]
            gene_col <- "gene"
        }
    } else if ("data.table" %in% class(gene_exp)) {
        if (!is.null(gene_col)) {
            gene_exp <-
                gene_exp %>%
                setkeyv(gene_col)
        } else {
            message("No gene names provided using row number...\n")
            gene_exp <- gene_exp[, "gene" := 1:.N]
            gene_col <- "gene"

            setkeyv(gene_exp, "gene")
        }
    } else if ("data.frame" %in% class(gene_exp)) {
        if (is.null(gene_col)) {
            genes_v <- row.names(gene_exp)
            gene_col <- "gene"
            message("No gene_col provided for data.frame, using row names",
                    "...")
            gene_exp <-
                gene_exp %>%
                as.data.table()

            gene_exp <- gene_exp[, "gene" := genes_v]
            setkeyv(gene_exp, "gene")
        } else {
            gene_exp <-
                gene_exp %>%
                as.data.table()

            setkeyv(gene_exp, gene_col)

        }
    }

    return(gene_exp[, setnames(.SD, gene_col, "gene")])
}
