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

    gene_exp_dt <- prepGeneExp(gene_exp, gene_col)
    if (is.null(gene_atac)) {
        gene_atac <- gene_col
    }

    atac_seq_dt <- prepAtac()

    return(gene_exp_dt)
}

#' prepGeneExp
#'
#' Prepare gene expression data for usage in monet
#'
#' @param gene_exp data or path to data.file.
#' @param gene_col gene name column.
#'
#' @importFrom data.table fread as.data.table setnames .N
#' @importFrom dplyr nth
prepGeneExp <- function(gene_exp = NULL, gene_col = NULL) {
    # TODO: if gene_col doesn't exist give warning
    if (!is.null(gene_col)) {
        gene_col <- testGeneCol(gene_col = gene_col, file_path = gene_exp)
    }

    if ("character" %in% class(gene_exp)) {
        checkFile(gene_exp)

        if (!is.null(gene_col)) {
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
                set_key_dt(gene_col)
        } else {
            message("No gene names provided using row number...\n")
            gene_exp <- gene_exp[, "gene" := 1:.N]
            gene_col <- "gene"

            set_key_dt(gene_exp, "gene")
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
            set_key_dt(gene_exp, "gene")
        } else {
            gene_exp <-
                gene_exp %>%
                as.data.table()

            set_key_dt(gene_exp, gene_col)

        }
    }

    return(gene_exp[, setnames(.SD, gene_col, "gene")])
}

#' prepAtac
#'
#' prepare ATAC seq data for usage in monet
#'
#' @param atac_seq
#' @param gene_col
#'
#' @importFrom data.table fread .N
#'
prepAtac <- function(atac_seq, gene_col = NULL) {
    if (!is.null(gene_col)) {
        gene_col <- testGeneCol(gene_col = gene_col, file_path = atac_seq)
    }

    if ("character" %in% class(atac_seq)) {
    # atac_seq is a file
        checkFile(atac_seq)

        if (!is.null(gene_col)) {
            gene_col <- testGeneCol(gene_col = gene_col, file_path = atac_seq)
            atac_seq <- fread(atac_seq, key = gene_col)

        } else {
            atac_seq <- fread(atac_seq)

            gene_col <- geneColCheck(atac_seq)
            if (gene_col == "index") {
                gene_col <- "gene"
                atac_seq <- atac_seq[, "gene" := 1:.N]
            }
            set_key_dt(atac_seq, gene_col)

        }

    } else if ("data.table" %in% class(atac_seq)) {
        if (is.null(gene_col)) {
            gene_col <- geneColCheck(atac_seq)
            if (gene_col == "index") {
                gene_col <- "gene"
                atac_seq <- atac_seq[, "gene" := 1:.N]
            }
        }

        atac_seq <- set_key_dt(atac_seq, gene_col)
    } else if ("data.frame" %in% class(atac_seq)) {
        if (is.null(gene_col)) {
            genes_v <- row.names(atac_seq)
            gene_col <- "gene"
            message("No gene_col provided for data.frame, using row names",
                    "...")
            atac_seq <-
                atac_seq %>%
                as.data.table()

            atac_seq <- atac_seq[, "gene" := genes_v]
            set_key_dt(atac_seq, "gene")
        } else {
            atac_seq <- as.data.table(atac_seq)

            set_key_dt(atac_seq, gene_col)

        }
    }

}
