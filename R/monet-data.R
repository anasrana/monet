#' monetData initialize
#'
#' Function to initialise monetData class
#'
#' @param gene_exp gene expression. Either path to gene expression file or
#'  data.frame of gene expression. Rows represent genes and columns represent
#' @param atac_seq ATAC-Seq.
#' @param gene_col column of gene names in gene_exp either in the file or
#'  data.frame
#' @param gene_atac column of gene names in atac_seq if different from
#'  gene_col.
#'
#' @return monetData class based on input to function.
#'
#' @importFrom data.table key .N
#'
#' @export
#'
#' @examples
#'
monetData_init <- function(gene_exp = NULL,
                           atac_seq = NULL,
                           gene_col = NULL,
                           gene_atac = NULL,
                           gene_exp0 = NULL) {


    gene_exp_dt <- prepGeneExp(gene_exp, gene_col)

    # If there is no gene_atac variable overwrite with gene_col
    if (is.null(gene_atac)) {
        gene_atac <- gene_col
    }

    if (is.null(gene_exp0)) {
        gene_exp0 <- gene_exp_dt[, -"gene"][, 1]
        gene_exp0 <- cbind(gene_exp_dt[, "gene"], gene_exp0)
        col_tf <- colnames(gene_exp0) != colnames(gene_exp_dt)
        col_nm <- c("gene", colnames(gene_exp_dt)[col_tf])
        gene_exp_dt <- gene_exp_dt[, ..col_nm]
    } else {
        gene_exp0 <- prepExpT0(gene_exp0, gene_col)
    }


    atac_seq_dt <- prepAtac(atac_seq, gene_atac)
    gene_names <- atac_seq_dt[, get(key(atac_seq_dt))]

    if (!("character" %in% class(gene_exp))) {
        gene_exp <- "NA"
    }
    if (!("character" %in% class(gene_atac))) {
        gene_atac <- "NA"
    }

    if (ncol(atac_seq_dt) != ncol(gene_exp_dt)) {

        name_int <- intersect(colnames(atac_seq_dt[, -c("gene")]),
            colnames(gene_exp_dt[, -c("gene")]))

        if (length(name_int) == 0) {
            stop("Mismatch in time points between RNAseq and ATACseq.\n",
                 "Ensure you provide gene expression at t = 0 separately.")
        } else {
            gene_exp_dt <- gene_exp_dt[, ..name_int]
            atac_seq_dt <- atac_seq_dt[, ..name_int]
        }
    }

    monet_input <- new("monetData",
                       gene_dt = gene_exp_dt,
                       atac_dt = atac_seq_dt,
                       gene_nm = gene_names,
                       gene_t0_dt = gene_exp0,
                       gene_path = gene_exp,
                       atac_path = atac_seq)

    return(monet_input)
}

#' Prepare gene expression at
prepExpT0 <- function(gene_exp0 = NULL, gene_col = NULL) {
    if ("character" %in% class(gene_exp0)) {
        checkFile(gene_exp0)
        if (!is.null(gene_col)) {
            gene_exp0 <- fread(gene_exp0, key = gene_col)
        } else if (file.exists(gene_exp0)) {
            gene_exp0 <- fread(gene_exp0)
            message("No gene names provided using row number ...\n")
            gene_exp0 <- gene_exp0[, "gene" := 1:.N]
            gene_col <- "gene"
            set_key_dt(gene_exp0, gene_col)
        }
    } else if ("data.table" %in% class(gene_exp0)) {
        if (!is.null(gene_col)) {
            gene_exp0 <-
                gene_exp0 %>%
                set_key_dt(gene_col)
        } else {
            message("No gene names provided using row number...\n")
            gene_exp0 <- gene_exp0[, "gene" := 1:.N]
            gene_col <- "gene"

            set_key_dt(gene_exp0, gene_col)
        }
    } else if ("data.frame" %in% class(gene_exp0)) {
        if (is.null(gene_col)) {
            genes_v <- row.names(gene_exp0)
            gene_col <- "gene"
            message("No gene_col provided for data.frame, using row names",
                    "...")
            gene_exp0 <-
                gene_exp0 %>%
                as.data.table()

            gene_exp0 <- gene_exp0[, "gene" := genes_v]
            set_key_dt(gene_exp0, "gene")
        } else {
            gene_exp0 <-
                gene_exp0 %>%
                as.data.table()

            set_key_dt(gene_exp0, gene_col)

        }
    }

    return(gene_exp0[, setnames(.SD, gene_col, "gene")])
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

    return(atac_seq[, setnames(.SD, gene_col, "gene")])

}
