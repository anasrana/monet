#' Filter genes of monetData class
#'
#' MONET will use this function to ensure RNAseq data and ATAC seq data has
#' data for the same genes. The simplest form subsets by genes found in both
#' tables.
#'
#' @param monet_data monetData class.
#' @param gene_filt vector. A vector containing genes to be filtered out.
#' @param noise_flt numeric. Threshold for signal to noise cutoff.
#' @param no_genes integer. No of genes to filter out based on high signal to
#'  noise in gene expression.
#' @param mu_th numeric. Default `=0.1`. Remove genes with average gene
#'  expression across time `< mu_th`
#' @param sd_th numeric. Default `=0.1`. Remove genes with standard deviation
#'  of gene expression across time.
#'
#' @export
#'
#' @importFrom data.table melt setkeyv
monet_data_filter <- function(monet_data,
                        gene_filt = NULL,
                        noise_flt = NULL,
                        no_genes = NULL,
                        mu_th = 0.1,
                        sd_th = 0.1) {
    if (!("monetDataRaw" %in% class(monet_data))) {
        stop("monet_data has to be a monetData class")
    }

    # perform initial filteration based on input information

    gen_l <- filter_stdd(monet_data, no_genes, noise_flt)

    gene_dt <- gen_l[["g_dt"]]
    gene_i_dt <- gen_l[["g_t0"]]

    if (!is.null(gene_filt)) {
        gene_dt <- gene_dt[gene_filt]
    }

    monet_data <- setGeneExp(monet_data, gene_dt)
    monet_data <- setGeneInit(monet_data, gene_i_dt)

    atac_dt <- getAtacSeq(monet_data)

    gene_list <- overlapGenes(gene_dt, atac_dt)

    monet_data <- new("monetData",
                      monet_data,
                      gene_list)

    monet_data <- setDataTest(monet_data, TRUE)

    return(monet_data)
}


overlapGenes <- function(gene_dt, atac_dt) {
    tmp <- merge(gene_dt, atac_dt, by = "gene")
    na.omit(tmp)[, get("gene")]
}

#' Filter and standardise gene expression
#'
#' @param monet_data monetData class.
#' @inheritParams monet_data_filter
#'
#' @importFrom data.table dcast
#'
filter_stdd <- function(monet_data,
                        no_genes = NULL,
                        noise_flt = NULL,
                        mu_th = 0.1,
                        sd_th = 0.1) {

    gene_dt <- getGeneExp(monet_data)
    gene_i_dt <- getGeneInit(monet_data)

    init_rn <- colnames(gene_i_dt)
    gene_rn <- colnames(gene_dt)

    gene_dt <- gene_dt[gene_i_dt, nomatch = NA]
    no_tpt <- ncol(gene_dt[, -"gene"])
    gene_melt <- geneDtMelt(gene_dt, mu_th, sd_th)

    if (!is.null(noise_flt)) {
        gene_melt <-
        gene_melt[sn > noise_flt]
    } else if (!is.null(no_genes)) {
        gene_melt <-
        gene_melt[order(-sd_v)
                ][1:(no_tpt * no_genes), ]
    }

    gene_dt <-
    geneStanderdise(gene_melt) %>%
        dcast(gene ~ variable)

    gene_i_dt <- gene_dt[, ..init_rn]
    gene_dt <- gene_dt[, ..gene_rn]

    return(list(g_t0 = gene_i_dt, g_dt = gene_dt))
}
