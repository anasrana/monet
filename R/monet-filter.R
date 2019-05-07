#' Filter genes of monetData class
#'
#' MONET will use this function to ensure RNAseq data and ATAC seq data has
#' data for the same genes. The simplest form subsets by genes found in both
#' tables.
#'
#' @param monet_data monetData class.
#' @param gene_filt vector. A vector containing genes to be filtered out.
#' @param noise_flt numeric. Threshold for signal to noise cutoff.
#' @param no_genes integer. No of genes to filter out based on signal to noise
#'  in gene expression.
#'
#' @export
#'
#' @importFrom data.table melt setkeyv
monet_data_filter <- function(monet_data,
                        gene_filt = NULL,
                        noise_flt = NULL,
                        no_genes = NULL) {
    if (!("monetDataRaw" %in% class(monet_data))) {
        stop("monet_data has to be a monetData class")
    }

    gene_dt <- getGeneExp(monet_data)
    atac_dt <- getAtacSeq(monet_data)

    if (is.null(gene_filt) & is.null(noise_flt) & is.null(no_genes)) {
        message("Finding overlapping genes between RNAseq and ATACseq")
        gene_list <- overlapGenes(gene_dt, atac_dt)
    } else if (!is.null(gene_filt)) {
        message("Checking gene list against RNAseq and ATACseq")
        gene_list <- overlapGenes(gene_dt[gene_filt], atac_dt[gene_filt])
    } else if (!is.null(noise_flt)) {
        gene_mlt <-
        melt(gene_dt,
             id.vars = "gene",
             variable.name = "t_point",
             value.name = "g_exp") %>%
        setkeyv("gene")

        gn_v <-
        gene_mlt[, .(mu = mean(g_exp), sig = sd(g_exp)), by = "gene"
               ][, .(sn = mu / sig), by = "gene"
               ][sn >= noise_flt][, get("gene")]

        gene_list <- overlapGenes(gene_dt[gn_v], atac_dt[gn_v])
    } else if (!is.null(no_genes)) {
        gene_mlt <-
        melt(gene_dt,
             id.vars = "gene",
             variable.name = "t_point",
             value.name = "g_exp") %>%
        setkeyv("gene")

        gn_v <-
        gene_mlt[, .(mu = mean(g_exp), sig = sd(g_exp)), by = "gene"
               ][, .(sn = mu / sig), by = "gene"
               ][order(-sn)][1:no_genes
               ][, get("gene")]

        gene_list <- overlapGenes(gene_dt[gn_v], atac_dt[gn_v])
    }

    monet_filt <- new("monetData",
                      monet_data,
                      gene_list)

    monet_filt <- setDataTest(monet_filt, TRUE)

    return(monet_filt)
}


overlapGenes <- function(gene_dt, atac_dt) {
    tmp <- merge(gene_dt, atac_dt, by = "gene")
    na.omit(tmp)[, get("gene")]
}
