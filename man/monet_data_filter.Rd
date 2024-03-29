% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/monet-filter.R
\name{monet_data_filter}
\alias{monet_data_filter}
\title{Filter genes of monetData class}
\usage{
monet_data_filter(monet_data, gene_filt = NULL, noise_flt = NULL,
  no_genes = NULL, mu_th = 0.1, sd_th = 0.1)
}
\arguments{
\item{monet_data}{monetData class.}

\item{gene_filt}{vector. A vector containing genes to be filtered out.}

\item{noise_flt}{numeric. Threshold for signal to noise cutoff.}

\item{no_genes}{integer. No of genes to filter out based on high signal to
noise in gene expression.}

\item{mu_th}{numeric. Default \code{=0.1}. Remove genes with average gene
expression across time \code{< mu_th}}

\item{sd_th}{numeric. Default \code{=0.1}. Remove genes with standard deviation
of gene expression across time.}
}
\description{
MONET will use this function to ensure RNAseq data and ATAC seq data has
data for the same genes. The simplest form subsets by genes found in both
tables.
}
