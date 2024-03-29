% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_fun.R
\name{geneDtMelt}
\alias{geneDtMelt}
\title{Melt gene_dt}
\usage{
geneDtMelt(gene_dt, mu_th = 0.1, sd_th = 0.1)
}
\arguments{
\item{gene_dt}{data.table. FPKM transformed gene expression data.table.}

\item{mu_th}{numeric. Threshold for mean values.}

\item{sd_th}{numeric. Threshold for sd values.}
}
\description{
Reshape the gene_dt object and prepare for filtering
}
