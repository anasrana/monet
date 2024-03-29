% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stan_monet.fit.R
\name{stan_monet.fit}
\alias{stan_monet.fit}
\title{Stan_monet}
\usage{
stan_monet.fit(standata, ..., algorithm = c("sampling", "optimizing",
  "fullrank", "mcmc"))
}
\arguments{
\item{algorithm}{character. Choose option for stan model.}

\item{model_par}{list conttaining model parameters. Use \link{monet_input_prep}
to create from monetData objects.}

\item{trnsfrm_ge}{character. Choose if gene expression data needs to be
transformed.}
}
\description{
Model fit for monet.
}
