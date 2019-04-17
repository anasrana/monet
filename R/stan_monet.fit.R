#' Stan_monet
#'
#' Model fit for monet.
#'
#' @param model_par list conttaining model parameters. Use [monet_input_prep]
#'  to create from monetData objects.
#' @param trnsfrm_ge character. Choose if gene expression data needs to be
#'  transformed.
#' @param algorithm character. Choose option for stan model.
#'
#' @export
#'
#' @importFrom rstan optimizing
stan_monet.fit <- function(standata,
                           ...,
                           algorithm =
                            c("sampling", "optimizing", "fullrank")) {

    algorithm <- match.arg(algorithm)
    stanfit <- stanmodels$monet

    if (algorithm == "optimizing") {
        optimizing_args <- list(...)
        if (is.null(optimizing_args$draws)) {
            optimizing_args$draws <- 100000L
        }
        optimizing_args$object <- stanfit
        optimizing_args$data <- standata
        optimizing_args$constrained <- TRUE
        out <- suppressWarnings(do.call(optimizing, args = optimizing_args))

        return(out)
    }
}
