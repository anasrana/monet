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
#' @importFrom rstan optimizing sampling
stan_monet.fit <- function(standata,
                           ...,
                           algorithm = c("sampling", "optimizing"),
                           cores = 1L) {
    algorithm <- match.arg(algorithm)
    stanfit <- stanmodels$monet


    cores = getOption("mc.cores", 1L)
    if (algorithm == "optimizing") {
        optim_arg <- list(...)
        if (is.null(optim_arg$draws)) {
            optim_arg$draws <- 100000L
        }
        optim_arg$object <- stanfit
        optim_arg$data <- standata
        optim_arg$constrained <- TRUE
        out <- suppressWarnings(do.call(optimizing, args = optim_arg))

        return(out)
    } else if (algorithm == "sampling") {
        sampling_pars <-
            prepSampling(stanfit_obj = stanfit,
                         data = standata,
                         user_dots = list(...),
                         show_messages = FALSE,
                         cores = cores)

        monet_fit <- do.call(sampling, sampling_pars)

        # TODO: include checks for results
        return(monet_fit)
    }
}
