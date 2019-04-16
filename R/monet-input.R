#'  Monet input prep
#'
#' Prepare input monet class for Rstan model.
#'
#' @param monet_class monet data class. Possibilities include `monetData` or
#'  `monetDataFilt` classes.
#' @param par_list list. This list obhject should contain any predefined
#'  parameters, this needs to be a names list. The following named parameters
#'  can be passed: `no_tf`, `de_b`, `sigma`, `sigma_x`, `alpha`
#'
#' @return list for monet stan model.
#'
#' @export
monet_input_prep <- function(monet_class,
                             trnsfrm_ge = c("asinh", "log", "none"),
                             par_list = list()) {
    # only possible list name
    list_nm <- c("no_tf", "de_b", "sigma", "sigma_x", "alpha")
    # predefined list
    pre_def_par <- list("no_tf" = 30,
                        "de_b" = 2,
                        "sigma" = 1,
                        "sigma_x" = 1,
                        "alpha" = 1)

    trnsfrm_ge <- match.arg(trnsfrm_ge)
    if (!all(names(par_list) %in% list_nm)) {
        not_found <- setdiff(names(par_list), list_nm)
        stop("Not all parameters provided are from the list:",
             list_nm)
    } else if (length(par_list) == 0) {
        par_list <- pre_def_par
    } else {
        new_nm <- setdiff(list_nm, names(par_list))
        par_list <- c(par_list, pre_def_par[new_nm])
    }

    if (class(monet_class) == "monetDataFilt") {
        gn_lst <- monet_class@gene_filt_list
        gene_dt <- monet_class@gene_exp[gn_lst][, -"gene"]
        atac_dt <- monet_class@atac_seq[gn_lst][, -"gene"]
        gene_init <- monet_class@gene_exp_init[gn_lst][[2]]
    } else if (class(monet_class) == "monetData") {
        if (monet_class@data_test != TRUE) {
            stop("monetData object can't be analysed, it needs to be filtered")
        } else {
            gn_lst <- monet_class@gene_names
            gene_dt <- monet_class@gene_exp[, -"gene"]
            atac_dt <- monet_class@atac_seq[, -"gene"]
            gene_init <- monet_class@gene_exp_init[[2]]
        }
    }

    # add monetData object
    par_list[["no_gns"]] <- length(gn_lst)
    par_list[["no_tpt"]] <- ncol(atac_dt)
    par_list[["a"]] <- t(as.matrix(atac_dt))

    if (trnsfrm_ge == "none") {
        par_list[["y"]] <- t(as.matrix(gene_dt))
        par_list[["y0"]] <- gene_init
    } else if (trnsfrm_ge == "log") {
        par_list[["y"]] <- log10(t(as.matrix(gene_dt)) + 1)
        par_list[["y0"]] <- log10(gene_init + 1)
    } else if (trnsfrm_ge == "asinh") {
        par_list[["y"]] <- asinh(t(as.matrix(gene_dt)))
        par_list[["y0"]] <- asinh(gene_init)
    }
    return(par_list)

}
