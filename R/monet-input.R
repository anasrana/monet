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
                             par_list = list()) {
    # only possible list name
    list_nm <- c("no_tf", "de_b", "sigma", "sigma_x", "alpha")
    # predefined list
    pre_def_par <- list("no_tf" = 30,
                        "de_b" = 2,
                        "sigma" = 1,
                        "sigma_x" = 1,
                        "alpha" = 1)

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


    if (class(monet_class) == "monetDataRaw") {
        warning("monetDataRaw object provided, default filtering used.",
                "\nTo use your own call 'monet_data_filter()' first.")

        monet_class <-
            monet_data_filter(monet_class, no_genes = 1000)
    }

    gn_lst <- getGeneList(monet_class)
    gene_dt <- getGeneExp(monet_class)[gn_lst][, -"gene"]
    atac_dt <- getAtacSeq(monet_class)[gn_lst][, -"gene"]
    gene_init <- getGeneInit(monet_class)[gn_lst][[2]]

    # add monetData object
    par_list[["no_gns"]] <- length(gn_lst)
    par_list[["no_tpt"]] <- ncol(atac_dt)
    par_list[["a"]] <- t(as.matrix(atac_dt))

    par_list[["y"]] <- t(as.matrix(gene_dt))
    par_list[["y0"]] <- gene_init

    return(par_list)

}
