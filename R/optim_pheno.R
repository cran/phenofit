#' optim_pheno
#'
#' Interface of optimization functions for double logistics and other parametric
#' curve fitting functions.
#'
#' @inheritParams smooth_wHANTS
#' @param prior A vector of initial values for the parameters for which optimal
#' values are to be found. `prior` is suggested giving a column name.
#' @param sFUN The name of fine curve fitting functions, can be one of `
#' 'FitAG', 'FitDL.Beck', 'FitDL.Elmore', 'FitDL.Gu' and 'FitDL.Klos',
#' 'FitDL.Zhang'`.
#' @param nptperyear Integer, number of images per year, passed to `wFUN`.
#' Only [wTSM()] needs `nptperyear`. If not specified,
#' `nptperyear` will be calculated based on `t`.
#' @param ylu `ymin, ymax`, which is used to force `ypred` in the
#' range of `ylu`.
#' @param tout Corresponding doy of prediction.
#' @param method The name of optimization method to solve fine fitting, one of
#' `'BFGS','CG','Nelder-Mead', 'L-BFGS-B', 'nlm', 'nlminb', 'ucminf'` and
#' `'spg','Rcgmin','Rvmmin', 'newuoa','bobyqa','nmkb','hjkb'`.
#'
#' @param verbose Whether to display intermediate variables?
#' @param ... other parameters passed to [I_optim()] or [I_optimx()].
#'
#' @return fFIT object, see [fFIT()] for details.
#'
#' @seealso [FitDL()]
#'
#' @example man/examples/ex-optim_pheno.R
#'
#' @import optimx
#' @export
optim_pheno <- function(
    prior, sFUN,
    y, t, tout, method,
    w, nptperyear, ylu,
    iters = 2, wFUN = wTSM,
    verbose = FALSE, ...)
{
    sFUN = gsub("\\.", "_", sFUN )
    FUN <- get(sFUN, mode = "function" )

    ntime = length(t)
    # add prior colnames
    parnames <- attr(FUN, 'par')
    if (is.vector(prior)) prior <- t(prior)
    colnames(prior) <- parnames
    npar = ncol(prior)

    # opt.df: par, value, fevals, niter, convcode
    J_VALUE    = npar + 1
    J_CONVCODE = npar + 4

    methods_optim  <- c('BFGS','CG','Nelder-Mead', 'L-BFGS-B', 'nlm', 'nlminb', 'ucminf')
    methods_optimx <- c('spg','Rcgmin','Rvmmin', 'newuoa','bobyqa','nmkb','hjkb')

    if (method %in% methods_optim){
        I_optimFUN <- I_optim
    } else if (method %in% methods_optimx){
        I_optimFUN <- I_optimx
    } else {
        stop(sprintf('optimization method (%s) is not supported!', method))
    }
    # To improve the performance of optimization, `t` needs to be normalized.
    tout_0 <- tout

    # check input parameters
    if (missing(w))          w   <- rep(1, length(y))
    if (missing(nptperyear)) nptperyear <- ceiling(365/mean(as.numeric(diff(t))))
    if (missing(ylu))        ylu <- range(y, na.rm = TRUE)
    A <- diff(ylu) # y amplitude

    fits <- list()
    ws   <- list() #record weights for every iteration
    ## 1. add weights updating methods here
    par   <- setNames(numeric(length(parnames))*NA, parnames)
    ypred <- rep(NA, length(tout))

    if (verbose){
        boundary <- list(...)[c('lower', 'upper')] %>% do.call(rbind, .) %>%
            set_colnames(parnames) %>% as_tibble()
        print(boundary)
    }

    for (i in 1:iters){
        ws[[i]] <- w
        # dot = list(...)
        # params = c(list(prior, FUN, y, t, method = method, w = w, verbose = verbose), dot[-4], use.julia = TRUE)
        # do.call(I_optimFUN, params)
        # save(list = ls(), file = "debug.rda")
        # pass verbose for optimx optimization methods selection
        opt.df  <- I_optimFUN(prior, FUN, y, t, method = method, w = w, verbose = verbose, ...)
        # print(opt.df)
        if (verbose){
            fprintf('Initial parameters:\n')
            print(as_tibble(prior))
            print(opt.df)
        }

        best <- which.min(opt.df[, J_VALUE])
        if (is_empty(best)){
            warning("optimize parameters failed!")
        }else{
            # test for convergence if maximum iterations where reached - restart from best with more iterations
            # If RMSE is small enough, then accept it.
            if (opt.df[best, J_CONVCODE] != 0) { # convcode
                # repeat with more maximum iterations if it did not converge
                par <- opt.df[best, 1:npar, drop = FALSE]  #best par generally
                opt <- I_optimFUN(par, FUN, y, t, method, w = w, verbose = verbose, ...)
            } else { #if (opt.df$convcode[best] == 0)
                opt <- opt.df[best, , drop = FALSE]
            }

            RMSE = sqrt(opt[1, J_VALUE])/ntime # SSE returned
            # check whether convergence, only 1 row now
            if (opt[1, J_CONVCODE] != 0 & RMSE > 0.1*A) {
                warning("Not convergent!")
            }else{
                par   <- opt[1, 1:npar, drop = FALSE]
                # put opt par into prior for the next iteration
                prior <- rbind(prior[1:(nrow(prior)-1), ], par, deparse.level = 0) %>%
                    unique.matrix()
                FUN(par, tout, ypred)
                # too much missing values
                # if (sum(w == 0)/length(w) > 0.5) ypred <- ypred*NA
                # to adapt wTS, set iter = i-1; #20180910
                # nptperyear, wfact = 0.5)
                w <- tryCatch(
                    wFUN(y, FUN(par, t), w, i, nptperyear, ...),
                    error = function(e) {
                        message(sprintf('[%s]: %s', sFUN, e$message))
                        return(w) #return original w
                    })
                ypred %<>% check_ylu(ylu) #values out of ylu are set as NA
            }
        }
        fits[[i]] <- ypred
    }
    fits %<>% set_names(paste0('iter', 1:iters))
    ws   %<>% set_names(paste0('iter', 1:iters)) # %>% as_tibble()

    # 02. uncertain part also could add here
    # returned object FUN need to be futher optimized
    structure(list(tout = tout_0, zs = fits, ws = ws,
        par  = par, fun = sFUN), class = 'fFIT')
}

# 1. Firstly, using \pkg{optimx} package to select the best performance
#   optimization function.
# 2. Secondly, even though optimx is extraordinary, but it running is not so
#   satisfied, so using self unified optimization function (with prefixion of
#   'p_', opt_FUN means phenology optimization function) to replace it through
#   `optim_p` function. The unified procedure was inspired by `optimx` package.

# 1. nlfr package
# start <- set_names(prior[1, ], attr(FUN, "par"))
# # lower = lower, upper = upper
# f_nlf <- nlfb(start, resfn = function(par, t, y, fun, ...){
#     fun(par, t) - y
# }, trace = FALSE, y = y, t=t, w = w, fun = FUN, ...)
# yfit  <- FUN(f_nlf$coefficients, t)
# plot(t, y, type = "b")
# lines(t, yfit)
# print(f_nlf)
