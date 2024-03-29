#' Fine fitting functions
#'
#' double logistics, piecewise logistics and many other functions to
#' curve fit VI time-series.
#' 
#' * `Logistic` The traditional simplest logistic function. It can
#'      be only used in half growing season, i.e. vegetation green-up or senescence
#'      period.
#' * `doubleLog.Zhang` Piecewise logistics, (Zhang Xiaoyang, RSE, 2003).
#' * `doubleAG` Asymmetric Gaussian.
#' * `doubleLog.Beck` Beck logistics.
#' * `doubleLog.Gu` Gu logistics.
#' * `doubleLog.Elmore` Elmore logistics.
#' * `doubleLog.Klos` Klos logistics.
#' 
#' All of those function have `par` and `formula` attributes for the
#' convenience for analytical D1 and D2
#' 
#' @param par A vector of parameters
#' @param t A `Date` or numeric vector
#' 
#' @references
#' 1. Beck, P.S.A., Atzberger, C., Hogda, K.A., Johansen, B., Skidmore, A.K.,
#'      2006. Improved monitoring of vegetation dynamics at very high latitudes:
#'      A new method using MODIS NDVI. Remote Sens. Environ.
#'      https://doi.org/10.1016/j.rse.2005.10.021.
#' 2. Elmore, A.J., Guinn, S.M., Minsley, B.J., Richardson, A.D., 2012.
#'      Landscape controls on the timing of spring, autumn, and growing season
#'      length in mid-Atlantic forests. Glob. Chang. Biol. 18, 656-674.
#'      https://doi.org/10.1111/j.1365-2486.2011.02521.x. \cr
#'
#' 3. Gu, L., Post, W.M., Baldocchi, D.D., Black, TRUE.A., Suyker, A.E., Verma,
#'      S.B., Vesala, TRUE., Wofsy, S.C., 2009. Characterizing the Seasonal Dynamics
#'      of Plant Community Photosynthesis Across a Range of Vegetation Types,
#'      in: Noormets, A. (Ed.), Phenology of Ecosystem Processes: Applications
#'      in Global Change Research. Springer New York, New York, NY, pp. 35-58.
#'      https://doi.org/10.1007/978-1-4419-0026-5_2. \cr
#' 4. Peter M. Atkinson, et al., 2012, RSE, 123:400-417
#'
#' 5. https://github.com/cran/phenopix/blob/master/R/FitDoubleLogGu.R
#' @example R/examples/ex-FitDL.R
#' @rdname logistics
#' @export
Logistic = function(par, t){
    mn   = par[1]
    mx   = par[2]
    sos  = par[3]
    rsp  = par[4]
    pred = (mx - mn)/(1 + exp(-rsp*(t - sos))) + mn
    # pred = c/(1 + exp(a + b * t)) + d
    return(pred)
}
attr(Logistic, 'par')     = c("mn", "mx", "sos", "rsp")
attr(Logistic, 'formula') = expression((mx - mn)/(1 + exp(-rsp*(t - sos))) + mn)

# piecewise function
#' @rdname logistics
#' @export
doubleLog.Zhang = function(par, t){
    t0  = par[1]
    mn  = par[2]
    mx  = par[3]
    sos = par[4]
    rsp = par[5]
    eos = par[6]
    rau = par[7]

    if (t0 - sos <= 1 || t0 - eos >= -1) return(rep(99.0, length(t)))
    # In order to make sure the shape of S curve, should be satisfy:
    # t0 < eos, t0 > sos

    # xpred1 = (mx - mn)/(1 + exp(-rsp*(t - sos))) + mn
    # xpred2 = (mx - mn)/(1 + exp(rau*(t - eos))) + mn
    # pred  = xpred1*(t <= t0) + xpred2*(t > t0)

    # above expressions cost 1.5 times as the below
    pred = (mx - mn)/c(1 + exp(-rsp*(t[t <= t0] - sos)),
                         1 + exp( rau*(t[t >  t0] - eos))) + mn
    return(pred)
}
attr(doubleLog.Zhang, 'par')     = c("t0", "mn", "mx", "sos", "rsp", "eos", "rau")
attr(doubleLog.Zhang, 'formula') = expression( mn + (mx - mn)/(1 + exp(-rsp*(t - sos))),
                                               mn + (mx - mn)/(1 + exp( rau*(t - eos))) )

#' @rdname logistics
#' @export
doubleLog.AG = function(par, t){
    t0  = par[1]
    mn  = par[2]
    mx  = par[3]
    rsp = par[4]
    a3  = par[5]
    rau = par[6]
    a5  = par[7]

    pred = mn + (mx - mn)*exp(- c( 
        ((t0 - t[t <= t0])*rsp) ^a3,
        ((t[t >  t0] - t0)*rau) ^a5) )
    return(pred)
}
# a3, a5 should be greater than 1
attr(doubleLog.AG, 'par')     = c("t0", "mn", "mx", "rsp", "a3", "rau", "a5")
attr(doubleLog.AG, 'formula') = expression( mn + (mx - mn)*exp(- ((t0 - t)*rsp) ^a3 ),
                                         mn + (mx - mn)*exp(- ((t - t0)*rau) ^a5 ))

#' @rdname logistics
#' @export
doubleLog.AG2 = function(par, t){
    t0   = par[1]
    mn_l = par[2]
    mn_r = par[3] 
    mx   = par[4]
    rsp  = par[5]
    a3   = par[6]
    rau  = par[7]
    a5   = par[8]

    t_l = t[t <= t0]
    t_r = t[t > t0]
    
    pred_l = mn_l + (mx - mn_l)*exp(- ((t0 - t_l)*rsp) ^a3)
    pred_r = mn_r + (mx - mn_r)*exp(- ((t_r - t0)*rau) ^a5)
    # pred = mn + (mx - mn)*exp(- c( ((t0 - t[t <= t0])*rsp) ^a3,
                    # ((t[t >  t0] - t0)*rau) ^a5) )
    # return(pred)
    c(pred_l, pred_r)
}
# a3, a5 should be greater than 1
attr(doubleLog.AG2, 'par')     = c("t0", "mn_l", "mn_r", "mx", "rsp", "a3", "rau", "a5")
attr(doubleLog.AG2, 'formula') = expression( 
    mn_l + (mx - mn_l)*exp(- ((t0 - t)*rsp) ^a3),
    mn_r + (mx - mn_r)*exp(- ((t - t0)*rau) ^a5)
)

#' @rdname logistics
#' @export
doubleLog.Beck = function(par, t) {
    mn  = par[1]
    mx  = par[2]
    sos = par[3]
    rsp = par[4]
    eos = par[5]
    rau = par[6]
    # if (sos >= eos) return(rep(9999, length(t)))
    # if (eos < sos) return(rep(99.0, length(t)))
    if (!all(is.finite(par)) || eos < sos) return(rep(99.0, length(t)))

    pred = mn + (mx - mn)*(1/(1 + exp(-rsp*(t - sos))) + 1/(1 + exp(rau*(t - eos))) - 1)
    return(pred)
}
attr(doubleLog.Beck, 'par')  = c("mn", "mx", "sos", "rsp", "eos", "rau")
attr(doubleLog.Beck, "formula") = expression(mn + (mx - mn) * (1 / (1 + exp(-rsp * (t - sos))) + 1 / (1 + exp(rau * (t - eos))) - 1))

#' @rdname logistics
#' @export
doubleLog.Elmore = function(par, t) {
    mn  = par[1]
    mx  = par[2]
    # m3  = par[3]
    # m4  = par[4]
    # m5  = par[5]
    # m6  = par[6]
    # m3l = m3/m4
    # m4l = 1/m4
    # m5l = m5/m6
    # m6l = 1/m6
    sos  = par[3] # SOS
    rsp  = par[4] # 1/rsp
    eos  = par[5] # EOS
    rau  = par[6] # 1/rau
    m7   = par[7]
    pred = mn + (mx - m7*t)*( 1/(1 + exp(-rsp*(t-sos))) - 1/(1 + exp(-rau*(t-eos))) )
    return(pred)
}
attr(doubleLog.Elmore, 'par')     = c("mn", "mx", "sos", "rsp", "eos", "rau", "m7")
attr(doubleLog.Elmore, 'formula') = expression( mn + (mx - m7*t)*( 1/(1 + exp(-rsp*(t-sos))) - 1/(1 + exp(-rau*(t-eos))) ) )

#' @rdname logistics
#' @export
doubleLog.Gu = function(par, t) {
    y0  = par[1]
    a1  = par[2]
    a2  = par[3]
    sos = par[4]
    rsp = par[5]
    eos = par[6]
    rau = par[7]
    # t1  = par[4]
    # t2  = par[5]
    # b1  = par[6]
    # b2  = par[7]
    c1  = par[8]
    c2  = par[9]
    # pred = y0 + (a1/(1 + exp(-(t - t1)/b1))^c1) - (a2/(1 + exp(-(t - t2)/b2))^c2)
    pred = y0 + (a1/(1 + exp(-rsp*(t - sos)))^c1) - (a2/(1 + exp(-rau*(t - eos)))^c2)
    return(pred)
}
attr(doubleLog.Gu, 'par')     = c('y0', 'a1', 'a2', 'sos', 'rsp', 'eos', 'rau', 'c1', 'c2')
attr(doubleLog.Gu, 'formula') = expression(y0 + (a1/(1 + exp(-rsp*(t - sos)))^c1) - (a2/(1 + exp(-rau*(t - eos)))^c2))

#' @rdname logistics
#' @export
doubleLog.Klos = function(par, t) {
    a1 = par[1]
    a2 = par[2]
    b1 = par[3]
    b2 = par[4]
    c  = par[5]
    B1 = par[6]
    B2 = par[7]
    m1 = par[8]
    m2 = par[9]
    q1 = par[10]
    q2 = par[11]
    v1 = par[12]
    v2 = par[13]
    pred = (a1*t + b1) + (a2*t^2 + b2*t + c) * 
        (1/(1 + q1 * exp(-B1 * (t - m1)))^v1 - 
         1/(1 + q2 * exp(-B2 * (t - m2)))^v2 )
    return(pred)
}
attr(doubleLog.Klos, 'par')     = c('a1', 'a2', 'b1', 'b2', 'c', 'B1', 'B2',
    'm1', 'm2', 'q1', 'q2', 'v1', 'v2')
attr(doubleLog.Klos, 'formula') = expression((a1*t + b1) + (a2*t^2 + b2*t + c) * (1/(1 + q1 * exp(-B1 * (t - m1)))^v1
        - 1/(1 + q2 * exp(-B2 * (t - m2)))^v2))

.qr.solve = function(a, b, tol = 1e-07, LAPACK = TRUE) {
    if (!is.qr(a)) a = qr(a, tol = tol, LAPACK = LAPACK)
    nc = ncol(a$qr)
    nr = nrow(a$qr)
    if (a$rank != min(nc, nr)) stop("singular matrix 'a' in solve")
    if (missing(b)) {
        if (nc != nr) stop("only square matrices can be inverted")
        b = diag(1, nc)
    }
    res = qr.coef(a, b)
    res[is.na(res)] = 0
    res
}
# vc = .qr.solve(opt$hessian)
# npar = nrow(vc)
# s2 = opt.df$cost[best]^2 / (n - npar)
# std.errors = sqrt(diag(vc) * s2)     # standard errors

# attach gradient and hessian analytical function to curve fitting functions
.dls = lapply(
    c("doubleLog.Beck", "doubleLog.Elmore", "doubleLog.Gu",
      "doubleLog.Klos", "doubleLog.Zhang", "doubleLog.AG"),
    function (FUN){
        # FUN = deparse(substitute(fun))
        fun = get(FUN)
        attr(fun, 'gradient') = gradf_t(fun) # gradient
        attr(fun, 'hessian')  = hessf_t(fun) # hessian
        # print(environment(fun))
        assign(FUN, fun, envir = environment(fun)) #environment("namespace:phenofit"))#
        # fun
    })
