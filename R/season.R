#' @title Growing season dividing
#' @name season
#'
#' @description
#' Divide growing seasons according to rough fitting (`rFUN`) result .
#'
#' For `season`, rough fitting is applied for whole.
#' For `season_mov` rough fitting is applied in every year, during which
#' `maxExtendMonth` is extended.
#'
#' @details
#' Before dividing growing season, `INPUT` should be added a year in head
#' and tail first by `add_HeadTail`.
#'
#' Finally, use [findpeaks()] to get local maximum and local minimum values.
#' Two local minimum define a growing season.
#' If two local minimum(maximum) are too closed, then only the smaller(biger) is left.
#'
#' @param INPUT A list object with the elements of `t`, `y`, `w`,
#' `Tn` (optional) and `ylu`, returned by [check_input()].
#' @param rFUN Rough curve fitting function, can be one of [wSG()],
#' [wWHIT()] and [wHANTS()].
#' @param wFUN weights updating function, can be one of [wTSM()],
#' [wChen()], [wBisquare()] and [wSELF()].
#' @param iters How many times curve fitting is implemented.
#' @param wmin Double, minimum weigth (i.e. weight of snow, ice and cloud).
#' @param lambda The smoothing parameter of [wWHIT()]. For
#' [season_mov()], if lambda is `NULL`, [init_lambda()]
#' will be used. Generally, it was set as 10000, 15, and 5 for daily, 8-day
#' and 16-day inputs respectively.
#' @param nf The parameter of [wHANTS()], number of frequencies to be
#' considered above the zero frequency.
#' @param frame The parameter of [wSG()], moving window size. Suggested by
#' TIMESAT, default `frame = floor(nptperyear/7)*2 + 1`.
#' @param minpeakdistance Numberic, in the unit of points (default as
#' `nptperyear/6`). The minimum distance of two peaks. If the distance of two
#' maximum extreme value less than `minpeakdistance`, only the real maximum
#' value will be left.
#' @param r_min Threshold is defined as the difference of peak value with
#' trough value. There are two threshold (left and right). The minimum threshold
#' should be greater than r_min.
#' @param r_max Similar as `r_min`, The maximum threshold should
#' be greater than `r_max`.
#' @param ypeak_min ypeak >= ypeak_min
#' @param rtrough_max `ytrough <= rtrough_max*A`, A is the amplitude of y.
#' @param MaxPeaksPerYear This parameter is used to adjust lambda in iterations.
#' If PeaksPerYear > MaxPeaksPerYear, then lambda = lambda*2.
#' @param MaxTroughsPerYear This parameter is used to adjust lambda in iterations.
#' If TroughsPerYear > MaxTroughsPerYear, then lambda = lambda*2.
#' @param calendarYear If true, only one static calendar growing season will be
#' returned.
#' @param IsPlot Boolean
#' @param plotdat (optional) A list or data.table, with `t`, `y` and `w`.
#' Only if `IsPlot=TRUE`, [plot_input()] will be used to plot.
#' Known that y and w in `INPUT` have been changed, we suggest using the
#' original data.table.
#' @param print Whether to print progress information
#' @param adj.param Adjust rough curve fitting function parameters automatically,
#' if too many or to less peak and trough values.
#' @param ... For [season_mov()], Other parameters passed to
#' [season()]; For [season()], other parameters passed to
#' [findpeaks()].
#'
#' @return
#' * `whit`: Rough fitting result.
#' * `dt`: Growing season dividing information.
#'
#' @seealso [findpeaks()].
#' @export
# $whit
# # A tibble: 574 x 5
#    t              y     w iter1 iter2
#    <date>     <dbl> <dbl> <dbl> <dbl>
#  1 2002-07-04  9.43 0.149 10.1  10.5
#   $dt
# # A tibble: 11 x 7
#    beg        peak       end        len     year season flag
#    <date>     <date>     <date>     <time> <dbl>  <int> <chr>
#  1 2003-01-09 2003-07-12 2004-02-18 406     2003      1 2003_1
season <- function(INPUT,
                   rFUN = wWHIT, wFUN = wTSM, iters = 2, wmin = 0.1,
                   lambda, nf  = 3, frame = floor(INPUT$nptperyear/5)*2 + 1,
                   minpeakdistance,
                   r_max = 0.2, r_min = 0.05,
                   ypeak_min   = 0.1, rtrough_max = 0.6,
                   MaxPeaksPerYear = 2, MaxTroughsPerYear = 3,
                   calendarYear = FALSE,
                   IsPlot  = FALSE, plotdat = INPUT, print = FALSE,
                   adj.param = TRUE,
                   ...)
{
    nptperyear <- INPUT$nptperyear
    south      <- INPUT$south
    t          <- INPUT$t
    y          <- INPUT$y

    nlen       <- length(t)

    # 1. How many years data
    date_year <- year(t) + ((month(t) >= 7)-1)*south
    info  <- table(date_year) # rm years with so limited obs
    years <- info[info > nptperyear*0.2] %>% {as.numeric(names(.))}
        #.[2:(length(.)-1)] # rm head and tail filled years
    # nyear     <- length(years)

    if (missing(frame)) frame <- floor(nptperyear/5) * 2 + 1
    if (missing(minpeakdistance)) minpeakdistance <- nptperyear/6

    if (all(is.na(y))) return(NULL)
    npt   <- sum(INPUT$w > wmin)
    nyear <- round(npt/nptperyear) # Needn't to ceil, matter most for parameter
                                   # adjustment

    ylu0  <- INPUT$ylu
    ylu   <- ylu0
    A0    <- diff(ylu0)
    # frame <- floor(nptperyear/7) * 2 + 1 #13, reference by TSM SG filter
    if (missing(lambda)) lambda <- max(nyear*frame, 15)

    ## 1. weighted curve fitting help to divide growing season
    iloop <- 1
    while (iloop <= 3){
        # sgfitw(INPUT$y, INPUT$w, nptperyear, INPUT$ylu, wFUN, iters, frame, d=2)
        # whitsmw(y, w, ylu, wFUN, iters = 1, lambda = 100, ..., d = 2, missval)
        param <- c(INPUT,
            wFUN = wFUN, wmin = wmin, iters = iters,
            lambda = lambda,  # param for whittaker
            nf     = nf,      # param for wHANTS,
            frame  = frame    # param for wSG
        )

        yfits <- do.call(rFUN, param)
        zs    <- yfits$zs
        ypred <- last(zs) #as.numeric(runmed(ypred, frame))
        alpha <- 0.01

        # default is three year data input, median will be much better
        # This module was primarily designed for `season_mov`. It also works for
        # group large than 3-year. 2-year median will be underestimated.
        if (nyear >= 2.5){ # considering NA values, nyear of 3-year will be smaller.

            ylu_min <- aggregate(ypred, list(year = year(t)), min)$x %>% median()
            ylu_max <- aggregate(ypred, list(year = year(t)), max)$x %>% median()

            # If multiple years are completely missing, ylu_min possiable equals ylu_max
            if (ylu_max - ylu_min > A0*0.2){
                ylu <- c(pmax(ylu_min, INPUT$ylu[1]), #quantile(ypred, alpha/2)
                         pmin(ylu_max, INPUT$ylu[2]))
            }
        }
        INPUT$ylu <- ylu
        A         <- diff(ylu)

        # ylu   <- quantile(ypred, c(alpha/2, 1 - alpha), na.rm = TRUE)
        ## Plateau peak will lead to failed to find local extreme values.
        #  To avoid fluctuating in peak of growing season or flat peak growing
        #  season, like fluxsite: ZM-Mon
        # max_slp <- 2*A/nptperyear
        # pek_slp <- abs(coefficients(lm(ypred[I]~I))[[2]])
        #
        # I <- which(ypred > (0.8*A + ylu[1]))
        # if (length(I)/length(y) > 0.3){
        #     ypred[I] <- median(ypred[I])
        # }
        #
        # local minimum values
        # peak values is small for minimum values, so can't use r_min here
        peaks <- findpeaks(-ypred,
                           r_max = r_max*A,
                           r_min = r_min*0,
                           minpeakdistance = minpeakdistance, zero = "-", nups = 0)
        pos_min   <- peaks$X
        pos_min[, 1] %<>% multiply_by(-1)
        ntrough_PerYear <- length(peaks$gregexpr)/nyear#max peaks
        # local maximum values,
        peaks <- findpeaks(ypred, zero = "+",
                           r_max = r_max*A,
                           r_min = r_min*A, # A
                           minpeakdistance = minpeakdistance,
                           minpeakheight = 0.1*A + ylu[1], nups = 1) #, ypeak_min
        pos_max <- peaks$X
        npeak_PerYear <- length(peaks$gregexpr)/nyear#max peaks

        if (print)
            cat(sprintf('iloop = %d: lambda = %.1f, ntrough_PerYear = %.2f, npeak_PerYear = %.2f\n',
                iloop, lambda, ntrough_PerYear, npeak_PerYear))
        # maxpeaksperyear <- 2

        ## This module will automatically update lambda, nf and wHANTS
        #  Not only wWHd, it has been extended to wHANT and wSG. 20180910

        if (adj.param){
            delta_frame <- ceiling(nptperyear/12) # adjust frame in the step of `month`
            if (npeak_PerYear > MaxPeaksPerYear | ntrough_PerYear > MaxTroughsPerYear){
                lambda <- lambda*2
                nf     <- max(1, nf - 1)
                frame  <- min(frame + delta_frame, nptperyear*2)
            }else if (npeak_PerYear < 0.8  | ntrough_PerYear < 0.8){
                lambda <- lambda/2
                nf     <- min(5, nf + 1)
                frame  <- max(frame - delta_frame, delta_frame)
            }else{
                break
            }
            iloop <- iloop + 1
        } else {
            break
        }
    }

    ## Prepare raw OUTPUT for error condition
    #  rough curve fitting time-series
    rfit = as.data.table(c(list(t = t, y = y), yfits$ws, yfits$zs))

    # 1.1 the local minimum value should small than rtrough_max*A
    if (!is.null(pos_min)) {
        pos_min <- pos_min[(val - ylu[1]) <= rtrough_max*A, ]
        pos_min[, type := -1]
    }
    if (!is.null(pos_max)) pos_max[, type :=  1]

    pos <- rbind(pos_min, pos_max) # dtaa.table returned
    if (!is.null(pos)) pos <- pos[order(pos), ] #c("val", "pos", "left", "right", "type")

    dt  <- di <- NULL
    res <- list(whit = rfit, dt = dt) # , pos = pos, di = di

    if (is.null(pos_max) || is.null(pos_min)){
        warning("Can't find a complete growing season before trim!")
        return(res)
    }
    # plot(ypred, type = "b"); grid()

    # rm peak value if peak value smaller than the nearest trough values
    I   <- !with(pos, (c(diff(val) > 0, FALSE) & c(diff(type) == -2, FALSE)) |
            (c(FALSE, diff(val) < 0) & c(FALSE, diff(type) == 2)))
    # whether remove y_peak first is important
    # I <- which(dt$y_peak >= ypeak_min)
    pos <- pos[I, ]


    # 1.2 remove both points (date or value of min and max too close)
    # I_del <- union(I_date, I_date+1)
    # I_del <- union(I_date + 1, I_val + 1)
    for (i in 1:2){
        if (i == 1){
            # remove dates too close first, then rm continue max or min value
            I_del  <- which(diff(pos$pos) < minpeakdistance/2) + 1 # for date
        } else if(i == 2){
            # remove values too close second, then rm continue max or min value again.
            I_del  <- which(abs(diff(pos$val)) < 0.05*A) + 1 #for value
        }
        if (length(I_del) > 0) pos <- pos[-I_del, ]
        # 1.3 remove replicated
        pos$flag <- cumsum(c(1, diff(pos$type) != 0))
        pos      <- ddply(pos, .(flag), rm_duplicate, y = ypred, threshold = r_min*A)[, 2:6]
    }

    pos$t    <- t[pos$pos]
    res$pos <- pos

    # print(nrow(pos))
    if (nrow(pos) < 2){ # at least two points, begin and end
        warning("Can't find a complete growing season before!")
        return(res)
    }

    ############################################################################
    ## 5. check head and tail break point, and reform breaks
    locals <- pos[, c("pos", "type")]
    ns     <- nrow(locals)
    # check the head and tail minimum values
    minlen <- nptperyear/3 #distance from peak point
    if (last(pos$type) == 1 && (nlen - nth(pos$pos, -2)) > minlen &&
        abs(last(ypred) - nth(pos$val, -2)) < 0.15*A )
        locals %<>% rbind.data.frame(., data.frame(pos = nlen, type = -1))
    if (pos$type[1] == 1 && pos$pos[2] > minlen && abs(ypred[1] - pos$val[2]) < 0.15*A)
        locals %<>% rbind.data.frame(data.frame(pos = 1, type = -1), .)

    # a complete growing season, from minimum to minimum
    I      <- which(locals$type == -1)
    locals <- locals[I[1]:I[length(I)], ]

    s  <- locals$pos; ns <- length(s)
    if (ns < 3) {
        warning("Can't find a complete growing season!")
        return(res)
    }

    locals %<>% mutate(val = ypred[pos], t = t[pos])

    pos_max <- subset(locals, type == 1)
    pos_min <- subset(locals, type == -1)
    ## 6. divide into multiple growing seasons
    di <- data.table(beg  = s[seq(1, ns-1, 2)],
                     peak = s[seq(2, ns, 2)],
                     end  = s[seq(3, ns, 2)])
    di %<>% fix_di(t = t) #fix whole year data missing

    dt <- map(di, ~t[.x]) %>% as.data.table() %>%
        .[, `:=`( y_beg  = ypred[di$beg],
                  y_peak = ypred[di$peak],
                  y_end  = ypred[di$end],
                  len    = as.integer(difftime(end, beg, units = "days") + 1),
                  year   = year(peak) )]

    # update 20180913
    # solve the problem of growing season too long, (e.g. US-Cop).
    I <- which(dt$y_peak >= ypeak_min)
    dt <- dt[I, ]
    # di <- di[I, ]
    # res$di <- di

    # remove di and add fix_dt
    if (calendarYear) {
        # need to remove incomplete year
        dt <- season_calendar(years, south) # [2:(nyear-1)]
    } else {
        dt <- dt[dt$len > 45 & dt$len < 650, ] # mask too long and short gs
        # phenofit:::fix_dt(dt) # c++ address operation, fix growing season overlap
        fix_dt(dt) # c++ address operation, fix growing season overlap
        # after fix_dt, growing season length will become shorter
        dt <- dt[dt$len > 45 & dt$len < 650, ] # mask too long and short gs
    }
     # get the growing season year, not only the calendar year
    if (south) dt[, year := year + as.integer(peak >= ymd(sprintf('%d0701', year))) - 1L]

    dt[, `:=`(season = as.numeric(1:.N), flag = sprintf("%d_%d", year, 1:.N)), .(year)]
    res <- list(whit = rfit, dt = dt)

    ## 7. plot
    if (IsPlot) plot_season(INPUT, res, plotdat, INPUT$ylu)
    return(res)
}

# rm duplicated max or min values
rm_duplicate <- function(d, y, threshold){
    d <- d[, 1:5]
    if (nrow(d) > 1){
        type <- d$type[1] #1(-1) represent max(min)
        # if range amplitude less than TRS, get median
        if (diff(range(d$val)) < threshold){
            I <- floor(median(d$pos))
            tibble(val = y[I], pos = I,
                       left = min(d$left),
                       right = max(d$right), type = type)
        }else{
            # else, get the local extreme value
            fun <- ifelse(type == 1, which.max, which.min)
            d[fun(d$val), ]
        }
    }else{
        d
    }
}

# fix across multi-year breaks points, when whole year data are missing
#
# This function is only for fluxsites data
fix_di <- function(di, t){
    for (i in 1:nrow(di)){
        I    <- di$beg[i]:di$end[i]
        # Try to remove NA values at head and tail. Fialed, Na value may not
        # at head or tail.
        # I_nona <- which(!is.na(y[I])) %>% {I[first(.):last(.)]}
        I_nona <- I#checking year brk is enough
        ti     <- t[I_nona]

        # check year brokens
        I_brkyear <- which(diff(ti) >= 365)
        nbrk      <- length(I_brkyear)

        if (nbrk > 0 & nbrk <= 2){
            if (nbrk == 1) {
                I_1 <- I_nona[1:I_brkyear]
                I_2 <- I_nona[(I_brkyear + 1):length(I_nona)]

                lst <- list(I_1, I_2)
            } else if (nbrk == 2) {
                I_1 <- I_nona[1:I_brkyear[1]]
                I_2 <- I_nona[(I_brkyear[1]+1):I_brkyear[2]]
                I_3 <- I_nona[(I_brkyear[2]+1):length(I_nona)]

                lst <- list(I_1, I_2, I_3)
            }
            #select the longest segment
            I_nona <- lst[[which.max(sapply(lst, length))]]
            di$beg[i] <- first(I_nona)
            di$end[i] <- last(I_nona)
        }
    }
    di#quickly return
}