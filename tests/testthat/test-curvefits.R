context("curvefits")

source('helper_MOD13A1.R')
wFUN = wTSM

# lambda   <- init_lambda(INPUT$y) # lambda for whittaker
# # param = listk(
# #     INPUT, nptperyear,
# #     FUN = whitsmw2, wFUN = wBisquare, iters = 2,
# #     lambda,
# #     IsPlot = IsPlot, plotdat = d,
# #     south = sp$lat[1] < 0,
# #     rymin_less = 0.6, ypeak_min = ypeak_min,
# #     max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5
# # )
brks2 <- season_mov(INPUT,
    options = list(rFUN = "smooth_wWHIT", wFUN = wFUN))

param <- list(
    INPUT, brks2,
    options = list(
        methods = c("AG", "Beck", "Elmore", "Gu", "Zhang"), #,"klos",
        wFUN = wFUN, nextend = 2, maxExtendMonth = 3, minExtendMonth = 1,
        minPercValid = 0.2, use.rough = TRUE, verbose = FALSE
    )
)

meth <- "AG"
test_curvefit <- function(meth){
    test_that(sprintf("`curvefits` with %s", meth), {
        expect_silent({
            suppressWarnings({
                param$options$methods <- meth
                fit  <- do.call(curvefits, param)
            })
        })
    })
}

test_curvefit("AG") # works with rough fitting result

param$use.rough <- FALSE
test_curvefit("AG")
test_curvefit("Beck")
test_curvefit("Elmore")
test_curvefit("Gu")
test_curvefit("Zhang")

# ## check the curve fitting parameters
# params <- getparam(fit)
# # print(str(params, 1))
# # print(params$AG)

# ## Get GOF information
# stat  <- ldply(fit$fits, function(fits_meth){
#     ldply(fits_meth, statistic.fFIT, .id = "flag")
# }, .id = "meth")
# fit$stat <- stat
# print(head(stat))
# # print(mean(stat$NSE, na.rm = T))
# ## visualization
# # svg("Figure1_phenofit_curve_fitting.svg", 11, 7)
# # Cairo::CairoPDF(file_pdf, 11, 6) #
# # dev.off()
# titlestr = "test logistics"
# g <- plot_phenofit(fit, d, titlestr)
# grid::grid.newpage(); grid::grid.draw(g)# plot to check the curve fitting
# test_that("multiplication works", {
#   expect_equal(2 * 2, 4)
# })
