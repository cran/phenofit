## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7, fig.height=5
)

## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set() 

## ----message=F----------------------------------------------------------------
library(phenofit)
library(data.table)
library(dplyr)
library(ggplot2)

## -----------------------------------------------------------------------------
d = MOD13A1$dt %>% subset(site == "CA-NS6" & date >= "2010-01-01" & date <= "2016-12-31") %>%
    .[, .(date, y = EVI/1e4, DayOfYear, QC = SummaryQA)]
d %<>% mutate(t = getRealDate(date, DayOfYear)) %>%
    cbind(d[, as.list(qc_summary(QC, wmin = 0.2, wmid = 0.5, wmax = 0.8))]) %>%
    .[, .(date, t, y, QC_flag, w)]
print(d)

## -----------------------------------------------------------------------------
lambda         <- 8
nptperyear     <- 23
minExtendMonth <- 0.5
maxExtendMonth <- 1
minPercValid   <- 0
wFUN           <- wTSM # wBisquare
wmin           <- 0.2
methods_fine <- c("AG", "Zhang", "Beck", "Elmore", "Gu")

## -----------------------------------------------------------------------------
INPUT <- check_input(d$t, d$y, d$w,
    QC_flag = d$QC_flag,
    nptperyear = nptperyear,
    maxgap = nptperyear / 4, wmin = 0.2
)

brks <- season_mov(INPUT,
    list(FUN = "smooth_wWHIT", wFUN = wFUN,
        maxExtendMonth = 3,
        wmin = wmin, r_min = 0.1
    ))
# plot_season(INPUT, brks)

## 2.4 Curve fitting
fit <- curvefits(INPUT, brks,
    list(
        methods = methods_fine, # ,"klos",, 'Gu'
        wFUN = wFUN,
        iters = 2,
        wmin = wmin,
        # constrain = FALSE,
        nextend = 2,
        maxExtendMonth = maxExtendMonth, minExtendMonth = minExtendMonth,
        minPercValid = minPercValid
    ))

## check the curve fitting parameters
l_param <- get_param(fit)
print(l_param$Beck)

dfit <- get_fitting(fit)
print(dfit)

## 2.5 Extract phenology
TRS <- c(0.1, 0.2, 0.5)
l_pheno <- get_pheno(fit, TRS = TRS, IsPlot = FALSE) # %>% map(~melt_list(., "meth"))
print(l_pheno$doy$Beck)

pheno <- l_pheno$doy %>% melt_list("meth")

## -----------------------------------------------------------------------------
# growing season dividing
plot_season(INPUT, brks, ylab = "EVI")
# Ipaper::write_fig({  }, "Figure4_seasons.pdf", 9, 4)

# fine curvefitting
g <- plot_curvefits(dfit, brks, title = NULL, cex = 1.5, ylab = "EVI", angle = 0)
grid::grid.newpage()
grid::grid.draw(g)
# Ipaper::write_fig(g, "Figure5_curvefitting.pdf", 8, 6, show = TRUE)

# extract phenology metrics, only the first 3 year showed at here
# write_fig({
l_pheno <- get_pheno(fit[1:7], method = "AG", TRS = TRS, IsPlot = TRUE, show.title = FALSE)
# }, "Figure6_phenology_metrics.pdf", 8, 6, show = TRUE)

## -----------------------------------------------------------------------------
# library(ggplot2)
# library(ggnewscale)

# # on the top of `Figure7_predata...`
# d_comp = fread("data-raw/dat_Figure7_comparison_with_others-CA-NS6.csv")
# d_comp = merge(d[, .(date, t)], d_comp[, .(date, TIMESAT, phenopix)]) %>%
#     merge(dfit[meth == "Beck", .(t, phenofit = ziter2)], by = "t") %>%
#     melt(c("date", "t"), variable.name = "meth")

# labels = c("good", "marginal", "snow", "cloud")
# theme_set(theme_grey(base_size = 16))
# cols_line = c(phenofit = "red", TIMESAT = "blue", phenopix = "black")
# p <- ggplot(dfit, aes(t, y)) +
#     geom_point(aes(color = QC_flag, fill = QC_flag, shape = QC_flag), size = 3) +
#     scale_shape_manual(values = qc_shapes[labels], guide = guide_legend(order = 1)) +
#     scale_color_manual(values = qc_colors[labels], guide = guide_legend(order = 1)) +
#     scale_fill_manual(values = qc_colors[labels], guide = guide_legend(order = 1)) +
#     new_scale_color() +
#     geom_line(data = d_comp, aes(t, value, color = meth)) +
#     # geom_line(data = d_comp[meth == "phenofit"], aes(t, value),
#     #           size = 1, show.legend = FALSE, color = "red") +
#     scale_color_manual(values = cols_line, guide = guide_legend(order = 2)) +
#     labs(x = "Time", y = "EVI") +
#     theme(
#         axis.title.x = element_text(margin = margin(t = 0, unit='cm')),
#         # plot.margin = margin(t = 0, unit='cm'),
#         legend.text = element_text(size = 13),
#         legend.position = "bottom",
#             legend.title  = element_blank(),
#             legend.margin = margin(t = -0.3, unit='cm'))
# p
# # write_fig(p, "Figure7_comparison_with_others.pdf", 10, 4, show = TRUE)

