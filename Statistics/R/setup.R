#' Global setup, variables and settings

here_statistics = function (...) here::here("Statistics", ...)
here_r = function (...) here_statistics("R", ...)

###############################################################################
# Plotting settings

ggplot2::theme_set(ggplot2::theme_classic())

###############################################################################
# Style

style = new.env()

#  4 categorical colors
# orange-red, green, orange, blue
style$pal4 = c("#D55E00", "#009E73", "#E69F00", "#56B4E9")
#  3 colors Red, green, grey
style$col_grey = "#999999"
style$col_trueneg = "#56B4E9"
style$col_truepos = "#D55E00"
style$pal3 <- c(style$col_grey, style$col_trueneg, style$col_truepos)

# Output formats
style$output_formats = c("pdf", "png")
style$dpi = 400

