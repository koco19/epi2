#---
#title: "OR Plot"
#author: "Mercè Garí"
#date: '2021-03-29'
#---

here_r = function (...) here::here("Statistics", ...)
here_output = function (...) here::here("Output", ...)

# Setup
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggthemes)
library(GGally)
library(colorspace)
library(readxl)
library(forcats)
library(scales)
library(gridExtra)
library(cowplot)
library(Cairo)

col_grey = "#999999"
col_trueneg = "#56B4E9" #"#0072B2"
col_truepos = "#D55E00"
black = "#000000"
pal <- c(col_grey, col_trueneg, col_truepos, black)

# Function to set the significant digits to 2
sigfig <- function(vec, n=4){ 
  formatC(signif(vec,digits=n), digits=n,format="fg", flag="#") 
}  

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}


d.orig <- read_excel(here_r("OR_R2Positive.xlsx"))
head(d.orig)
str(d.orig)

# Not imputed model data
not.imp <- d.orig[-1,1:5]
# Imputed model data
imp <- d.orig[-1, c(1:2, 6:8)]

# Names of data
names(not.imp)
names(not.imp) <- c("Variable", "Value", "estimate", "low", "high")
names(imp) <- names(not.imp)

# Create a single data frame of the imputed and not imputed data
d <- bind_rows(mutate(not.imp, Model="Not Imputed"),
               mutate(imp, Model="Imputed"))
head(d)
str(d)
table(d$Model)

# Remove NAs
d <- d %>%
  filter(!is.na(Variable)) %>%
  mutate(estimate = as.numeric(estimate),
         low = as.numeric(low),
         high = as.numeric(high))


# Select the variables to show in the final plot
d$Variable
id.vars <- c(1:3, 12:28, 31:54, 57:58, 74:76)
name.vars <- d$Variable[id.vars]
name.vars

# OR plot
d %>%
  filter(Variable %in% name.vars) %>% 
  filter(!is.na(estimate)) %>%
  mutate(Parameter = paste(Variable, Value, sep = " : ")) %>%
  mutate(Parameter = fct_inorder(factor(Parameter)),
         Parameter = fct_rev(Parameter)) %>%
  ggplot(aes(y=estimate, x=Parameter, color=Model)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_linerange(aes(ymin=low, ymax=high), 
                 position=position_dodge(width=0.5)) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values=pal[c(3,2)]) +
  geom_hline(yintercept = 1, lty=3, color="grey30") +
  scale_y_continuous(trans = "log10", 
                     breaks=c(0.00001, 0.01, 1, 100),
                     labels=c("0.00001", "0.01", "1", "100")) +
   # scale_y_continuous(trans = log_trans(), breaks = base_breaks(),
   #                    labels = prettyNum) + 
  ylab("Estimate (95% CI)") + xlab("") +
  theme_bw()

# Save plot as pdf and png
ggsave(here_output(
  file="Figure-5.pdf"), 
       device=cairo_pdf,
       width=12, height=8)
ggsave(here_output(
  file="Figure-5.png"), 
     width=12, height=8)

