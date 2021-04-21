#---
#title: "Plot Registered cases Munich"
#author: "Mercè Garí"
#date: '2021-03-16'
#---


here_r = function (...) here::here("Statistics", ...)
here_output = function (...) here::here("Output", ...)

# Load packages
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
library(stringr)


col_grey = "#999999"
col_trueneg = "#56B4E9" #"#0072B2"
col_truepos = "#D55E00"
black = "#000000"
pal <- c(col_grey, col_trueneg, col_truepos, black)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# Load the data from the LMU nowcast
# https://corona.stat.uni-muenchen.de/nowcast/

d.orig <- read.csv(here_r("results_nowcast_2021-03-16.csv"))

head(d.orig)  
dim(d.orig)

ggplot(d.orig, aes(x=date, y=reported)) +
         geom_bar(stat="identity")

plot.framework <- d.orig %>%
  mutate(month = str_sub(date, start = 6, end=7),
         year = str_sub(date, start=1, end=4)) %>%
  group_by(month, year) %>%
  summarize(n = sum(reported)) %>%
  mutate(date = paste0(year, "-", month)) %>%
  filter(!date %in% c("2020-02", "2021-03", "2021-02")) %>%
  ggplot(aes(x=date, y=n)) +
  geom_bar(stat="identity", width=0.75) +
  theme_classic() +
   xlab("") + ylab("Num. Cases in Munich") +
   geom_rect(aes(xmin = 1.5, xmax = 4.5, ymin = -Inf, ymax = Inf),
                    fill = "#E69F00", alpha = 0.02) +
   geom_rect(aes(xmin = 4.5, xmax = 8.5, ymin = -Inf, ymax = Inf),
                    fill = "#56B4E9", alpha = 0.02) +
   geom_rect(aes(xmin = 8.5, xmax = 11.5, ymin = -Inf, ymax = Inf),
                    fill = "#CC79A7", alpha = 0.02) +
  annotate(geom="text", x=3, y=16000, label="Baseline study", size=4) +
   annotate(geom="text", x=3, y=14750, 
            label="(Questionnaires,\nSARS-CoV-2 serology)", size=3) +
  annotate(geom="text", x=6.5, y=16000, 
           label="Questionnaire follow-up", size=4) +
    annotate(geom="text", x=6.5, y=14750, 
           label="(Leisure time activities,\nrisk perception, social contacts)", size=3) +
  annotate(geom="text", x=10, y=16000, 
           label="1st follow-up", size=4) +
   annotate(geom="text", x=10, y=15000, 
           label="(SARS-CoV-2 serology)", size=3) +
  ylim(c(0, 16000))
  
plot.framework
ggsave(here_output(
  file="Figure_1.pdf"), 
       device=cairo_pdf,
       width=8, height=4)

ggsave(here_output(
  file="Figure_1.png"), 
       width=8, height=4)

