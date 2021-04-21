#---
#title: "Prevalence Plot"
#author: "Mercè Garí"
#date: '2021-02-10'
#---

here_r = function (...) here::here("Statistics", ...)
here_weights = function (...) here::here("SamplingWeights", ...)
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
sigfig <- function(vec, n=3){ 
  formatC(signif(vec,digits=n), digits=n,format="fg", flag="#") 
}  


d <- read.csv(here_weights("Estimates_Cal_R1_Test.csv"))
head(d)
str(d)

plot <- d %>%
  filter(Calculation %in% "Weighted") %>%
  filter(Categories %in% c("Prevalence R2", "Negative R1 Positive R2")) %>%
  mutate(Categories = case_when(
    Categories == "Prevalence R2" ~ "Prevalence",
    Categories == "Negative R1 Positive R2" ~ "Incidence")) %>%
  # mutate(Estimates = ifelse(Estimates <= 0, 0.9, Estimates),
  #        Lower_95_CB = ifelse(Lower_95_CB <= 0, 0.9, Lower_95_CB)) %>%
  ggplot(aes(y=Estimates, 
             x=Categories, 
             color=Adjustment)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_linerange(aes(ymin=Lower_95_CB,
                     ymax=Upper_95_CB), 
                 position=position_dodge(width=0.5)) +
  coord_flip() +
  scale_linetype_manual(values=c(1,6)) +
  scale_color_manual(values=pal[c(3,2,1,4)]) +
  xlab("") + ylab("Percentage") +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.justification=c(1,0), legend.position=c(1,0.05)) +
  geom_text(data=.%>% filter(Adjustment == "Adjusted"),
            aes(label=paste0(signif(Estimates, 2), " (", 
                             signif(Lower_95_CB, 2), " - ", 
                             signif(Upper_95_CB, 2), ")"), 
                hjust=1.8, vjust=2),
            show.legend=FALSE, size=3) +
  geom_text(data=.%>% filter(Adjustment == "Unadjusted"),
            aes(label=paste0(signif(Estimates, 2), " (", 
                             signif(Lower_95_CB, 2), " - ", 
                             signif(Upper_95_CB, 2), ")"), 
                hjust=1.8, vjust=-1.1),
            show.legend=FALSE, size=3) +
  ylim(c(0,4.5))
plot


ggsave(here_output(
  file="Figure_4.pdf"), 
       device=cairo_pdf,
       width=5, height=3)
ggsave(here_output(
  file="Figure_4.png"), 
     width=5, height=3)

