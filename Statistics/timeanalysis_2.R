rm(list=ls(all=TRUE))

# specify location to read files
setwd("C:/COVID Sentinel Kohorte Munich/KoCo19 Data/KOCO19V2")
library(ggplot2)
library(rmarkdown)
library(readxl)
library(magrittr)
library(summarytools)
library(formattable)
library(ggplot2)
library(zoo)
library(FSA)
library(reshape2)
library(plyr)
library(dplyr)
library(lubridate)
library(ggpubr)
library(arsenal)
library(stringi)
library(stringr)
library(tidyr)
library(foreign)
library(stats4)
library(glmmTMB)
library(bbmle)
library(countreg)
library(cowplot) 
library(lme4) 
library(sjPlot) 
library(sjmisc) 
library(effects)
library(sjstats) 
library(tidyquant)
library(binom)
library(geepack)
library(gee)
library(BAS)
library(BMS)
library(DTComPair)
library(ggpmisc)
library(pscl)
library(hurdlr)
library(PupillometryR)
library(broom.mixed)
library(emmeans)
library(epitools)
library(patchwork)
library(svglite)
library(visdat)
library(JointAI)
library(pbkrtest)
library(openxlsx)
library(brms)
library(bayestestR)
library(parameters)
library(vcd)
library(tidybayes)

koco_r2_complete<-read.csv("koco_r2_matched.csv",header = TRUE)

# Exlcude intermediate observations
koco_r2_complete<- koco_r2_complete[which(koco_r2_complete$DBS_Result!="Intermediate"),]

time_dat<- subset(koco_r2_complete, select = c("ind_id","DBS_Result","Date_DBS_prick"))
colnames(time_dat)<-c("ind_id","R2_Result","date_R2")
koco_weights<- read.csv("KoCo_weights.csv",header = TRUE)

time_dat<- merge.data.frame(time_dat,koco_weights[c("ind_id","w_ind_samp")])

time_dat$week <- round(difftime(time1 = time_dat$date_R2, time2 = as.Date("2020-11-01",format="%Y-%m-%d"),
                                units = "weeks"),0) 

time_dat$week_cat<- ifelse(time_dat$week<=1,"Week1",
                           ifelse(time_dat$week==2,"Week2",
                                  ifelse(time_dat$week==3,"Week3",
                                         ifelse(time_dat$week==4,"Week4","Week5+"))))
time_dat$week_cat<- factor(time_dat$week_cat,
                           levels = c("Week1","Week2","Week3","Week4","Week5+"))

# Including Sampling Weights

d1<- reshape2::dcast(time_dat[which(time_dat$R2_Result=="Positive"),], week_cat~., value.var = "w_ind_samp",
                     fun.aggregate = sum)
d2<- reshape2::dcast(time_dat, week_cat~., value.var = "w_ind_samp",
                     fun.aggregate = sum)

colnames(d1)[2]<-"Pos_weights"
colnames(d2)[2]<-"Total_weights"

d<- merge.data.frame(d1,d2,by="week_cat")

d$Pos_weights<- round(d$Pos_weights,0)
d$Total_weights<- round(d$Total_weights,0)


fit1<-
  brm(data = d, family = binomial,
      Pos_weights | trials(Total_weights) ~  week_cat,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 10), class = b)),
      chains = 3,warmup = 5000, 
      iter = 10000,thin = 10,
      control = list(adapt_delta = 0.95),
      seed = 1234)
mcmc_plot(fit1, 
          type = "trace")

summary(fit1)

f<-cbind(d,fitted(fit1))
f$prev_est<- f$Estimate/f$Total_weights
f$prev_2.5<- f$Q2.5/f$Total_weights
f$prev_97.5.5<- f$Q97.5/f$Total_weights

prev_tab_weighted_bin<- f

d$log_Total_weights<- log(d$Total_weights)
fit2 <- brm(data = d, family = poisson,
            Pos_weights ~ offset(log_Total_weights) + week_cat,
            prior = c(prior(normal(0, 10), class = Intercept),
                      prior(normal(0, 10), class = b)),
            chains = 3,warmup = 5000, 
            iter = 10000,thin = 10,
            control = list(adapt_delta = 0.95),
            seed = 1234)

summary(fit2)


conditional_effects(fit2, conditions = data.frame(size = 1))
conditional_effects(fit1, conditions = data.frame(size = 1))


d<-cbind(d,fitted(fit2))
d$prev_est<- d$Estimate/d$Total_weights
d$prev_2.5<- d$Q2.5/d$Total_weights
d$prev_97.5.5<- d$Q97.5/d$Total_weights

prev_tab_weighted_pois<- d


# without including sampling weights

d<- reshape2::dcast(time_dat,week_cat~R2_Result)
d$Sum<- d$Positive+d$Negative
d$logSum<- log(d$Sum)

fit1 <- brm(
  Positive|trials(Sum)  ~ week_cat, data =d, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit1, 
          type = "trace")


d1<-cbind(d,fitted(fit1))
d1$prev_est<- d1$Estimate/d1$Sum
d1$prev_2.5<- d1$Q2.5/d1$Sum
d1$prev_97.5.5<- d1$Q97.5/d1$Sum

prev_tab_crude_bin<-d1

conditional_effects(fit1, conditions = data.frame(size = 1))
d$log_Sum<- log(d$Sum)
fit2 <- brm(data = d, family = poisson,
            Positive ~ offset(log_Sum) + week_cat,
            prior = c(prior(normal(0, 10), class = Intercept),
                      prior(normal(0, 10), class = b)),
            chains = 3,warmup = 5000, 
            iter = 10000,thin = 10,
            control = list(adapt_delta = 0.95),
            seed = 1234)


f<-cbind(d,fitted(fit2))
f$prev_est<- f$Estimate /f$Sum
f$prev_2.5<- f$Q2.5 /f$Sum
f$prev_97.5.5<- f$Q97.5/f$Sum

prev_tab_crude_pois<- f


# Table Comparison

prev_tab_crude_pois<-prev_tab_crude_pois[c(1,11,12,13)]
prev_tab_crude_bin<-prev_tab_crude_bin[c(1,10,11,12)]
prev_tab_weighted_pois<-prev_tab_weighted_pois[c(1,9,10,11)]
prev_tab_weighted_bin<-prev_tab_weighted_bin[c(1,8,9,10)]


prev_tab_crude_pois$cat<-"Poisson\nCrude"
prev_tab_crude_bin$cat<- "Binomial\nCrude"
prev_tab_weighted_pois$cat<-"Poisson\nWeighted"
prev_tab_weighted_bin$cat<-"Binomial\nWeighted"

prev<- rbind.fill(prev_tab_crude_pois, prev_tab_crude_bin,
                  prev_tab_weighted_pois,prev_tab_weighted_bin)


ggplot(prev, aes(x=week_cat,y=prev_est))+
  geom_point(aes(colour=cat),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=cat), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Prevalence")+xlab("Recruitment (Excluding DBS Intermediates) Round 2")+
  labs(colour="Method")+theme(legend.position="bottom")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

prev$outcome<- "DBS"

write.csv(prev,"prev_DBS_week.csv")