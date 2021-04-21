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
library(sp)
library(tmap)
library(ape)
library(mgcv)
library(geojsonR)
library(geojsonio)
library(geojsonsf)
library(mapproj)
library(spdep)
load("C:/COVID Sentinel Kohorte Munich/KoCo19 Data/Random Routes/Data-for-prevalence-maps.RData")
ls()
# Reading the shape file
spdf <- geojson_read("http://data.insideairbnb.com/germany/bv/munich/2020-06-20/visualisations/neighbourhoods.geojson",  what = "sp")
# Viewing the shape file
ggplot() +
  geom_polygon(data = spdf, aes( x = long, y = lat, group = group), fill="#69b3a2", color="white") +
  theme_void() +
  coord_map()
bei.adj.q <- poly2nb(spdf,row.names = rownames(spdf@data))
W.bin <- nb2listw(bei.adj.q, style = "B")
W.rs <- nb2listw(bei.adj.q, style = "W")
W.bin
W.rs
W.mat <- nb2mat(bei.adj.q, style = "B")
W.mat_cont <- nb2mat(bei.adj.q, style = "W")


koco_r2_complete<-read.csv("koco_r2_matched.csv",header = TRUE)

sp_dat<- subset(koco_r2_complete, select = c(ind_id,hh_id,R2_Result))

koco_weights<- read.csv("KoCo_weights.csv",header = TRUE)
dist_dat<- read.csv("Munich_Districts.csv",header = TRUE)

pop_dens<- read.xlsx("Muc_popn_density_Wiki.xlsx",colNames = TRUE)
pop_dens$density<- as.numeric(pop_dens$density)

c_dens<- read.xlsx("Constituencywisedensity.xlsx",colNames = TRUE)

rr_dat<-read.csv("KoCo19_Haushalte4Modeler_wRRstartConstituency_20200910.csv",header = TRUE,
                 na.strings = c(""," "))
colnames(rr_dat)[1]<-"hh_id"

rr_dat$city_district<- ifelse(is.na(rr_dat$city_district),rr_dat$suburb,rr_dat$city_district)

sp_dat<- merge.data.frame(sp_dat,koco_weights[c("ind_id","hh_id","w_ind_samp")], id="ind_id")

hh_dat<- reshape2::dcast(sp_dat, hh_id~R2_Result, value.var = 'w_ind_samp', fun.aggregate = sum)
hh_dat$sum<- hh_dat$Negative+hh_dat$Positive

hh_dat<- merge.data.frame(hh_dat,rr_dat[c("hh_id","city_district")], by="hh_id", all.x = TRUE)

hh_dat$dist<- as.numeric(substr(hh_dat$city_district,13,14))

hh_dat<- merge.data.frame(hh_dat,dist_dat, by="dist",all.x = TRUE)

hh_dat<- merge.data.frame(hh_dat, rr_dat[c("hh_id","constituency")],by="hh_id")
hh_dat$c1<- str_sub(hh_dat$constituency, start= -2)
c_dens$c1<-str_sub(c_dens$constituency, start= -2)
c_dens$dist<-ifelse(c_dens$constituency<10000,str_sub(c_dens$constituency, start= 1,end=1),
                  str_sub(c_dens$constituency, start= 1,end=2))
hh_dat[which(hh_dat$hh_id %in% c("2374T00","3706K00","2664P00")),]$dist<-4
hh_dat<- merge.data.frame(hh_dat,c_dens[c("c1","dist","density_Mar2020","density_Dec2020")], by=c("c1","dist"), all.x = TRUE)


# Aggregate by district

t1<- reshape2::dcast(hh_dat, c1+dist+ density_Mar2020~., value.var = 'Positive', fun.aggregate = sum)
t2<- reshape2::dcast(hh_dat, c1+dist+ density_Mar2020~., value.var = 'sum', fun.aggregate = sum)

tab_dist<- merge.data.frame(t1,t2, by=c("c1","dist" ,"density_Mar2020"))
colnames(tab_dist)[c(4,5)]<-c("Positives","Total")
tab_dist$crude_prev<- tab_dist$Positives/tab_dist$Total

ggplot(tab_dist,aes(x=density_Mar2020,y=crude_prev))+ geom_point()+geom_smooth()

tab_dist$Positives<- round(tab_dist$Positives,0)
tab_dist$Total<- round(tab_dist$Total,0)

formula<- "Positives ~ offset(log(Total)) + density_Mar2020"
f1<- glm(formula = formula, data = tab_dist,family = "quasipoisson")
f2<- glm(formula = formula, data = tab_dist,family = "poisson")
f3<- glm(cbind(Positives ,Total-Positives )~ density_Mar2020,data = tab_dist, family = "binomial")
f4 <- glm.nb(formula = formula, data = tab_dist)

summary(f1)
summary(f2)
summary(f3)
summary(f4)

tab_dist$pop_cat<- ifelse(tab_dist$density_Mar2020<50,"<50",
                          ifelse(tab_dist$density_Mar2020>= 50 & tab_dist$density_Mar2020<100,"50-100",
                                 ifelse(tab_dist$density_Mar2020>=100 & tab_dist$density_Mar2020<150,"100-150","150+")))

tab_dist$logSum<- log(tab_dist$Total)
tab_dist$pop_cat<- factor(tab_dist$pop_cat, levels = c("<50","50-100","100-150","150+"))

fit_zip_w <- brm(data = tab_dist, 
                   bf(Positives ~ offset(logSum) + density_Mar2020 + (1|dist),
                      zi~density_Mar2020),
                   prior = c(prior(normal(0, 10), class = Intercept),
                             prior(normal(0, 10), class = b)),
                   chains = 3,warmup = 5000, 
                   iter = 10000,thin = 1,
                   control = list(adapt_delta = 0.95),
                   seed = 1234, family = zero_inflated_poisson())

mcmc_plot(fit_zip_w, 
          type = "trace")

summary(fit_zip_w)
fit_zip_w$

M4 <- zeroinfl(Positives ~ offset(logSum) + density_Mar2020|
                 density_Mar2020,
               dist = 'negbin',
               data =tab_dist)
summary(M4)

confint(M4)

M5 <- zeroinfl(Positives ~ offset(logSum) + density_Mar2020|
                 density_Mar2020,
               dist = 'poisson',
               data =tab_dist)
summary(M5)

confint(M5)

t1<- reshape2::dcast(tab_dist,pop_cat~., value.var = "Positives", fun.aggregate = sum)
t2<- reshape2::dcast(tab_dist,pop_cat~., value.var = "Total", fun.aggregate = sum)

t<- merge.data.frame(t1,t2,by="pop_cat")
colnames(t)<-c("dens_cat","Positives","Total")
t$Positives<- round(t$Positives, 0)
t$Total<- round(t$Total,0)
t$dens_cat<- factor(t$dens_cat, levels = c("<50","50-100","100-150","150+"))
fit1<-
  brm(data = t, family = binomial,
      Positives | trials(Total) ~  dens_cat,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 10), class = b)),
      chains = 3,warmup = 5000, 
      iter = 10000,thin = 1,
      control = list(adapt_delta = 0.95),
      seed = 1234)
mcmc_plot(fit1, 
          type = "trace")

summary(fit1)

f<-cbind(t,fitted(fit1))
f$prev_est<- f$Estimate/f$Total
f$prev_2.5<- f$Q2.5/f$Total
f$prev_97.5.5<- f$Q97.5/f$Total

prev_tab_weighted_bin<- f

t$log_Total_weights<- log(t$Total)
fit2 <- brm(data = t, family = poisson,
            Positives ~ offset(log_Total_weights) + dens_cat,
            prior = c(prior(normal(0, 10), class = Intercept),
                      prior(normal(0, 10), class = b)),
            chains = 3,warmup = 5000, 
            iter = 10000,thin = 1,
            control = list(adapt_delta = 0.95),
            seed = 1234)

summary(fit2)


conditional_effects(fit2, conditions = data.frame(size = 1))
conditional_effects(fit1, conditions = data.frame(size = 1))


d<-cbind(t,fitted(fit2))
d$prev_est<- d$Estimate/d$Total
d$prev_2.5<- d$Q2.5/d$Total
d$prev_97.5.5<- d$Q97.5/d$Total

prev_tab_weighted_pois<- d


# Aggregate by district Crude
hh_1<- reshape2::dcast(sp_dat, hh_id~R2_Result)
hh_1$sum<- hh_1$Negative+hh_1$Positive
colnames(hh_1)[c(2,3,4)]<-c("Neg_c","Pos_c","Sum_c")
hh_1<- merge.data.frame(hh_1, hh_dat, by="hh_id")

tab_dist1<- reshape2::dcast(hh_1, c1+dist+ density_Mar2020~., value.var = 'Pos_c', fun.aggregate = sum)
tab_dist2<- reshape2::dcast(hh_1, c1+dist+ density_Mar2020~., value.var = 'Sum_c', fun.aggregate = sum)

tab_dist1<- merge.data.frame(tab_dist1,tab_dist2, by=c("c1","dist" ,"density_Mar2020"))
colnames(tab_dist1)[c(4,5)]<-c("Positives","Total")
tab_dist1$crude_prev<- tab_dist1$Positives/tab_dist1$Total

ggplot(tab_dist1,aes(x=density_Mar2020,y=crude_prev))+ geom_point()+geom_smooth()

formula<- "Positives ~ offset(log(Total)) + density_Mar2020"
f1<- glm(formula = formula, data = tab_dist1,family = "quasipoisson")
f2<- glm(formula = formula, data = tab_dist1,family = "poisson")
f3<- glm(cbind(Positives ,Total-Positives )~ density_Mar2020,data = tab_dist1, family = "binomial")
f4 <- glm.nb(formula = formula, data = tab_dist1)

summary(f1)
summary(f2)
summary(f3)
summary(f4)

tab_dist1$pop_cat<- ifelse(tab_dist1$density_Mar2020<50,"<50",
                          ifelse(tab_dist1$density_Mar2020>= 50 & tab_dist1$density_Mar2020<100,"50-100",
                                 ifelse(tab_dist1$density_Mar2020>=100 & tab_dist1$density_Mar2020<150,"100-150","150+")))

tab_dist1$logSum<- log(tab_dist1$Total)
tab_dist1$pop_cat<- factor(tab_dist1$pop_cat, levels = c("<50","50-100","100-150","150+"))


fit_zip <- brm(brm(data = tab_dist1, 
                     bf(Positives ~ offset(logSum) + pop_cat + (1|dist),
                        zi~pop_cat),
                     prior = c(prior(normal(0, 10), class = Intercept),
                               prior(normal(0, 10), class = b)),
                     chains = 3,warmup = 5000, 
                     iter = 10000,thin = 10,
                     control = list(adapt_delta = 0.95),
                     seed = 1234), family = zero_inflated_poisson())

mcmc_plot(fit_zip, 
          type = "trace")

summary(fit_zip)


t1<- reshape2::dcast(tab_dist1,pop_cat~., value.var = "Positives", fun.aggregate = sum)
t2<- reshape2::dcast(tab_dist1,pop_cat~., value.var = "Total", fun.aggregate = sum)

t<- merge.data.frame(t1,t2,by="pop_cat")
colnames(t)<-c("dens_cat","Positives","Total")
t$dens_cat<- factor(t$dens_cat, levels = c("<50","50-100","100-150","150+"))
fit1<-
  brm(data = t, family = binomial,
      Positives | trials(Total) ~  dens_cat,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 10), class = b)),
      chains = 3,warmup = 5000, 
      iter = 10000,thin = 1,
      control = list(adapt_delta = 0.95),
      seed = 1234)
mcmc_plot(fit1, 
          type = "trace")

summary(fit1)

f<-cbind(t,fitted(fit1))
f$prev_est<- f$Estimate/f$Total
f$prev_2.5<- f$Q2.5/f$Total
f$prev_97.5.5<- f$Q97.5/f$Total

prev_tab_crude_bin<- f

t$log_Total_weights<- log(t$Total)
fit2 <- brm(data = t, family = poisson,
            Positives ~ offset(log_Total_weights) + dens_cat,
            prior = c(prior(normal(0, 10), class = Intercept),
                      prior(normal(0, 10), class = b)),
            chains = 3,warmup = 5000, 
            iter = 10000,thin = 1,
            control = list(adapt_delta = 0.95),
            seed = 1234)

summary(fit2)


conditional_effects(fit2, conditions = data.frame(size = 1))
conditional_effects(fit1, conditions = data.frame(size = 1))


d<-cbind(t,fitted(fit2))
d$prev_est<- d$Estimate/d$Total
d$prev_2.5<- d$Q2.5/d$Total
d$prev_97.5.5<- d$Q97.5/d$Total

prev_tab_crude_pois<- d


prev_tab_crude_bin$dist<-"Binomial"
prev_tab_crude_bin$sample<-"Crude"
prev_tab_weighted_bin$dist<-"Binomial"
prev_tab_weighted_bin$sample<-"Sample weighted"
prev_tab_crude_pois$dist<-"Poisson"
prev_tab_crude_pois$sample<-"Crude"
prev_tab_weighted_pois$dist<-"Poisson"
prev_tab_weighted_pois$sample<-"Sample weighted"

prev<- rbind( prev_tab_crude_bin[c(1,8,9,10,11,12)],
              prev_tab_crude_pois[c(1,9,10,11,12,13)],
              prev_tab_weighted_bin[c(1,8,9,10,11,12)],
              prev_tab_weighted_pois[c(1,9,10,11,12,13)])

ggplot(prev, aes(y=prev_est,x=dens_cat))+
  geom_point(aes(colour=dist),position = position_dodge(width = 0.9))+ 
  geom_errorbar(aes(ymin =prev_2.5, ymax =prev_97.5.5, colour=dist), position = position_dodge(width = 0.9),
                width = 0.2)+
  facet_grid(.~sample,scales="free")+
  ylab("Prevalence")+ xlab("Population density (persons/100 sq km)")+theme_bw()+
  labs(colour="Outcome distribution")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+theme(legend.position="bottom")

ggsave("prevalence vs density_discrete.jpeg",dpi=720, units="in", height=5, width = 7)
