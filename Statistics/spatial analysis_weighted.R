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

pop_dens<- read.xlsx("Muc_popn_density.xlsx",colNames = TRUE)
pop_dens$density<- as.numeric(pop_dens$density)

rr_dat<-read.csv("KoCo19_Haushalte4Modeler_wRRstartConstituency_20200910.csv",header = TRUE,
                 na.strings = c(""," "))
colnames(rr_dat)[1]<-"hh_id"

rr_dat$city_district<- ifelse(is.na(rr_dat$city_district),rr_dat$suburb,rr_dat$city_district)

sp_dat<- merge.data.frame(sp_dat,koco_weights[c("hh_id","w_ind_samp")], by="hh_id")

hh_dat<- reshape2::dcast(sp_dat, hh_id~R2_Result, value.var = 'w_ind_samp', fun.aggregate = sum)
hh_dat$sum<- hh_dat$Negative+hh_dat$Positive

hh_dat<- merge.data.frame(hh_dat,rr_dat[c("hh_id","city_district")], by="hh_id", all.x = TRUE)

hh_dat$dist<- as.numeric(substr(hh_dat$city_district,13,14))

hh_dat<- merge.data.frame(hh_dat,dist_dat, by="dist",all.x = TRUE)


# Prevalence by districts
tab_dist<- reshape2::dcast(hh_dat, dist+name+ lat+long~., value.var = "sum",fun.aggregate = sum)
colnames(tab_dist)[5]<-"Tested_Ind"

t1<-reshape2::dcast(hh_dat, dist+name~., value.var = "Positive",fun.aggregate = sum)
colnames(t1)[3]<- "Positive_Ind"
tab_dist<- merge.data.frame(tab_dist,t1,by=c("dist","name"))
tab_dist<- merge.data.frame(tab_dist,pop_dens[c(2,3)],by="dist")

tab_dist$prev_crude<- tab_dist$Positive_Ind/tab_dist$Tested_Ind

tab_dist$Positive_Ind<- round(tab_dist$Positive_Ind,0)
tab_dist$Tested_Ind<- round(tab_dist$Tested_Ind,0)

moran.mc(x=tab_dist$prev_crude, listw=W.rs, nsim=10000)
moran.mc(x=tab_dist$prev_crude, listw=W.bin, nsim=10000)
moran.mc(x=tab_dist$density, listw=W.bin, nsim=10000)
moran.mc(x=tab_dist$density, listw=W.bin, nsim=10000)

cor.test(tab_dist$density,tab_dist$prev_crude,method="spearman", exact = FALSE)
cor.test(tab_dist$density,tab_dist$prev_crude,method="kendall",exact = FALSE)
cor.test(tab_dist$density,tab_dist$prev_crude,method="pearson")

colnames(tab_dist)[2]<-"neighbourhood"
tab_dist$neighbourhood[15]<-"Tudering-Riem" 
tab_dist$neighbourhood[17]<-"Obergiesing" 

prevdata.sp <- merge(x=spdf, y=tab_dist, by="neighbourhood")

library(rgdal)
library(CARBayes)

formula <- Positive_Ind ~ offset(log(Tested_Ind)) + density

f1<- glm(formula = formula, data = prevdata.sp@data,family = "quasipoisson")
f2<- glm(formula = formula, data = prevdata.sp@data,family = "poisson")
f3<- glm(cbind(Positive_Ind ,Tested_Ind-Positive_Ind )~ density,data = prevdata.sp@data, family = "binomial")


summary(f1)
summary(f2)
summary(f3)


moran.mc(x=residuals(f1), listw=W.bin, nsim=10000)
moran.mc(x=residuals(f1), listw=W.rs, nsim=10000)
moran.mc(x=residuals(f2), listw=W.bin, nsim=10000)
moran.mc(x=residuals(f2), listw=W.rs, nsim=10000)
moran.mc(x=residuals(f3), listw=W.bin, nsim=10000)
moran.mc(x=residuals(f3), listw=W.rs, nsim=10000)
set.seed(1234)
s1 <- S.CARleroux(
  Positive_Ind ~ offset(log(Tested_Ind)) + density,
  data = prevdata.sp@data,
  family = "poisson",
  W = W.mat,
  rho = 1,        # Fit intrinsic CAR model
  burnin = 10000,
  n.sample = 25000
)

set.seed(1234)
s1.0 <- S.CARleroux(
  Positive_Ind ~ offset(log(Tested_Ind)) ,
  data = prevdata.sp@data,
  family = "poisson",
  W = W.mat,
  rho = 1,        # Fit intrinsic CAR model
  burnin = 10000,
  n.sample = 25000
)
set.seed(1234)
s2 <- S.CARleroux(
  Positive_Ind ~ density,
  data = prevdata.sp@data ,
  family = "binomial",
  trials = prevdata.sp@data$Tested_Ind,
  W = W.mat,
  rho = 1,        # Fit intrinsic CAR model
  burnin = 10000,
  n.sample = 25000
)
set.seed(1234)
s2.0 <- S.CARleroux(
  Positive_Ind ~ 1,
  data = prevdata.sp@data ,
  family = "binomial",
  trials = prevdata.sp@data$Tested_Ind,
  W = W.mat,
  rho = 1,        # Fit intrinsic CAR model
  burnin = 10000,
  n.sample = 25000
)
print(s1$summary.results)
print(s1.0$summary.results)

print(s2$summary.results)
print(s2.0$summary.results)



y.fit <- s1$samples$fitted
SIR <- t(t(y.fit) / prevdata.sp@data$Tested_Ind)
map<- NULL
SIR.50 <- apply(SIR, 2, median)
SIR.025 <- apply(SIR, 2, quantile, 0.025)
SIR.975 <- apply(SIR, 2, quantile, 0.975)
map_Pois_adj<- cbind.data.frame(s1$fitted.values,prevdata.sp@data$Tested_Ind,SIR.50,SIR.025,SIR.975)
map_Pois_adj$Outcome_dist<-"Poisson"
map_Pois_adj$smoothing<-"Global"
y.fit <- s1.0$samples$fitted
SIR <- t(t(y.fit) / prevdata.sp@data$Tested_Ind)
map<- NULL
SIR.50 <- apply(SIR, 2, median)
SIR.025 <- apply(SIR, 2, quantile, 0.025)
SIR.975 <- apply(SIR, 2, quantile, 0.975)
map_Pois<- cbind.data.frame(s1.0$fitted.values,prevdata.sp@data$Tested_Ind,SIR.50,SIR.025,SIR.975)
map_Pois$Outcome_dist<-"Poisson"
map_Pois$smoothing<-"Global"
y.fit <- s2$samples$fitted
SIR <- t(t(y.fit) / prevdata.sp@data$Tested_Ind)
map<- NULL
SIR.50 <- apply(SIR, 2, median)
SIR.025 <- apply(SIR, 2, quantile, 0.025)
SIR.975 <- apply(SIR, 2, quantile, 0.975)
map_Binom_adj<- cbind.data.frame(s2$fitted.values,prevdata.sp@data$Tested_Ind,SIR.50,SIR.025,SIR.975)
map_Binom_adj$Outcome_dist<-"Binomial"
map_Binom_adj$smoothing<-"Global"
y.fit <- s2.0$samples$fitted
SIR <- t(t(y.fit) / prevdata.sp@data$Tested_Ind)
map<- NULL
SIR.50 <- apply(SIR, 2, median)
SIR.025 <- apply(SIR, 2, quantile, 0.025)
SIR.975 <- apply(SIR, 2, quantile, 0.975)
map_Binom<- cbind.data.frame(s2.0$fitted.values,prevdata.sp@data$Tested_Ind,SIR.50,SIR.025,SIR.975)
map_Binom$Outcome_dist<-"Binomial"
map_Binom$smoothing<-"Global"

# CAR Dissimilarity
Z.dens <- as.matrix(dist(cbind(prevdata.sp$density,
                               prevdata.sp$density), method = "manhattan", diag = TRUE,
                         upper = TRUE)) * W.mat/2

s1_dis <-  S.CARdissimilarity(Positive_Ind ~ offset(log(Tested_Ind))+ density,
                              data = prevdata.sp@data, W = W.mat, Z = list(Z.dens = Z.dens), 
                              family = "poisson",
                              burnin = 10000, n.sample = 25000)
print(s1_dis$summary.results)
y.fit <- s1_dis$samples$fitted
SIR <- t(t(y.fit) / prevdata.sp@data$Tested_Ind)
map<- NULL
SIR.50 <- apply(SIR, 2, median)
SIR.025 <- apply(SIR, 2, quantile, 0.025)
SIR.975 <- apply(SIR, 2, quantile, 0.975)
map_Pois_dis_adj<-cbind.data.frame(s1_dis$fitted.values,prevdata.sp@data$Tested_Ind,SIR.50,SIR.025,SIR.975)
map_Pois_dis_adj$Outcome_dist<-"Poisson"
map_Pois_dis_adj$smoothing<-"Local"

s1.0_dis <-  S.CARdissimilarity(Positive_Ind ~ offset(log(Tested_Ind)),
                                data = prevdata.sp@data, W = W.mat, Z = list(Z.dens = Z.dens), 
                                family = "poisson",
                                burnin = 10000, n.sample = 25000)
print(s1.0_dis$summary.results)
y.fit <- s1.0_dis$samples$fitted
SIR <- t(t(y.fit) / prevdata.sp@data$Tested_Ind)
map<- NULL
SIR.50 <- apply(SIR, 2, median)
SIR.025 <- apply(SIR, 2, quantile, 0.025)
SIR.975 <- apply(SIR, 2, quantile, 0.975)
map_Pois_dis<-cbind.data.frame(s1.0_dis$fitted.values,prevdata.sp@data$Tested_Ind,SIR.50,SIR.025,SIR.975)
map_Pois_dis$Outcome_dist<-"Poisson"
map_Pois_dis$smoothing<-"Local"

s2_dis <-  S.CARdissimilarity(Positive_Ind ~  density,
                              data = prevdata.sp@data, W = W.mat, Z = list(Z.dens = Z.dens), 
                              family = "binomial",trials = prevdata.sp@data$Tested_Ind,
                              burnin = 10000, n.sample = 25000)
print(s2_dis$summary.results)
y.fit <- s1_dis$samples$fitted
SIR <- t(t(y.fit) / prevdata.sp@data$Tested_Ind)
map<- NULL
SIR.50 <- apply(SIR, 2, median)
SIR.025 <- apply(SIR, 2, quantile, 0.025)
SIR.975 <- apply(SIR, 2, quantile, 0.975)
map_Binom_dis_adj<-cbind.data.frame(s2_dis$fitted.values,prevdata.sp@data$Tested_Ind,SIR.50,SIR.025,SIR.975)
map_Binom_dis_adj$Outcome_dist<-"Binomial"
map_Binom_dis_adj$smoothing<-"Local"

s2.0_dis <-  S.CARdissimilarity(Positive_Ind ~ 1,
                                data = prevdata.sp@data, W = W.mat, Z = list(Z.dens = Z.dens), 
                                family = "binomial",trials = prevdata.sp@data$Tested_Ind,
                                burnin = 10000, n.sample = 25000)
print(s2.0_dis$summary.results)
y.fit <- s2.0_dis$samples$fitted
SIR <- t(t(y.fit) / prevdata.sp@data$Tested_Ind)
map<- NULL
SIR.50 <- apply(SIR, 2, median)
SIR.025 <- apply(SIR, 2, quantile, 0.025)
SIR.975 <- apply(SIR, 2, quantile, 0.975)
map_Binom_dis<-cbind.data.frame(s2.0_dis$fitted.values,prevdata.sp@data$Tested_Ind,SIR.50,SIR.025,SIR.975)
map_Binom_dis$Outcome_dist<-"Binomial"
map_Binom_dis$smoothing<-"Local"

colnames(map_Binom)[c(1,2)]<-colnames(map_Binom_adj)[c(1,2)]<- colnames(map_Binom_dis)[c(1,2)]<-
  colnames(map_Binom_dis_adj)[c(1,2)]<-c("Est_Count","Tested_Ind")
colnames(map_Pois)[c(1,2)]<-colnames(map_Pois_adj)[c(1,2)]<- colnames(map_Pois_dis)[c(1,2)]<-
  colnames(map_Pois_dis_adj)[c(1,2)]<-c("Est_Count","Tested_Ind")
map_weighted<- rbind(map_Binom,map_Binom_adj,map_Binom_dis,map_Binom_dis_adj,
                  map_Pois,map_Pois_adj,map_Pois_dis,map_Pois_dis_adj)
map_weighted$prev_type<-rep(c("Not adjusted","Adj. Popn. Density"),each=25, times=4)
map_weighted$outcome_type<- "Weighted"

write.csv(map_weighted,"spatial_prev_weighted.csv")

map_crude<- read.csv("spatial_prev_crude.csv",header = TRUE)
map_crude$X<-NULL
map_crude$outcome_type<-"Crude"
map<- rbind(map_crude,map_weighted)
map$Borough<- rep(tab_dist$neighbourhood,times=16)
write.csv(map,"District_Seroprev.csv")

tab_dist$type<-"Weighted"
tab_un<- read.csv("Simple_seroprev_unwght.csv",header = TRUE)
tab_un$type<-"Unweighted"
tab_un$X<-NULL
tab_dist<- rbind(tab_dist,tab_un)

ggplot(tab_dist, aes(y=prev_crude,x=density))+
  geom_point(aes(colour=type))+ geom_smooth(aes(group=type,colour=type, fill=type))+
  ylab("Crude prevalence (Borough level)")+ xlab("Population per sq-km (Borough level)")+theme_bw()+
  labs(colour="Prevalence Type")+guides(fill=FALSE)+
  scale_colour_discrete(labels=c("Without including sample weights","Including sample weights"))+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+theme(legend.position="bottom")

ggsave("prevalence vs density.jpeg",dpi=720, units="in", height=5, width = 7)

map$smoothing<- factor(map$smoothing, levels = c("Global","Local"),
                       labels = c("Spatial AR - Global","Spatial AR - Local"))
g1<-ggplot(map[which(map$smoothing=="Spatial AR - Global"),], aes(y=SIR.50,x=Borough))+
  geom_point(aes(colour=Outcome_dist),position = position_dodge(width = 0.9))+ 
  geom_errorbar(aes(ymin =SIR.025, ymax =SIR.975, colour=Outcome_dist), position = position_dodge(width = 0.9),
                width = 0.2)+
  facet_grid(prev_type~outcome_type,scales="free")+
  ylab("Prevalence (Borough level)")+ xlab("Borough/District in Munich")+theme_bw()+
  labs(colour="Outcome distribution")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+theme(legend.position="bottom")+
  coord_flip()
g1
ggsave("prevalence vs Borough_spatialAC_Global.jpeg",dpi=720, units="in", height=7, width = 9)

g2<-ggplot(map[which(map$smoothing=="Spatial AR - Local"),], aes(y=SIR.50,x=Borough))+
  geom_point(aes(colour=Outcome_dist),position = position_dodge(width = 0.9))+ 
  geom_errorbar(aes(ymin =SIR.025, ymax =SIR.975, colour=Outcome_dist), position = position_dodge(width = 0.9),
                width = 0.2)+
    facet_grid(prev_type~outcome_type,scales="free")+
  ylab("Prevalence (Borough level)")+ xlab("Borough/District in Munich")+theme_bw()+
  labs(colour="Outcome distribution")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+theme(legend.position="bottom")+
  coord_flip()
g2
ggsave("prevalence vs Borough_spatialAC_Local.jpeg",dpi=720, units="in", height=7, width = 9)
