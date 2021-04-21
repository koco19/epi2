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

bdat_l<-read.csv("bdat_R2.csv",header = TRUE)

bdat<- read.csv("bdat_total.csv",header = TRUE)

koco_r2_dat<-read.csv("koco_r2_dat.csv",header = TRUE)
koco_r2_complete<-read.csv("koco_r2_matched.csv",header = TRUE)
bdat_l<-merge.data.frame(bdat_l,koco_r2_complete[c("ind_id","Roche_Result_new")],by="ind_id")
colnames(bdat_l)[27]<-"R1_Result"
# number of people with free time activity answered
koco_risk<-0
koco_freetime<-0
for(i in 1:nrow(bdat))
{ t1<- bdat[i, c(6,18)]
  s1<- length(t1[which(is.na(t1))])
  t2<- bdat[i, c(20)]
  s2<- length(t2[which(is.na(t2))])
  koco_risk[i]<- ifelse(s1==2,"No entry","Entry")
  koco_freetime[i]<- ifelse(s2==1,"No entry","Entry")
  
}
bdat<- cbind(bdat,koco_freetime,koco_risk)

bdat_freetime<- bdat[which(bdat$koco_freetime=="Entry" &
                             bdat$ind_id %in% bdat_l$ind_id==TRUE),]
bdat_risk<- bdat[which(bdat$koco_risk=="Entry" &
                         bdat$ind_id %in% bdat_l$ind_id==TRUE),]

length(unique(substr(bdat_freetime$ind_id,1,4)))
length(unique(substr(bdat_risk$ind_id,1,4)))

x<-which(bdat_risk$ind_id %in% bdat_freetime$ind_id)
x1<- setdiff(bdat_freetime$ind_id,bdat_risk$ind_id)
# Risk R1
bdat$risk_R1_c<- ifelse(bdat$risk_R1_c=="High","High","Not High")
bdat$risk_R2_c<- ifelse(bdat$risk_R2_c=="High","High","Not High")
d<- reshape2::dcast(bdat, Sex+Agegroup+risk_R1_c~status)
d[is.na(d)]<-"Missing"
d$Sum<- d$`Never Positive`+d$`Atleast once Positive`
d$Pos<- d$`Atleast once Positive`

d$risk_R1_c<- factor(d$risk_R1_c,levels = c("High","Not High","Missing"))

fit1 <- brm(
  Pos|trials(Sum)  ~ Sex*risk_R1_c+Agegroup*risk_R1_c, data =d, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit1, 
          type = "trace")

d<-cbind(d,fitted(fit1))
d$prev_est<- d$Estimate/d$Sum
d$prev_2.5<- d$Q2.5/d$Sum
d$prev_97.5.5<- d$Q97.5/d$Sum

d$Risk<-"Risk1"
write.csv(d,"Risk1_prev.csv")

d$Agegroup<- factor(d$Agegroup, levels = c("<35","35-65",">=65"  ))
d$risk_R1_c<- factor(d$risk_R1_c ,levels = c("Not High","High","Missing"))
ggplot(d, aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=risk_R1_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=risk_R1_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Prevalence")+xlab("Agegroup")+
  labs(colour="Risk perception\nfor self-estimated health-related\nrisk-taking behavior (Summer)")+theme(legend.position="bottom")+
  facet_grid(.~Sex)+scale_y_continuous(labels = scales::percent_format(accuracy = 1))

ggsave("agegroup vs sex vs risk_R1.png")

# Risk R1 for R2
bdat_l$risk_R1_c<- ifelse(bdat_l$risk_R1_c=="High","High","Not High")
bdat_l$risk_R2_c<- ifelse(bdat_l$risk_R2_c=="High","High","Not High")
d<- reshape2::dcast(bdat_l, Sex+Agegroup+risk_R1_c+ R1_Result~R2_Result)
d[is.na(d)]<-"Missing"
d$Sum<- d$Negative+d$Positive
d$Pos<- d$Positive

d$risk_R1_c<- factor(d$risk_R1_c,levels = c("High","Not High","Missing"))

fit1 <- brm(
  Pos|trials(Sum)  ~ Sex*risk_R1_c+Agegroup*risk_R1_c + R1_Result, data =d, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit1, 
          type = "trace")

d<-cbind(d,fitted(fit1))
d$prev_est<- d$Estimate/d$Sum
d$prev_2.5<- d$Q2.5/d$Sum
d$prev_97.5.5<- d$Q97.5/d$Sum

d$Risk<-"Risk1"
write.csv(d,"Risk1_prev_R2.csv")

d$Agegroup<- factor(d$Agegroup, levels = c("<35","35-65",">=65"  ))
d$risk_R1_c<- factor(d$risk_R1_c ,levels = c("Not High","High","Missing"))
f1<-ggplot(d[which(d$R1_Result=="Negative"),], aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=risk_R1_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=risk_R1_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Percentage")+xlab("Agegroup")+
  labs(colour="Risk perception for self-estimated health-\nrelated risk-taking behavior (Summer 2020)")+theme(legend.position="bottom")+
  facet_grid(.~Sex, scales = "free")+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
f1
ggsave("agegroup vs sex vs risk_R1 vs R1res.png")


# Risk R2
d1<- reshape2::dcast(bdat, Sex+Agegroup+risk_R2_c~status)
d1[is.na(d1)]<-"Missing"
d1$Sum<- d1$`Never Positive`+d1$`Atleast once Positive`
d1$Pos<- d1$`Atleast once Positive`


d1$risk_R2_c<- factor(d1$risk_R2_c,levels = c("High","Not High","Missing"))

fit2 <- brm(
  Pos|trials(Sum)  ~ Sex*risk_R2_c+Agegroup*risk_R2_c, data =d1, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit2, 
          type = "trace")


d1<-cbind(d1,fitted(fit2))
d1$prev_est<- d1$Estimate/d1$Sum
d1$prev_2.5<- d1$Q2.5/d1$Sum
d1$prev_97.5.5<- d1$Q97.5/d1$Sum
d1$Risk<-"Risk2"
write.csv(d1,"Risk2_prev.csv")

d1$Agegroup<- factor(d1$Agegroup, levels = c("<35","35-65",">=65"  ))
d1$risk_R2_c<- factor(d1$risk_R2_c ,levels = c("Not High","High","Missing"))
ggplot(d1, aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=risk_R2_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=risk_R2_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Prevalence")+xlab("Agegroup")+
  labs(colour="Risk Perception (Winter)")+theme(legend.position="bottom")+
  facet_grid(.~Sex)+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggsave("agegroup vs sex vs risk_R2.png")

# Risk R2 for R2 result
d1<- reshape2::dcast(bdat_l, Sex+Agegroup+risk_R2_c+R1_Result~R2_Result)
d1[is.na(d1)]<-"Missing"
d1$Sum<- d1$Positive+d1$Negative
d1$Pos<- d1$Positive


d1$risk_R2_c<- factor(d1$risk_R2_c,levels = c("High","Not High","Missing"))

fit2 <- brm(
  Pos|trials(Sum)  ~ Sex*risk_R2_c+Agegroup*risk_R2_c+R1_Result , data =d1, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit2, 
          type = "trace")


d1<-cbind(d1,fitted(fit2))
d1$prev_est<- d1$Estimate/d1$Sum
d1$prev_2.5<- d1$Q2.5/d1$Sum
d1$prev_97.5.5<- d1$Q97.5/d1$Sum
d1$Risk<-"Risk2"
write.csv(d1,"Risk2_prev_R2.csv")

d1$Agegroup<- factor(d1$Agegroup, levels = c("<35","35-65",">=65"  ))
d1$risk_R2_c<- factor(d1$risk_R2_c ,levels = c("Not High","High","Missing"))
ggplot(d1[which(d1$R1_Result=="Negative"),], aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=risk_R2_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=risk_R2_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Percentage")+xlab("Agegroup")+
  labs(colour="Risk Perception (Winter)")+theme(legend.position="bottom")+
  facet_grid(.~Sex, scales = "free")+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggsave("agegroup vs sex vs risk_R2 vs R1res.png")

# Inf Risk R1

d<- reshape2::dcast(bdat, Sex+Agegroup+Inf_risk_R1_c~status)
d[is.na(d)]<-"Missing"
d$Sum<- d$`Never Positive`+d$`Atleast once Positive`
d$Pos<- d$`Atleast once Positive`

d$Inf_risk_R1_c<- factor(d$Inf_risk_R1_c,levels = c("High","Low","Missing"))

fit1 <- brm(
  Pos|trials(Sum)  ~ Sex*Inf_risk_R1_c+Agegroup*Inf_risk_R1_c, data =d, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit1, 
          type = "trace")

d<-cbind(d,fitted(fit1))
d$prev_est<- d$Estimate/d$Sum
d$prev_2.5<- d$Q2.5/d$Sum
d$prev_97.5.5<- d$Q97.5/d$Sum

d$Risk<-"Inf_Risk1"
write.csv(d,"Inf_Risk1_prev.csv")
d$Agegroup<- factor(d$Agegroup, levels = c("<35","35-65",">=65"  ))
d$Inf_risk_R1_c<- factor(d$Inf_risk_R1_c ,levels = c("Low","High","Missing"))
ggplot(d, aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=Inf_risk_R1_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=Inf_risk_R1_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Prevalence")+xlab("Agegroup")+
  labs(colour="Infection Risk (Summer)")+theme(legend.position="bottom")+
  facet_grid(.~Sex)+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggsave("agegroup vs sex vs inf_riskR1.png")

# Inf Risk R1 for R2 result

d<- reshape2::dcast(bdat_l, Sex+Agegroup+Inf_risk_R1_c + R1_Result~R2_Result)
d[is.na(d)]<-"Missing"
d$Sum<- d$Positive+d$Negative
d$Pos<- d$Positive

d$Inf_risk_R1_c<- factor(d$Inf_risk_R1_c,levels = c("High","Low","Missing"))

fit1 <- brm(
  Pos|trials(Sum)  ~ Sex*Inf_risk_R1_c+Agegroup*Inf_risk_R1_c + R1_Result, data =d, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit1, 
          type = "trace")

d<-cbind(d,fitted(fit1))
d$prev_est<- d$Estimate/d$Sum
d$prev_2.5<- d$Q2.5/d$Sum
d$prev_97.5.5<- d$Q97.5/d$Sum

d$Risk<-"Inf_Risk1"
write.csv(d,"Inf_Risk1_prev_R2.csv")
d$Agegroup<- factor(d$Agegroup, levels = c("<35","35-65",">=65"  ))
d$Inf_risk_R1_c<- factor(d$Inf_risk_R1_c ,levels = c("Low","High","Missing"))
ggplot(d[which(d$R1_Result=="Negative"),], aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=Inf_risk_R1_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=Inf_risk_R1_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Percentage")+xlab("Agegroup")+
  labs(colour="Infection Risk (Summer)")+theme(legend.position="bottom")+
  facet_grid(R1_Result~Sex, scales="free")+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggsave("agegroup vs sex vs inf_riskR1 vs R1Res.png")

# Inf_Risk R2
d1<- reshape2::dcast(bdat, Sex+Agegroup+Inf_risk_R2_c~status)
d1[is.na(d1)]<-"Missing"
d1$Sum<- d1$`Never Positive`+d1$`Atleast once Positive`
d1$Pos<- d1$`Atleast once Positive`


d1$Inf_risk_R2_c<- factor(d1$Inf_risk_R2_c,levels = c("High","Low","Missing"))

fit2 <- brm(
  Pos|trials(Sum)  ~ Sex*Inf_risk_R2_c+Agegroup*Inf_risk_R2_c, data =d1, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit2, 
          type = "trace")


d1<-cbind(d1,fitted(fit2))
d1$prev_est<- d1$Estimate/d1$Sum
d1$prev_2.5<- d1$Q2.5/d1$Sum
d1$prev_97.5.5<- d1$Q97.5/d1$Sum
d1$Risk<-"Inf_Risk2"
write.csv(d1,"Inf_Risk2_prev.csv")

d1$Agegroup<- factor(d1$Agegroup, levels = c("<35","35-65",">=65"  ))
d1$Inf_risk_R2_c<- factor(d1$Inf_risk_R2_c ,levels = c("Low","High","Missing"))
ggplot(d1, aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=Inf_risk_R2_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=Inf_risk_R2_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Prevalence")+xlab("Agegroup")+
  labs(colour="Infection Risk (Winter)")+theme(legend.position="bottom")+
  facet_grid(.~Sex)+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggsave("agegroup vs sex vs Inf_risk_R2.png")

# Inf_Risk R2 for R2 Result
d1<- reshape2::dcast(bdat_l, Sex+Agegroup+Inf_risk_R2_c + R1_Result~R2_Result)
d1[is.na(d1)]<-"Missing"
d1$Sum<- d1$Positive+d1$Negative
d1$Pos<- d1$Positive


d1$Inf_risk_R2_c<- factor(d1$Inf_risk_R2_c,levels = c("High","Low","Missing"))

fit2 <- brm(
  Pos|trials(Sum)  ~ Sex*Inf_risk_R2_c+Agegroup*Inf_risk_R2_c +R1_Result, data =d1, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit2, 
          type = "trace")


d1<-cbind(d1,fitted(fit2))
d1$prev_est<- d1$Estimate/d1$Sum
d1$prev_2.5<- d1$Q2.5/d1$Sum
d1$prev_97.5.5<- d1$Q97.5/d1$Sum
d1$Risk<-"Inf_Risk2"
write.csv(d1,"Inf_Risk2_prev_R2.csv")

d1$Agegroup<- factor(d1$Agegroup, levels = c("<35","35-65",">=65"  ))
d1$Inf_risk_R2_c<- factor(d1$Inf_risk_R2_c ,levels = c("Low","High","Missing"))
ggplot(d1[which(d1$R1_Result=="Negative"),], aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=Inf_risk_R2_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=Inf_risk_R2_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Percentage")+xlab("Agegroup")+
  labs(colour="Infection Risk (Winter)")+theme(legend.position="bottom")+
  facet_grid(.~Sex,scales="free")+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggsave("agegroup vs sex vs Inf_risk_R2 vs R1Result.png")

# Inf Grade R1

d<- reshape2::dcast(bdat, Sex+Agegroup+Inf_grade_R1_c~status)
d[is.na(d)]<-"Missing"
d$Sum<- d$`Never Positive`+d$`Atleast once Positive`
d$Pos<- d$`Atleast once Positive`

d$Inf_grade_R1_c<- factor(d$Inf_grade_R1_c,levels = c("High","Low","Missing"))

fit1 <- brm(
  Pos|trials(Sum)  ~ Sex*Inf_grade_R1_c+Agegroup*Inf_grade_R1_c, data =d, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit1, 
          type = "trace")

d<-cbind(d,fitted(fit1))
d$prev_est<- d$Estimate/d$Sum
d$prev_2.5<- d$Q2.5/d$Sum
d$prev_97.5.5<- d$Q97.5/d$Sum

d$Risk<-"Inf_Grade1"
write.csv(d,"Inf_Grade1_prev.csv")

d$Agegroup<- factor(d$Agegroup, levels = c("<35","35-65",">=65"  ))
d$Inf_grade_R1_c<- factor(d$Inf_grade_R1_c ,levels = c("Low","High","Missing"))
ggplot(d, aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=Inf_grade_R1_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=Inf_grade_R1_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Prevalence")+xlab("Agegroup")+
  labs(colour="Infection Grade (Summer)")+theme(legend.position="bottom")+
  facet_grid(.~Sex)+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggsave("agegroup vs sex vs inf_gradeR1.png")

# Inf Grade R1 for R2 Result

d<- reshape2::dcast(bdat_l, Sex+Agegroup+Inf_grade_R1_c +R1_Result~R2_Result)
d[is.na(d)]<-"Missing"
d$Sum<- d$Positive+d$Negative
d$Pos<- d$Positive

d$Inf_grade_R1_c<- factor(d$Inf_grade_R1_c,levels = c("High","Low","Missing"))

fit1 <- brm(
  Pos|trials(Sum)  ~ Sex*Inf_grade_R1_c+Agegroup*Inf_grade_R1_c + R1_Result, data =d, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit1, 
          type = "trace")

d<-cbind(d,fitted(fit1))
d$prev_est<- d$Estimate/d$Sum
d$prev_2.5<- d$Q2.5/d$Sum
d$prev_97.5.5<- d$Q97.5/d$Sum

d$Risk<-"Inf_Grade1"
write.csv(d,"Inf_Grade1_prev_R2.csv")

d$Agegroup<- factor(d$Agegroup, levels = c("<35","35-65",">=65"  ))
d$Inf_grade_R1_c<- factor(d$Inf_grade_R1_c ,levels = c("Low","High","Missing"))
ggplot(d[which(d$R1_Result=="Negative"),], aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=Inf_grade_R1_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=Inf_grade_R1_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Percentage")+xlab("Agegroup")+
  labs(colour="Infection Grade (Summer)")+theme(legend.position="bottom")+
  facet_grid(.~Sex,scales="free")+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggsave("agegroup vs sex vs inf_gradeR1 vs R1Res.png")


# Inf_Grade R2
d1<- reshape2::dcast(bdat, Sex+Agegroup+Inf_grade_R2_c~status)
d1[is.na(d1)]<-"Missing"
d1$Sum<- d1$`Never Positive`+d1$`Atleast once Positive`
d1$Pos<- d1$`Atleast once Positive`


d1$Inf_grade_R2_c<- factor(d1$Inf_grade_R2_c,levels = c("High","Low","Missing"))

fit2 <- brm(
  Pos|trials(Sum)  ~ Sex*Inf_grade_R2_c+Agegroup*Inf_grade_R2_c, data =d1, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit2, 
          type = "trace")


d1<-cbind(d1,fitted(fit2))
d1$prev_est<- d1$Estimate/d1$Sum
d1$prev_2.5<- d1$Q2.5/d1$Sum
d1$prev_97.5.5<- d1$Q97.5/d1$Sum
d1$Risk<-"Inf_grade2"
write.csv(d1,"Inf_grade2_prev.csv")

d1$Agegroup<- factor(d1$Agegroup, levels = c("<35","35-65",">=65"  ))
d1$Inf_grade_R2_c<- factor(d1$Inf_grade_R2_c ,levels = c("Low","High","Missing"))
ggplot(d1, aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=Inf_grade_R2_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=Inf_grade_R2_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Prevalence")+xlab("Agegroup")+
  labs(colour="Infection Grade (Winter)")+theme(legend.position="bottom")+
  facet_grid(.~Sex)+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggsave("agegroup vs sex vs Inf_grade_R2.png")


# Inf_Grade R2 for R2 Result
d1<- reshape2::dcast(bdat_l, Sex+Agegroup+Inf_grade_R2_c + R1_Result~R2_Result)
d1[is.na(d1)]<-"Missing"
d1$Sum<- d1$Positive+d1$Negative
d1$Pos<- d1$Positive


d1$Inf_grade_R2_c<- factor(d1$Inf_grade_R2_c,levels = c("High","Low","Missing"))

fit2 <- brm(
  Pos|trials(Sum)  ~ Sex*Inf_grade_R2_c+Agegroup*Inf_grade_R2_c + R1_Result, data =d1, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit2, 
          type = "trace")


d1<-cbind(d1,fitted(fit2))
d1$prev_est<- d1$Estimate/d1$Sum
d1$prev_2.5<- d1$Q2.5/d1$Sum
d1$prev_97.5.5<- d1$Q97.5/d1$Sum
d1$Risk<-"Inf_grade2"
write.csv(d1,"Inf_grade2_prev_R2.csv")

d1$Agegroup<- factor(d1$Agegroup, levels = c("<35","35-65",">=65"  ))
d1$Inf_grade_R2_c<- factor(d1$Inf_grade_R2_c ,levels = c("Low","High","Missing"))
ggplot(d1[which(d1$R1_Result=="Negative"),], aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=Inf_grade_R2_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=Inf_grade_R2_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Percentage")+xlab("Agegroup")+
  labs(colour="Infection Grade (Winter)")+theme(legend.position="bottom")+
  facet_grid(.~Sex,scales="free")+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggsave("agegroup vs sex vs Inf_grade_R2 vs R1Res.png")


# AKT_C

d<- reshape2::dcast(bdat, Sex+Agegroup+AKT_c~status)
d[is.na(d)]<-"Missing"
d$Sum<- d$`Never Positive`+d$`Atleast once Positive`
d$Pos<- d$`Atleast once Positive`

d$AKT_c<- factor(d$AKT_c,levels = c("<6",">=6","Missing"))

fit1 <- brm(
  Pos|trials(Sum)  ~ Sex*AKT_c+Agegroup*AKT_c, data =d, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit1, 
          type = "trace")

d<-cbind(d,fitted(fit1))
d$prev_est<- d$Estimate/d$Sum
d$prev_2.5<- d$Q2.5/d$Sum
d$prev_97.5.5<- d$Q97.5/d$Sum

d$Risk<-"AKT_C"
write.csv(d,"AKT_c_prev.csv")


d$Agegroup<- factor(d$Agegroup, levels = c("<35","35-65",">=65"  ))
d$AKT_c<- factor(d$AKT_c ,levels = c("<6",">=6","Missing"))
ggplot(d, aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=AKT_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=AKT_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Prevalence")+xlab("Agegroup")+
  labs(colour="Activity level in Summer 2020")+theme(legend.position="bottom")+
  facet_grid(.~Sex)+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggsave("agegroup vs sex vs activityR2.png")

# AKT_C for R2 Result

d<- reshape2::dcast(bdat_l, Sex+Agegroup+AKT_c+R1_Result~R2_Result)
d[is.na(d)]<-"Missing"
d$Sum<- d$Positive+d$Negative
d$Pos<- d$Positive

d$AKT_c<- factor(d$AKT_c,levels = c("<6",">=6","Missing"))

fit1 <- brm(
  Pos|trials(Sum)  ~ Sex*AKT_c+Agegroup*AKT_c + R1_Result, data =d, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit1, 
          type = "trace")

d<-cbind(d,fitted(fit1))
d$prev_est<- d$Estimate/d$Sum
d$prev_2.5<- d$Q2.5/d$Sum
d$prev_97.5.5<- d$Q97.5/d$Sum

d$Risk<-"AKT_C"
write.csv(d,"AKT_c_prev_R2.csv")


d$Agegroup<- factor(d$Agegroup, levels = c("<35","35-65",">=65"  ))
d$AKT_c<- factor(d$AKT_c ,levels = c("<6",">=6","Missing"))
f2<-ggplot(d[which(d$R1_Result=="Negative"),], aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=AKT_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=AKT_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Percentage")+xlab("Agegroup")+
  labs(colour="Activity Level (2021)")+theme(legend.position="bottom")+
  labs(colour="Activity level (Summer 2020)")+
  scale_color_discrete(labels=c("Not High","High","Missing"))+
  facet_grid(.~Sex,scales="free")+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
f2
ggsave("agegroup vs sex vs activityR2 vs R1Res.png")


# FEB_c
d1<- reshape2::dcast(bdat, Sex+Agegroup+FEB_c~status)
d1[is.na(d1)]<-"Missing"
d1$Sum<- d1$`Never Positive`+d1$`Atleast once Positive`
d1$Pos<- d1$`Atleast once Positive`


d1$FEB_c<- factor(d1$FEB_c,levels = c("<12",">=12","Missing"))

fit2 <- brm(
  Pos|trials(Sum)  ~ Sex*FEB_c+Agegroup*FEB_c, data =d1, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit2, 
          type = "trace")


d1<-cbind(d1,fitted(fit2))
d1$prev_est<- d1$Estimate/d1$Sum
d1$prev_2.5<- d1$Q2.5/d1$Sum
d1$prev_97.5.5<- d1$Q97.5/d1$Sum
d1$Risk<-"FEB_c"
write.csv(d1,"FEB_c_prev.csv")

d1$Agegroup<- factor(d1$Agegroup, levels = c("<35","35-65",">=65"  ))
d1$FEB_c<- factor(d1$FEB_c,levels = c("<12",">=12","Missing"))
ggplot(d1, aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=FEB_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=FEB_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Prevalence")+xlab("Agegroup")+
  labs(colour="Activity Level (2020 February)")+theme(legend.position="bottom")+
  facet_grid(.~Sex)+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggsave("agegroup vs sex vs activityFeb2020.png")

# FEB_c for R2 Result
d1<- reshape2::dcast(bdat_l, Sex+Agegroup+FEB_c+R1_Result~R2_Result)
d1[is.na(d1)]<-"Missing"
d1$Sum<- d1$Positive+d1$Negative
d1$Pos<- d1$Positive


d1$FEB_c<- factor(d1$FEB_c,levels = c("<12",">=12","Missing"))

fit2 <- brm(
  Pos|trials(Sum)  ~ Sex*FEB_c+Agegroup*FEB_c+R1_Result, data =d1, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit2, 
          type = "trace")


d1<-cbind(d1,fitted(fit2))
d1$prev_est<- d1$Estimate/d1$Sum
d1$prev_2.5<- d1$Q2.5/d1$Sum
d1$prev_97.5.5<- d1$Q97.5/d1$Sum
d1$Risk<-"FEB_c"
write.csv(d1,"FEB_c_prevR2.csv")

d1$Agegroup<- factor(d1$Agegroup, levels = c("<35","35-65",">=65"  ))
d1$FEB_c<- factor(d1$FEB_c,levels = c("<12",">=12","Missing"))
ggplot(d1[which(d1$R1_Result=="Negative"),], aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=FEB_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=FEB_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Percentage")+xlab("Agegroup")+
  labs(colour="Activity Level (2020 February)")+theme(legend.position="bottom")+
  facet_grid(.~Sex,scales = "free")+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggsave("agegroup vs sex vs activityFeb2020 vs R1Res.png")

# Sum Contacts R1

bdat$sum_1_c<- ifelse(bdat$sum_contact_R1>9,"High","Not High")
bdat$sum_2_c<- ifelse(bdat$sum_contact_R2>9,"High","Not High")

d<- reshape2::dcast(bdat, Sex+Agegroup+sum_1_c~status)
d[is.na(d)]<-"Missing"
d$Sum<- d$`Never Positive`+d$`Atleast once Positive`
d$Pos<- d$`Atleast once Positive`

d$sum_1_c<- factor(d$sum_1_c,levels = c("High","Not High","Missing"))

fit1 <- brm(
  Pos|trials(Sum)  ~ Sex*sum_1_c+Agegroup*sum_1_c, data =d, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit1, 
          type = "trace")

d<-cbind(d,fitted(fit1))
d$prev_est<- d$Estimate/d$Sum
d$prev_2.5<- d$Q2.5/d$Sum
d$prev_97.5.5<- d$Q97.5/d$Sum

d$Risk<-"Sum of Contacts R1"
write.csv(d,"sum_contacts_R1_prev.csv")

d$Agegroup<- factor(d$Agegroup, levels = c("<35","35-65",">=65"  ))
d$sum_1_c<- factor(d$sum_1_c ,levels = c("Not High","High","Missing"),
                    labels = c("<=9",">9","Missing"))
ggplot(d, aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=sum_1_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=sum_1_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Prevalence")+xlab("Agegroup")+
  labs(colour="Sum of contacts (Summer)")+theme(legend.position="bottom")+
  facet_grid(.~Sex)+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggsave("agegroup vs sex vs sumcontactsR1.png")


# Sum Contacts R1 for R2 Result

bdat_l$sum_1_c<- ifelse(bdat_l$sum_contact_R1>9,"High","Not High")
bdat_l$sum_2_c<- ifelse(bdat_l$sum_contact_R2>9,"High","Not High")

d<- reshape2::dcast(bdat_l, Sex+Agegroup+sum_1_c + R1_Result~R2_Result)
d[is.na(d)]<-"Missing"
d$Sum<- d$Positive+d$Negative
d$Pos<- d$Positive

d$sum_1_c<- factor(d$sum_1_c,levels = c("Not High","High","Missing"))

fit1 <- brm(
  Pos|trials(Sum)  ~ Sex*sum_1_c+Agegroup*sum_1_c + R1_Result, data =d, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit1, 
          type = "trace")

d<-cbind(d,fitted(fit1))
d$prev_est<- d$Estimate/d$Sum
d$prev_2.5<- d$Q2.5/d$Sum
d$prev_97.5.5<- d$Q97.5/d$Sum

d$Risk<-"Sum of Contacts R1"
write.csv(d,"sum_contacts_R1_prev_R2.csv")

d$Agegroup<- factor(d$Agegroup, levels = c("<35","35-65",">=65"  ))
d$sum_1_c<- factor(d$sum_1_c ,levels = c("Not High","High","Missing"),
                   labels = c("<=9",">9","Missing"))
f3<-ggplot(d[which(d$R1_Result=="Negative"),], aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=sum_1_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=sum_1_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Percentage")+xlab("Agegroup")+
  labs(colour="Sum of contacts (Summer 2020)")+theme(legend.position="bottom")+
  scale_color_discrete(labels=c("Not High","High","Missing"))+
  facet_grid(.~Sex,scales="free")+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
f3
ggsave("agegroup vs sex vs sumcontactsR1 vs R1_Result.png")

(f1 / f3 / f2)+plot_annotation(tag_levels = 'A')

ggsave("behav.plot.jpeg",dpi=300,height=9, width=7,units="in")
# Sum Contacts R2
d1<- reshape2::dcast(bdat, Sex+Agegroup+sum_2_c~status)
d1[is.na(d1)]<-"Missing"
d1$Sum<- d1$`Never Positive`+d1$`Atleast once Positive`
d1$Pos<- d1$`Atleast once Positive`


d1$sum_2_c<- factor(d1$sum_2_c,levels = c("High","Not High","Missing"))

fit2 <- brm(
  Pos|trials(Sum)  ~ Sex*sum_2_c+Agegroup*sum_2_c, data =d1, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit2, 
          type = "trace")


d1<-cbind(d1,fitted(fit2))
d1$prev_est<- d1$Estimate/d1$Sum
d1$prev_2.5<- d1$Q2.5/d1$Sum
d1$prev_97.5.5<- d1$Q97.5/d1$Sum
d1$Risk<-"Sum of Contacts R2"
write.csv(d1,"sum_contacts_R2_prev.csv")


d1$Agegroup<- factor(d1$Agegroup, levels = c("<35","35-65",">=65"  ))
d1$sum_2_c<- factor(d1$sum_2_c ,levels = c("Not High","High","Missing"),
                   labels = c("<=9",">9","Missing"))
ggplot(d1, aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=sum_2_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=sum_2_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Prevalence")+xlab("Agegroup")+
  labs(colour="Sum of contacts (Winter)")+theme(legend.position="bottom")+
  facet_grid(.~Sex)+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggsave("agegroup vs sex vs sumcontactsR2.png")


# Sum Contacts R2 for R2_Result
d1<- reshape2::dcast(bdat_l, Sex+Agegroup+sum_2_c + R1_Result~R2_Result)
d1[is.na(d1)]<-"Missing"
d1$Sum<- d1$Positive +d1$Negative
d1$Pos<- d1$Positive


d1$sum_2_c<- factor(d1$sum_2_c,levels = c("High","Not High","Missing"))

fit2 <- brm(
  Pos|trials(Sum)  ~ Sex*sum_2_c+Agegroup*sum_2_c +R1_Result, data =d1, 
  chains = 3,warmup = 5000, 
  iter = 10000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit2, 
          type = "trace")


d1<-cbind(d1,fitted(fit2))
d1$prev_est<- d1$Estimate/d1$Sum
d1$prev_2.5<- d1$Q2.5/d1$Sum
d1$prev_97.5.5<- d1$Q97.5/d1$Sum
d1$Risk<-"Sum of Contacts R2"
write.csv(d1,"sum_contacts_R2_prev_R2.csv")


d1$Agegroup<- factor(d1$Agegroup, levels = c("<35","35-65",">=65"  ))
d1$sum_2_c<- factor(d1$sum_2_c ,levels = c("Not High","High","Missing"),
                    labels = c("<=9",">9","Missing"))
ggplot(d1[which(d1$R1_Result=="Negative"),], aes(x=Agegroup,y=prev_est))+
  geom_point(aes(colour=sum_2_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=sum_2_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Percentage")+xlab("Agegroup")+
  labs(colour="Sum of contacts (Winter)")+theme(legend.position="bottom")+
  facet_grid(.~Sex,scales = "free")+scale_y_continuous(labels = scales::percent_format(accuracy = 1))
ggsave("agegroup vs sex vs sumcontactsR2 vs R1Res.png")


# merge data on socio economics and smoking to prevalence
bdat<- merge.data.frame(bdat, koco_r2_dat[c("ind_id","NetIncome_monthly_1","smokestatus")],
                        by="ind_id")
bdat_l<- merge.data.frame(bdat_l, koco_r2_dat[c("ind_id","NetIncome_monthly_1","smokestatus")],
                        by="ind_id")


bdat$NetIncome_monthly_2<- ifelse(bdat$NetIncome_monthly_1 =="<=2500"|
                                    bdat$NetIncome_monthly_1=="2500-4000","<4000",">=4000")
bdat_l$NetIncome_monthly_2<- ifelse(bdat_l$NetIncome_monthly_1 =="<=2500"|
                                    bdat_l$NetIncome_monthly_1=="2500-4000","<4000",">=4000")

bdat$smokestatus1<- ifelse(bdat$smokestatus== "Current smoker","Current smoker","Exsmoker/Nonsmoker")
bdat_l$smokestatus1<- ifelse(bdat_l$smokestatus== "Current smoker","Current smoker","Exsmoker/Nonsmoker")

# Smoking and Income vs Overall Prevalence
d1<- reshape2::dcast(bdat,NetIncome_monthly_2+smokestatus1~status)
d1[is.na(d1)]<-"Missing"
d1$Sum<- d1$`Never Positive`+d1$`Atleast once Positive`
d1$Pos<- d1$`Atleast once Positive`

fit2 <- brm(
  Pos|trials(Sum)  ~ NetIncome_monthly_2*smokestatus1, data =d1, 
  chains = 3,warmup = 8000, 
  iter = 20000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit2, 
          type = "trace")


d1<-cbind(d1,fitted(fit2))
d1$prev_est<- d1$Estimate/d1$Sum
d1$prev_2.5<- d1$Q2.5/d1$Sum
d1$prev_97.5.5<- d1$Q97.5/d1$Sum
d1$Risk<-"Smoking_Income"
write.csv(d1,"smoking_income_prev.csv")

d1$smokestatus1<- factor(d1$smokestatus1, levels = c("Exsmoker/Nonsmoker","Current smoker","Missing"),
                        labels = c("Ex-smoker/\nNonsmoker","Current smoker","Missing"))
d1$NetIncome_monthly_2<- factor(d1$NetIncome_monthly_2, levels = c("<4000",">=4000","Missing"))
ggplot(d1, aes(x=smokestatus1,y=prev_est))+
  geom_point(aes(colour=NetIncome_monthly_2),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=NetIncome_monthly_2), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Prevalence")+xlab("Smoking status")+
  labs(colour="Monthly Household Income")+theme(legend.position="bottom")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

ggsave("smoking vs householdincome.png")

# Smoking vs Income on R2 result
d1<- reshape2::dcast(bdat_l,NetIncome_monthly_2+smokestatus1 + R1_Result~R2_Result)
d1[is.na(d1)]<-"Missing"
d1$Sum<- d1$Positive+d1$Negative
d1$Pos<- d1$Positive

fit2 <- brm(
  Pos|trials(Sum)  ~ NetIncome_monthly_2*smokestatus1 + R1_Result, data =d1, 
  chains = 3,warmup = 8000, 
  iter = 20000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit2, 
          type = "trace")


d1<-cbind(d1,fitted(fit2))
d1$prev_est<- d1$Estimate/d1$Sum
d1$prev_2.5<- d1$Q2.5/d1$Sum
d1$prev_97.5.5<- d1$Q97.5/d1$Sum
d1$Risk<-"Smoking_Income"
write.csv(d1,"smoking_income_prev_R2.csv")

d1$smokestatus1<- factor(d1$smokestatus1, levels = c("Exsmoker/Nonsmoker","Current smoker","Missing"),
                         labels = c("Ex-smoker/\nNonsmoker","Current smoker","Missing"))
d1$NetIncome_monthly_2<- factor(d1$NetIncome_monthly_2, levels = c("<4000",">=4000","Missing"))
ggplot(d1[which(d1$R1_Result=="Negative"),], aes(x=smokestatus1,y=prev_est))+
  geom_point(aes(colour=NetIncome_monthly_2),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=NetIncome_monthly_2), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Percentage")+xlab("Smoking status")+#facet_grid(R1_Result~.,scales = "free")+
  labs(colour="Monthly Household Income")+theme(legend.position="bottom")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

ggsave("smoking vs householdincome vs R1_Result.png")

d1<- reshape2::dcast(bdat,smokestatus+sum_1_c~status)
d1[is.na(d1)]<-"Missing"
d1$Sum<- d1$`Never Positive`+d1$`Atleast once Positive`
d1$Pos<- d1$`Atleast once Positive`

fit2 <- brm(
  Pos|trials(Sum)  ~ smokestatus*sum_1_c, data =d1, 
  chains = 3,warmup = 8000, 
  iter = 20000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit2, 
          type = "trace")


d1<-cbind(d1,fitted(fit2))
d1$prev_est<- d1$Estimate/d1$Sum
d1$prev_2.5<- d1$Q2.5/d1$Sum
d1$prev_97.5.5<- d1$Q97.5/d1$Sum
d1$Risk<-"Smoking_Contacts_R1"
write.csv(d1,"smoking_contacts_prev.csv")
d1$smokestatus<- factor(d1$smokestatus, levels = c("Non smoker","Past smoker","Current smoker","Missing"))
d1$sum_1_c<- factor(d1$sum_1_c, levels = c("Not High","High","Missing"))
ggplot(d1, aes(x=smokestatus,y=prev_est))+
  geom_point(aes(colour=sum_1_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=sum_1_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Prevalence")+xlab("Smoking status")+
  labs(colour="Sum of Contacts (Summer)")+theme(legend.position="bottom")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))
  
ggsave("smoking vs sum of contacts.png")


d1<- reshape2::dcast(bdat_l,smokestatus+sum_1_c + R1_Result~R2_Result)
d1[is.na(d1)]<-"Missing"
d1$Sum<- d1$Positive+d1$Negative
d1$Pos<- d1$Positive

fit2 <- brm(
  Pos|trials(Sum)  ~ smokestatus*sum_1_c + R1_Result, data =d1, 
  chains = 3,warmup = 8000, 
  iter = 20000,
  control = list(adapt_delta = 0.95),family = binomial,seed=1234
)
mcmc_plot(fit2, 
          type = "trace")


d1<-cbind(d1,fitted(fit2))
d1$prev_est<- d1$Estimate/d1$Sum
d1$prev_2.5<- d1$Q2.5/d1$Sum
d1$prev_97.5.5<- d1$Q97.5/d1$Sum
d1$Risk<-"Smoking_Contacts_R1"
write.csv(d1,"smoking_contacts_prev_R2.csv")
d1$smokestatus<- factor(d1$smokestatus, levels = c("Non smoker","Past smoker","Current smoker","Missing"))
d1$sum_1_c<- factor(d1$sum_1_c, levels = c("Not High","High","Missing"))
ggplot(d1[which(d1$R1_Result=="Negative"),], aes(x=smokestatus,y=prev_est))+
  geom_point(aes(colour=sum_1_c),size=1.25,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin = prev_2.5, ymax =prev_97.5.5, colour=sum_1_c), position = position_dodge(width = 0.9),
                width = 0.2)+
  ylab("Percentage")+xlab("Smoking status")+#facet_grid(R1_Result~., scales="free")+
  labs(colour="Sum of Contacts (Summer)")+theme(legend.position="bottom")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

ggsave("smoking vs sum of contacts vs R1_Result.png")


