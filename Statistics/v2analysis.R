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
# Read Datasets
dat <- read.csv("second_round_one_meas.csv",header = T)
d1<-dat


# select KoCo19 data

dat <- dat[which(dat$study_per_ID %in% c("KoCo19-Filterpaper-Nov","KoCo19-Tend-Filter-Nov")==TRUE),]

# subset DBS data
id_dat <- subset(dat,select = c(study_per_ID,HaushaltsId,tln_ID,blut_ID,Q1_b_Datum_Besuch_blut_1_2,
                               Q1_b_Datum_Besuch,date_prick_self,date_entry_TROPI,
                               Q6_d_Geschlecht,Q6_f1_Geburtsdatum_DDMONYYYY,submissionStartDate,
                               DBS_value_1,DBS_value_2,DBS_value_3,DBS_value_4,DBS_value_5,DBS_value_6,
                               finalDBS,result_DBS))

l <- length(id_dat$tln_ID)
min_DBS <- rep(0, times=l)
for(i in 1:l)
{ t <-  id_dat[i,12:17]
  s<- length(t[which(is.na(t))])
  
  min_DBS[i]<- ifelse(s<6, min(t,na.rm = T),NA)  
}

id_dat <- cbind(id_dat,min_DBS)

id_dat$diff<- id_dat$finalDBS - id_dat$min_DBS

id_dat$DBS_Result<- ifelse(id_dat$min_DBS <0.09,"Negative",
                           ifelse(id_dat$min_DBS>=0.09 &
                                    id_dat$min_DBS<0.12,"Intermediate","Positive"))
table(id_dat$DBS_Result)

# Make Date of Selp Prick if na to be 2 days lower than date of entry

id_dat$date_prick_self<- ifelse(is.na(id_dat$date_prick_self)==FALSE,as.Date(id_dat$date_prick_self, format = "%Y-%m-%d"),
                                as.Date(id_dat$date_entry_TROPI,format = "%Y-%m-%d") -2)
id_dat$date_prick_self<- as.Date(id_dat$date_prick_self, format="%Y-%m-%d")

table(id_dat$date_prick_self)

id_dat$week_prick<- strftime(id_dat$date_prick_self, format = "%V")

dbs_dat<- id_dat[c(2,3,4,7,8,20,22)]

colnames(dbs_dat)<-c("hh_id","ind_id","BloodID_DBS","Date_DBS_prick","Date_DBS_Tropi","DBS_quant","DBS_Result")



#dbs_dat$duplicated<- duplicated(dbs_dat$ind_id)


# People who couldnot do DBS
no_dbs<- d1[which(d1$odkForm %in% c("FU-KoCo19-4" ,"FU-KoCo19-5")==TRUE),]

ndbs<- subset(no_dbs,select = c(HaushaltsId,BewohnerId_mult,tln_ID,blut_ID,Q1_b_Datum_Besuch,
                                IgA_value_1,              
                                IgG_value_1,
                                Roche8000_value_1,
                                DBS_value_1,DBS_value_2,DBS_value_3))

#ndbs$match<- ifelse(ndbs$BewohnerId_mult==ndbs$tln_ID,"Match","Mismatch")

ndbs$new_DBS <- ifelse(ndbs$BewohnerId_mult %in% id_dat$tln_ID==TRUE,"Selfprick_Blood","Non_selfprick_blood") 

# new non self prick
new_dbs<- ndbs[which(ndbs$new_DBS=="Non_selfprick_blood"),]
new_dbs1<- new_dbs[c(1,2,4,5,9,10,11,12)]

l <- length(new_dbs1$BewohnerId_mult)
DBS_quant <- rep(0, times=l)
for(i in 1:l)
{ t <-  new_dbs1[i,c(5,6,7)]
  
  s<- length(t[which(is.na(t))])

  DBS_quant[i]<- ifelse(s<3, min(t,na.rm = T),NA)  
}

new_dbs1 <- cbind(new_dbs1,DBS_quant)
new_dbs1$DBS_Result<- ifelse(new_dbs1$DBS_quant <0.09,"Negative",
                             ifelse(new_dbs1$DBS_quant>=0.09 &
                                      new_dbs1$DBS_quant<0.12,"Intermediate","Positive"))
new_dbs1<- new_dbs1[which(!is.na(new_dbs1$DBS_quant)),]

new_dbs2<- new_dbs1[c(1,2,3,4,8,9,10)]

colnames(new_dbs2)<- c("hh_id","ind_id","BloodID_DBS","Date_DBS_prick",
                       "Pricktype" ,"DBS_quant" , "DBS_Result"  )
new_dbs2$Date_DBS_Tropi<- new_dbs2$Date_DBS_prick
new_dbs2$Pricktype<-"Tent_prick"

new_dbs2$Date_DBS_prick<- as.Date(new_dbs2$Date_DBS_prick, format="%Y-%m-%d")
new_dbs2$Date_DBS_Tropi<- as.Date(new_dbs2$Date_DBS_Tropi,format="%Y-%m-%d")
dbs_dat$Date_DBS_prick<- as.Date(dbs_dat$Date_DBS_prick,format="%Y-%m-%d")
dbs_dat$Date_DBS_Tropi<- as.Date(dbs_dat$Date_DBS_Tropi,format="%Y-%m-%d")

dbs_dat_new<- merge.data.frame(dbs_dat,new_dbs2, all = TRUE)

dbs_dat_new$Pricktype<- ifelse(is.na(dbs_dat_new$Pricktype),"Selfprick",dbs_dat_new$Pricktype)

# Blood samples DBS

blood<- ndbs[c(1,2,4,5,6,7,8,12)]

colnames(blood)<-c("hh_id","ind_id","BloodID_R2","Date_Blooddraw","IgA_quant_R2","IgG_quant_R2","Roche_quant_R2","VisitTent")

blood$Roche_R2_new<- ifelse(blood$Roche_quant_R2>=0.4218,"Positive","Negative")
blood$IgG_R2_new<- ifelse(blood$IgG_quant_R2>=1.015,"Positive","Negative")
blood$IgA_R2_new<- ifelse(blood$IgA_quant_R2>=1.085,"Positive","Negative")

blood<- blood[which(!is.na(blood$IgA_quant_R2)| !is.na(blood$IgG_quant_R2) |! is.na(blood$Roche_quant_R2)),]

blood$duplicated<- duplicated(blood$ind_id)

blood_dup<- blood[which(blood$duplicated=="TRUE"),]

blood2<- blood[which(blood$ind_id %in% blood_dup$ind_id==TRUE),]

blood<- blood[which(blood$duplicated=="FALSE"),]

dbs_dat_new1<- merge.data.frame(dbs_dat_new[which(!is.na(dbs_dat_new$ind_id)),],blood,by=c("ind_id","hh_id"),all.x = TRUE)


ctable(dbs_dat_new1$DBS_Result,dbs_dat_new1$Roche_R2_new,useNA = "no")

# Read KOCO Baseline Data
kb <- read.csv("KOCO_R1.csv",header = T)

kbase<- read.csv("Koco_baseline.csv",header = TRUE)
colnames(kb)[84]<-"IgA_Result_Old"

# Merge koco baseline data to dbs data

kb$duplicated<- duplicated(kb$ind_id,kb$hh_id)

koco_r2_dat<-merge.data.frame(kb[c(2:84)],dbs_dat_new1[which(!is.na(dbs_dat_new1$ind_id)),], by=c("ind_id","hh_id"), all.x = TRUE)


# Consider people having blood with blood results

koco_r2_dat$R2_Result<- ifelse(is.na(koco_r2_dat$DBS_Result) |
                                 koco_r2_dat$DBS_Result=="Intermediate",koco_r2_dat$Roche_R2_new,
                               koco_r2_dat$DBS_Result)

koco_r2_dat$date_R2<- ifelse( is.na(koco_r2_dat$DBS_Result) |
                               koco_r2_dat$DBS_Result=="Intermediate", as.Date(as.character(koco_r2_dat$Date_Blooddraw), format="%Y-%m-%d"),
                             as.Date(as.character(koco_r2_dat$Date_DBS_prick),format="%Y-%m-%d"))

koco_r2_dat$date_R2<- as.Date(koco_r2_dat$date_R2,format="%Y-%m-%d")

koco_r2_dat$R2_Result<- factor(koco_r2_dat$R2_Result)
table(koco_r2_dat$DBS_Result)
ctable(koco_r2_dat$Roche_R2_new, koco_r2_dat$DBS_Result)
ctable(koco_r2_dat$Roche_R2_new, koco_r2_dat$R2_Result)
ctable(koco_r2_dat$DBS_Result, koco_r2_dat$R2_Result)


koco_r2_complete<- koco_r2_dat[which(!is.na(koco_r2_dat$R2_Result) &
                                       koco_r2_dat$R2_Result!="Intermediate"),]
write.csv(koco_r2_dat,"koco_r2_dat.csv")
write.csv(koco_r2_complete,"koco_r2_matched.csv")
koco_r2_complete$Roche_Result_new<-factor(koco_r2_complete$Roche_Result_new)

table(koco_r2_complete$DBS_Result)
ctable(koco_r2_complete$Roche_R2_new, koco_r2_complete$DBS_Result)
ctable(koco_r2_complete$Roche_R2_new, koco_r2_complete$R2_Result)
ctable(koco_r2_complete$DBS_Result, koco_r2_complete$R2_Result)
ctable(koco_r2_complete$Roche_Result_new, koco_r2_complete$R2_Result)

# additional covariates
# Additional Covariates

dat_risk<-read.csv("d_risk.csv",header=TRUE)
dat_kon<-read.csv("d_kon.csv",header=TRUE)
dat_fz<-read.csv("d_fz.csv",header=TRUE)

dat_risk<-dat_risk[c(2,6,11,7,9,8,10,
                     15,20,16,18,17,19)]

dat_risk$risk_R1_c<- ifelse(dat_risk$risk_R1<=2,"Very Low",
                            ifelse(dat_risk$risk_R1>2 & dat_risk$risk_R1<=4,"Low",
                                   ifelse(dat_risk$risk_R1>4 & dat_risk$risk_R1<=6,"Moderate","High")))
dat_risk$risk_R2_c<- ifelse(dat_risk$risk_R2<=2,"Very Low",
                            ifelse(dat_risk$risk_R2>2 & dat_risk$risk_R2<=4,"Low",
                                   ifelse(dat_risk$risk_R2>4 & dat_risk$risk_R2<=6,"Moderate","High")))

dat_risk$Inf_risk_R1_c<- ifelse(dat_risk$Inf_risk_R1<4,"Low","High")
dat_risk$Inf_risk_R2_c<- ifelse(dat_risk$Inf_risk_R2<4,"Low","High")

dat_risk$Inf_grade_R1_c<- ifelse(dat_risk$Inf_grade_R1<4,"Low","High")
dat_risk$Inf_grade_R2_c<- ifelse(dat_risk$Inf_grade_R2<4,"Low","High")

dat_risk$risk_R1_c<- factor(dat_risk$risk_R1_c,levels = c("High","Low","Moderate","Very Low"),
                            labels = c("High","Moderate","Low","Very Low"))
dat_risk$Inf_risk_R1_c<- factor(dat_risk$Inf_risk_R1_c,levels = c("High","Low"),
                            labels = c("High","Low"))
dat_risk$Inf_grade_R1_c<- factor(dat_risk$Inf_grade_R1_c,levels = c("High","Low"),
                            labels = c("High","Low"))
dat_risk$risk_R2_c<- factor(dat_risk$risk_R2_c,levels = c("High","Low","Moderate","Very Low"),
                            labels = c("High","Moderate","Low","Very Low"))
dat_risk$Inf_risk_R2_c<- factor(dat_risk$Inf_risk_R2_c,levels = c("High","Low"),
                            labels = c("High","Low"))
dat_risk$Inf_grade_R2_c<- factor(dat_risk$Inf_grade_R2_c,levels = c("High","Low"),
                                labels = c("High","Low"))

dat_kon<- dat_kon[c(2,13,24)]
dat_fz<-dat_fz[c(2,3,4)]

dat_fz$AKT_c<- ifelse(dat_fz$FZ_AKT<6,"<6",">=6")
dat_fz$FEB_c<- ifelse(dat_fz$FZ_FEB<12,"<12",">=12")

koco_r2_complete<- merge.data.frame(koco_r2_complete,dat_risk,by="ind_id",all.x = TRUE)
koco_r2_complete<- merge.data.frame(koco_r2_complete,dat_kon,by="ind_id",all.x = TRUE)
koco_r2_complete<- merge.data.frame(koco_r2_complete,dat_fz,by="ind_id",all.x = TRUE)

koco_r2_dat<- merge.data.frame(koco_r2_dat,dat_risk,by="ind_id",all.x = TRUE)
koco_r2_dat<- merge.data.frame(koco_r2_dat,dat_kon,by="ind_id",all.x = TRUE)
koco_r2_dat<- merge.data.frame(koco_r2_dat,dat_fz,by="ind_id",all.x = TRUE)

# Basic Tables
my_controls <- tableby.control(
  test = T,
  total = T,
  numeric.test = "kwt", 
  ordered.test = "trend",
  numeric.stats = c("meansd", "medianq1q3", "range", "Nmiss2"),
  cat.stats = c("countpct", "Nmiss2"),
  stats.labels = list(
    meansd = "Mean (SD)",
    medianq1q3 = "Median (Q1, Q3)",
    range = "Min - Max",
    Nmiss2 = "Missing"
  )
)
t1<-tableby(R2_Result~   fe(Sex)+    chisq(Agegroup)+ fe(Birth_Country)+
              chisq(EducationLevel)+ chisq(EmploymentStatus)+
              chisq(smokestatus)+ fe(Employment_ind_lastyear)+
              fe(any_risk_emp)+ 
              chisq(Household_Type_1)+ chisq(NetIncome_monthly_1)+
              chisq(Housing_Type)+ chisq(LivingArea_perInhabitant_cat)+
              chisq(Household_Inhabitants)+
              chisq(Symptoms_past2w_Fever)+
              chisq(Symptoms_past2w_Chills)+
              chisq(Symptoms_past2w_Limbpain)+
              chisq(Symptoms_past2w_Runningnose)+
              chisq(Symptoms_past2w_Sorethroat)+
              chisq(Symptoms_past2w_Headache)+
              chisq(Symptoms_past2w_Drycough)+
              chisq(Symptoms_past2w_Breathlessness)+
              chisq(Symptoms_past2w_Fatigue)+
              chisq(Symptoms_past2w_diarrhoea)+
              chisq(any_symptoms)+
              chisq(Symptoms_past2w_SenseLoss.Smell.Taste.)+
              num_symptoms_pos+
              fe(Facemask_public)+ chisq(HealthStatus_Self_Overall_1)+
              fe(HealthStatus_Self_Allergies_Respiratory)+
              fe(HealthStatus_Self_Diabetes)+
              fe(HealthStatus_Self_LungDisease)+
              fe(HealthStatus_Self_CardiovascularDisease)+
              fe(HealthStatus_Self_Cancer)+
              fe(HealthStatus_Self_Obesity)+
              fe(HealthStatus_Self_Allergies_Skin)+
              fe(HealthStatus_Self_Autoimmune)+
              fe(HealthStatus_Self_HIV)+
              fe(HealthStatus_Self_OtherChronic)+
              fe(Medication_Immunosuppressive)+
              fe(Medication_Cortisone)+
              fe(Medication_Bloodpressure)+
              fe(Roche_Result_new)+
              fe(IgG_result_new)+
              fe(IgA_Result_new)+
              fe(Vaccination_Influenza2019)+
              risk_R1+fe(risk_R1_c)+
              Inf_risk_R1+ fe(Inf_risk_R1_c)+
              Inf_grade_R1+fe(Inf_grade_R1_c)+
              risk_R2+fe(risk_R2_c)+
              Inf_risk_R2+fe(Inf_risk_R2_c)+
              Inf_grade_R2+fe(Inf_grade_R2_c)+
              sum_contact_R1+sum_contact_R2+
              FZ_AKT+FZ_FEB+
              fe(AKT_c)+fe(FEB_c)
              
              ,data =koco_r2_complete,
            control = my_controls , simulate.p.value = TRUE , B=10000)

write2html(summary(t1, pfootnote=TRUE),"Table1.html")


ctable(koco_r2_complete$Roche_Result_new,koco_r2_complete$R2_Result)

ctable(koco_r2_complete$Roche_Result_new,koco_r2_complete$DBS_Result,useNA = "no")


koco_r2_complete$status<- ifelse(koco_r2_complete$Roche_Result_new =="Positive" &
                                   koco_r2_complete$R2_Result=="Positive","Pos-Pos",
                                 ifelse(koco_r2_complete$Roche_Result_new=="Positive" &
                                          koco_r2_complete$R2_Result=="Negative","Pos-Neg",
                                        ifelse(koco_r2_complete$Roche_Result_new=="Negative" &
                                                 koco_r2_complete$R2_Result=="Positive","Neg-Pos","Neg-Neg")))


t2<-tableby(status~   chisq(Sex)+    chisq(Agegroup)+ chisq(Birth_Country)+
              chisq(EducationLevel)+ chisq(EmploymentStatus)+
              chisq(smokestatus)+ chisq(Employment_ind_lastyear)+
              chisq(any_risk_emp)+ 
              chisq(Household_Type_1)+ chisq(NetIncome_monthly_1)+
              chisq(Housing_Type)+ chisq(LivingArea_perInhabitant_cat)+
              chisq(Household_Inhabitants)+
              chisq(Symptoms_past2w_Fever)+
              chisq(Symptoms_past2w_Chills)+
              chisq(Symptoms_past2w_Limbpain)+
              chisq(Symptoms_past2w_Runningnose)+
              chisq(Symptoms_past2w_Sorethroat)+
              chisq(Symptoms_past2w_Headache)+
              chisq(Symptoms_past2w_Drycough)+
              chisq(Symptoms_past2w_Breathlessness)+
              chisq(Symptoms_past2w_Fatigue)+
              chisq(Symptoms_past2w_diarrhoea)+
              chisq(any_symptoms)+
              chisq(Symptoms_past2w_SenseLoss.Smell.Taste.)+
              num_symptoms_pos+
              chisq(Facemask_public)+ chisq(HealthStatus_Self_Overall_1)+
              chisq(HealthStatus_Self_Allergies_Respiratory)+
              chisq(HealthStatus_Self_Diabetes)+
              chisq(HealthStatus_Self_LungDisease)+
              chisq(HealthStatus_Self_CardiovascularDisease)+
              chisq(HealthStatus_Self_Cancer)+
              chisq(HealthStatus_Self_Obesity)+
              chisq(HealthStatus_Self_Allergies_Skin)+
              chisq(HealthStatus_Self_Autoimmune)+
              chisq(HealthStatus_Self_HIV)+
              chisq(HealthStatus_Self_OtherChronic)+
              chisq(Medication_Immunosuppressive)+
              chisq(Medication_Cortisone)+
              chisq(Medication_Bloodpressure)+
              chisq(Roche_Result_new)+
              chisq(IgG_result_new)+
              chisq(IgA_Result_new)+
              chisq(Vaccination_Influenza2019)+
              risk_R1+chisq(risk_R1_c)+
              Inf_risk_R1+ chisq(Inf_risk_R1_c)+
              Inf_grade_R1+chisq(Inf_grade_R1_c)+
              risk_R2+chisq(risk_R2_c)+
              Inf_risk_R2+chisq(Inf_risk_R2_c)+
              Inf_grade_R2+chisq(Inf_grade_R2_c)+
              sum_contact_R1+sum_contact_R2+
              FZ_AKT+FZ_FEB+
              chisq(AKT_c)+chisq(FEB_c)
            ,data =koco_r2_complete,
            control = my_controls , simulate.p.value = TRUE , B=10000)

write2html(summary(t2, pfootnote=TRUE),"Table2.html")


koco_r2_complete$status1<- factor(koco_r2_complete$status,levels = c("Pos-Pos","Pos-Neg","Neg-Neg","Neg-Pos"),
                                  labels = c("Stable/Pos-Neg","Stable/Pos-Neg","Stable/Pos-Neg","New Positives"))

t3<-tableby(status1~   fe(Sex)+    chisq(Agegroup)+ fe(Birth_Country)+
              chisq(EducationLevel)+ chisq(EmploymentStatus)+
              chisq(smokestatus)+ fe(Employment_ind_lastyear)+
              fe(any_risk_emp)+ 
              chisq(Household_Type_1)+ chisq(NetIncome_monthly_1)+
              chisq(Housing_Type)+ chisq(LivingArea_perInhabitant_cat)+
              chisq(Household_Inhabitants)+
              fe(Symptoms_past2w_Fever)+
              fe(Symptoms_past2w_Chills)+
              fe(Symptoms_past2w_Limbpain)+
              fe(Symptoms_past2w_Runningnose)+
              fe(Symptoms_past2w_Sorethroat)+
              fe(Symptoms_past2w_Headache)+
              fe(Symptoms_past2w_Drycough)+
              fe(Symptoms_past2w_Breathlessness)+
              fe(Symptoms_past2w_Fatigue)+
              fe(Symptoms_past2w_diarrhoea)+
              fe(any_symptoms)+
              fe(Symptoms_past2w_SenseLoss.Smell.Taste.)+
              num_symptoms_pos+
              fe(Facemask_public)+ chisq(HealthStatus_Self_Overall_1)+
              fe(HealthStatus_Self_Allergies_Respiratory)+
              fe(HealthStatus_Self_Diabetes)+
              fe(HealthStatus_Self_LungDisease)+
              fe(HealthStatus_Self_CardiovascularDisease)+
              fe(HealthStatus_Self_Cancer)+
              fe(HealthStatus_Self_Obesity)+
              fe(HealthStatus_Self_Allergies_Skin)+
              fe(HealthStatus_Self_Autoimmune)+
              fe(HealthStatus_Self_HIV)+
              fe(HealthStatus_Self_OtherChronic)+
              fe(Medication_Immunosuppressive)+
              fe(Medication_Cortisone)+
              fe(Medication_Bloodpressure)+
              fe(Roche_Result_new)+
              fe(IgG_result_new)+
              fe(IgA_Result_new)+
              fe(Vaccination_Influenza2019)+
              risk_R1+fe(risk_R1_c)+
              Inf_risk_R1+ fe(Inf_risk_R1_c)+
              Inf_grade_R1+fe(Inf_grade_R1_c)+
              risk_R2+fe(risk_R2_c)+
              Inf_risk_R2+fe(Inf_risk_R2_c)+
              Inf_grade_R2+fe(Inf_grade_R2_c)+
              sum_contact_R1+sum_contact_R2+
              FZ_AKT+FZ_FEB+
              fe(AKT_c)+fe(FEB_c)
            ,data =koco_r2_complete,
            control = my_controls , simulate.p.value = TRUE , B=10000)

write2html(summary(t3, pfootnote=TRUE),"Table3.html")

koco_r2_complete$status2<- factor(koco_r2_complete$status,levels = c("Pos-Pos","Pos-Neg","Neg-Neg","Neg-Pos"),
                                  labels = c("Positive","Positive","Never Positive","Positive"))
t4<-tableby(status2~   fe(Sex)+    chisq(Agegroup)+ fe(Birth_Country)+
              chisq(EducationLevel)+ chisq(EmploymentStatus)+
              chisq(smokestatus)+ fe(Employment_ind_lastyear)+
              fe(any_risk_emp)+ 
              chisq(Household_Type_1)+ chisq(NetIncome_monthly_1)+
              chisq(Housing_Type)+ chisq(LivingArea_perInhabitant_cat)+
              chisq(Household_Inhabitants)+
              fe(Symptoms_past2w_Fever)+
              fe(Symptoms_past2w_Chills)+
              fe(Symptoms_past2w_Limbpain)+
              fe(Symptoms_past2w_Runningnose)+
              fe(Symptoms_past2w_Sorethroat)+
              fe(Symptoms_past2w_Headache)+
              fe(Symptoms_past2w_Drycough)+
              fe(Symptoms_past2w_Breathlessness)+
              fe(Symptoms_past2w_Fatigue)+
              fe(Symptoms_past2w_diarrhoea)+
              fe(any_symptoms)+
              fe(Symptoms_past2w_SenseLoss.Smell.Taste.)+
              num_symptoms_pos+
              fe(Facemask_public)+ chisq(HealthStatus_Self_Overall_1)+
              fe(HealthStatus_Self_Allergies_Respiratory)+
              fe(HealthStatus_Self_Diabetes)+
              fe(HealthStatus_Self_LungDisease)+
              fe(HealthStatus_Self_CardiovascularDisease)+
              fe(HealthStatus_Self_Cancer)+
              fe(HealthStatus_Self_Obesity)+
              fe(HealthStatus_Self_Allergies_Skin)+
              fe(HealthStatus_Self_Autoimmune)+
              fe(HealthStatus_Self_HIV)+
              fe(HealthStatus_Self_OtherChronic)+
              fe(Medication_Immunosuppressive)+
              fe(Medication_Cortisone)+
              fe(Medication_Bloodpressure)+
              fe(Roche_Result_new)+
              fe(IgG_result_new)+
              fe(IgA_Result_new)+
              fe(Vaccination_Influenza2019)+
              risk_R1+fe(risk_R1_c)+
              Inf_risk_R1+ fe(Inf_risk_R1_c)+
              Inf_grade_R1+fe(Inf_grade_R1_c)+
              risk_R2+fe(risk_R2_c)+
              Inf_risk_R2+fe(Inf_risk_R2_c)+
              Inf_grade_R2+fe(Inf_grade_R2_c)+
              sum_contact_R1+sum_contact_R2+
              FZ_AKT+FZ_FEB+
              fe(AKT_c)+fe(FEB_c)
            ,data =koco_r2_complete,
            control = my_controls , simulate.p.value = TRUE , B=10000)

write2html(summary(t4, pfootnote=TRUE),"TableEverPositive.html")


koco_r2_dat$statusmis<- ifelse(koco_r2_dat$R2_Result=="Intermediate" | is.na(koco_r2_dat$R2_Result),"Missing",
                               "Not Missing")
t5<-tableby(statusmis~   fe(Sex)+    chisq(Agegroup)+ fe(Birth_Country)+
              chisq(EducationLevel)+ chisq(EmploymentStatus)+
              chisq(smokestatus)+ fe(Employment_ind_lastyear)+
              fe(any_risk_emp)+ 
              chisq(Household_Type_1)+ chisq(NetIncome_monthly_1)+
              chisq(Housing_Type)+ chisq(LivingArea_perInhabitant_cat)+
              chisq(Household_Inhabitants)+
              fe(Symptoms_past2w_Fever)+
              fe(Symptoms_past2w_Chills)+
              fe(Symptoms_past2w_Limbpain)+
              fe(Symptoms_past2w_Runningnose)+
              fe(Symptoms_past2w_Sorethroat)+
              fe(Symptoms_past2w_Headache)+
              fe(Symptoms_past2w_Drycough)+
              fe(Symptoms_past2w_Breathlessness)+
              fe(Symptoms_past2w_Fatigue)+
              fe(Symptoms_past2w_diarrhoea)+
              fe(any_symptoms)+
              fe(Symptoms_past2w_SenseLoss.Smell.Taste.)+
              num_symptoms_pos+
              fe(Facemask_public)+ chisq(HealthStatus_Self_Overall_1)+
              fe(HealthStatus_Self_Allergies_Respiratory)+
              fe(HealthStatus_Self_Diabetes)+
              fe(HealthStatus_Self_LungDisease)+
              fe(HealthStatus_Self_CardiovascularDisease)+
              fe(HealthStatus_Self_Cancer)+
              fe(HealthStatus_Self_Obesity)+
              fe(HealthStatus_Self_Allergies_Skin)+
              fe(HealthStatus_Self_Autoimmune)+
              fe(HealthStatus_Self_HIV)+
              fe(HealthStatus_Self_OtherChronic)+
              fe(Medication_Immunosuppressive)+
              fe(Medication_Cortisone)+
              fe(Medication_Bloodpressure)+
              fe(Roche_Result_new)+
              fe(IgG_result_new)+
              fe(IgA_Result_new)+
              fe(Vaccination_Influenza2019)+
              risk_R1+fe(risk_R1_c)+
              Inf_risk_R1+ fe(Inf_risk_R1_c)+
              Inf_grade_R1+fe(Inf_grade_R1_c)+
              risk_R2+fe(risk_R2_c)+
              Inf_risk_R2+fe(Inf_risk_R2_c)+
              Inf_grade_R2+fe(Inf_grade_R2_c)+
              sum_contact_R1+sum_contact_R2+
              FZ_AKT+FZ_FEB+
              fe(AKT_c)+fe(FEB_c)
            ,data =koco_r2_dat,
            control = my_controls , simulate.p.value = TRUE , B=10000)

write2html(summary(t5, pfootnote=TRUE),"TableMissing.html")


# Full Data for ever Positive
koco_r2_dat$status<- ifelse(is.na(koco_r2_dat$R2_Result)==TRUE &
                                    koco_r2_dat$Roche_Result_new=="Positive","Atleast once Positive",
                                  ifelse(is.na(koco_r2_dat$R2_Result)==TRUE &
                                           koco_r2_dat$Roche_Result_new=="Negative","Never Positive",
                                         ifelse(koco_r2_dat$R2_Result=="Positive" |
                                                  koco_r2_dat$Roche_Result_new=="Positive","Atleast once Positive","Never Positive")))

t6<-tableby(status~   fe(Sex)+    chisq(Agegroup)+ fe(Birth_Country)+
              chisq(EducationLevel)+ chisq(EmploymentStatus)+
              chisq(smokestatus)+ fe(Employment_ind_lastyear)+
              fe(any_risk_emp)+ 
              chisq(Household_Type_1)+ chisq(NetIncome_monthly_1)+
              chisq(Housing_Type)+ chisq(LivingArea_perInhabitant_cat)+
              chisq(Household_Inhabitants)+
              fe(Symptoms_past2w_Fever)+
              fe(Symptoms_past2w_Chills)+
              fe(Symptoms_past2w_Limbpain)+
              fe(Symptoms_past2w_Runningnose)+
              fe(Symptoms_past2w_Sorethroat)+
              fe(Symptoms_past2w_Headache)+
              fe(Symptoms_past2w_Drycough)+
              fe(Symptoms_past2w_Breathlessness)+
              fe(Symptoms_past2w_Fatigue)+
              fe(Symptoms_past2w_diarrhoea)+
              fe(any_symptoms)+
              fe(Symptoms_past2w_SenseLoss.Smell.Taste.)+
              num_symptoms_pos+
              fe(Facemask_public)+ chisq(HealthStatus_Self_Overall_1)+
              fe(HealthStatus_Self_Allergies_Respiratory)+
              fe(HealthStatus_Self_Diabetes)+
              fe(HealthStatus_Self_LungDisease)+
              fe(HealthStatus_Self_CardiovascularDisease)+
              fe(HealthStatus_Self_Cancer)+
              fe(HealthStatus_Self_Obesity)+
              fe(HealthStatus_Self_Allergies_Skin)+
              fe(HealthStatus_Self_Autoimmune)+
              fe(HealthStatus_Self_HIV)+
              fe(HealthStatus_Self_OtherChronic)+
              fe(Medication_Immunosuppressive)+
              fe(Medication_Cortisone)+
              fe(Medication_Bloodpressure)+
              fe(Roche_Result_new)+
              fe(IgG_result_new)+
              fe(IgA_Result_new)+
              fe(Vaccination_Influenza2019)+
              risk_R1+fe(risk_R1_c)+
              Inf_risk_R1+ fe(Inf_risk_R1_c)+
              Inf_grade_R1+fe(Inf_grade_R1_c)+
              risk_R2+fe(risk_R2_c)+
              Inf_risk_R2+fe(Inf_risk_R2_c)+
              Inf_grade_R2+fe(Inf_grade_R2_c)+
              sum_contact_R1+sum_contact_R2+
              FZ_AKT+FZ_FEB+
              fe(AKT_c)+fe(FEB_c)
            ,data =koco_r2_dat,
            control = my_controls , simulate.p.value = TRUE , B=10000)

write2html(summary(t6, pfootnote=TRUE),"TableEverPos_Full.html")

x<- subset(koco_r2_dat,select = c(ind_id,hh_id,Roche_Result_new,Roche_R2_new,DBS_Result,R2_Result,status,statusmis))


dbs_dat_new1$delay<- as.Date(dbs_dat_new1$Date_Blooddraw) -  as.Date(dbs_dat_new1$Date_DBS_prick)

dbs_dat_new1$delaycat<- ifelse(dbs_dat_new1$delay<30,"<30",
                               ifelse(dbs_dat_new1$delay>=30 & dbs_dat_new1$delay<60,"30-60",">=60"))

table(dbs_dat_new1$delaycat)
dbs_dat_new1$delaycat<- factor(dbs_dat_new1$delaycat, levels = c("<30","30-60",">=60"))
dbs_dat_new1$Roche_R2_new1<- ifelse(dbs_dat_new1$Roche_R2_new=="Positive","Positive","Negative")
dbs_dat_new1$cat<- ifelse(dbs_dat_new1$DBS_Result==dbs_dat_new1$Roche_R2_new1,"concordant","discordant")

ggplot(dbs_dat_new1[which(!is.na(dbs_dat_new1$Roche_quant_R2)),], aes(x=Roche_quant_R2,
                                                                      y=DBS_quant))+
  geom_point(aes(shape=delaycat, colour=cat))+
  geom_vline(xintercept = 0.4218)+
  geom_hline(yintercept = 0.09,linetype=2)+
  geom_hline(yintercept = 0.12,linetype=2)+
  scale_y_log10()+scale_x_log10()

reshape2::dcast(dbs_dat_new1[which(!is.na(dbs_dat_new1$Roche_quant_R2)),],
                delaycat+cat+Roche_R2_new~DBS_Result,margins = TRUE)  

reshape2::dcast(dbs_dat_new1[which(!is.na(dbs_dat_new1$Roche_quant_R2) &
                                     dbs_dat_new1$ind_id %in% koco_r2_dat$ind_id==TRUE),],
                delaycat+cat+Roche_R2_new~DBS_Result,margins = TRUE)  

disc_dat<-dbs_dat_new1[which(!is.na(dbs_dat_new1$Roche_quant_R2) &
                             dbs_dat_new1$ind_id %in% koco_r2_dat$ind_id==TRUE &
                               dbs_dat_new1$cat=="discordant"),]

base_disc<- subset( koco_r2_dat[which( koco_r2_dat$ind_id %in% disc_dat$ind_id==TRUE),], select = c(ind_id,Roche_Result_new,R_quant))
colnames(base_disc)<-c("ind_id","Roche_Baseline","Roche_qBaseline")

disc_dat<- merge.data.frame(base_disc,disc_dat,by="ind_id",all = TRUE)

disc_dat$Roche_R2_new1<-NULL

write.xlsx(disc_dat,"discordant_dat.xlsx")

koco_r2_dat$DBS_Result_new<- ifelse(koco_r2_dat$DBS_quant <0.09,"Negative",
                             ifelse(koco_r2_dat$DBS_quant>=0.09 &
                                      koco_r2_dat$DBS_quant<0.14,"Intermediate","Positive"))
koco_r2_dat$DBS_Result<- factor(koco_r2_dat$DBS_Result,levels = c("Negative","Intermediate","Positive"))
koco_r2_dat$DBS_Result_new<- factor(koco_r2_dat$DBS_Result_new,levels = c("Negative","Intermediate","Positive"))

ctable(koco_r2_dat$DBS_Result,koco_r2_dat$DBS_Result_new)

112/(112+4287)


koco_r2_dat$DBS_Result_new2<- ifelse(koco_r2_dat$DBS_quant <0.09,"Negative",
                                    ifelse(koco_r2_dat$DBS_quant>=0.09 &
                                             koco_r2_dat$DBS_quant<0.137,"Intermediate","Positive"))
koco_r2_dat$DBS_Result_new2<- factor(koco_r2_dat$DBS_Result_new2,levels = c("Negative","Intermediate","Positive"))
ctable(koco_r2_dat$DBS_Result,koco_r2_dat$DBS_Result_new2)


# Koco_V2_Long complete 2 measurements

koco_l<- subset(koco_r2_complete,select = 
                     c(ind_id,hh_id,
                       Sex,Agegroup,Age,Birth_Country,
                         EducationLevel,EmploymentStatus,
                         smokestatus,Employment_ind_lastyear,
                         any_risk_emp, 
                         Household_Type_1,NetIncome_monthly_1,
                         Housing_Type,LivingArea_perInhabitant_cat,
                         Household_Inhabitants,
                         Symptoms_past2w_Fever,
                         Symptoms_past2w_Chills,
                         Symptoms_past2w_Limbpain,
                         Symptoms_past2w_Runningnose,
                         Symptoms_past2w_Sorethroat,
                         Symptoms_past2w_Headache,
                         Symptoms_past2w_Drycough,
                         Symptoms_past2w_Breathlessness,
                         Symptoms_past2w_Fatigue,
                         Symptoms_past2w_diarrhoea,
                         any_symptoms,
                         Symptoms_past2w_SenseLoss.Smell.Taste.,
                         num_symptoms_pos,
                         Facemask_public,HealthStatus_Self_Overall_1,
                         HealthStatus_Self_Allergies_Respiratory,
                         HealthStatus_Self_Diabetes,
                         HealthStatus_Self_LungDisease,
                         HealthStatus_Self_CardiovascularDisease,
                         HealthStatus_Self_Cancer,
                         HealthStatus_Self_Obesity,
                         HealthStatus_Self_Allergies_Skin,
                         HealthStatus_Self_Autoimmune,
                         HealthStatus_Self_HIV,
                         HealthStatus_Self_OtherChronic,
                         Medication_Immunosuppressive,
                         Medication_Cortisone,
                         Medication_Bloodpressure,
                         Roche_Result_new,
                         IgG_result_new,
                         IgA_Result_new,VisitDate_Baseline,R2_Result, date_R2,
                         status1,status2,
                         risk_R1,risk_R1_c,
                         Inf_risk_R1,Inf_risk_R1_c,
                         Inf_grade_R1,Inf_grade_R1_c,
                         risk_R2,risk_R2_c,
                         Inf_risk_R2,Inf_risk_R2_c,
                         Inf_grade_R2,Inf_grade_R2_c,
                         sum_contact_R1,sum_contact_R2,
                         FZ_AKT,FZ_FEB,
                         AKT_c,FEB_c))
colnames(koco_l)[c(45,48)]<-c("R1_Result","date_R1")

koco_l$timeelapsed<-round(as.numeric( difftime(as.Date(koco_l$date_R2, format="%Y-%m-%d"),
                               as.Date(koco_l$date_R1, format="%Y-%m-%d"),units="weeks")),0)

koco_l1<- subset(koco_l, select = c(ind_id,hh_id,
                                    Sex,Agegroup,Age,Birth_Country,
                                    EducationLevel,EmploymentStatus,
                                    smokestatus,Employment_ind_lastyear,
                                    any_risk_emp, 
                                    Household_Type_1,NetIncome_monthly_1,
                                    Housing_Type,LivingArea_perInhabitant_cat,
                                    Household_Inhabitants,
                                    Symptoms_past2w_Fever,
                                    Symptoms_past2w_Chills,
                                    Symptoms_past2w_Limbpain,
                                    Symptoms_past2w_Runningnose,
                                    Symptoms_past2w_Sorethroat,
                                    Symptoms_past2w_Headache,
                                    Symptoms_past2w_Drycough,
                                    Symptoms_past2w_Breathlessness,
                                    Symptoms_past2w_Fatigue,
                                    Symptoms_past2w_diarrhoea,
                                    any_symptoms,
                                    Symptoms_past2w_SenseLoss.Smell.Taste.,
                                    num_symptoms_pos,
                                    Facemask_public,HealthStatus_Self_Overall_1,
                                    HealthStatus_Self_Allergies_Respiratory,
                                    HealthStatus_Self_Diabetes,
                                    HealthStatus_Self_LungDisease,
                                    HealthStatus_Self_CardiovascularDisease,
                                    HealthStatus_Self_Cancer,
                                    HealthStatus_Self_Obesity,
                                    HealthStatus_Self_Allergies_Skin,
                                    HealthStatus_Self_Autoimmune,
                                    HealthStatus_Self_HIV,
                                    HealthStatus_Self_OtherChronic,
                                    Medication_Immunosuppressive,
                                    Medication_Cortisone,
                                    Medication_Bloodpressure,
                                    IgG_result_new,
                                    IgA_Result_new,status1,status2,
                                    R1_Result,date_R1))
koco_l1$Visit<- 1
koco_l1$timeelapsed<-0

koco_l2<- subset(koco_l, select = c(ind_id,hh_id,
                                    Sex,Agegroup,Age,Birth_Country,
                                    EducationLevel,EmploymentStatus,
                                    smokestatus,Employment_ind_lastyear,
                                    any_risk_emp, 
                                    Household_Type_1,NetIncome_monthly_1,
                                    Housing_Type,LivingArea_perInhabitant_cat,
                                    Household_Inhabitants,
                                    Symptoms_past2w_Fever,
                                    Symptoms_past2w_Chills,
                                    Symptoms_past2w_Limbpain,
                                    Symptoms_past2w_Runningnose,
                                    Symptoms_past2w_Sorethroat,
                                    Symptoms_past2w_Headache,
                                    Symptoms_past2w_Drycough,
                                    Symptoms_past2w_Breathlessness,
                                    Symptoms_past2w_Fatigue,
                                    Symptoms_past2w_diarrhoea,
                                    any_symptoms,
                                    Symptoms_past2w_SenseLoss.Smell.Taste.,
                                    num_symptoms_pos,
                                    Facemask_public,HealthStatus_Self_Overall_1,
                                    HealthStatus_Self_Allergies_Respiratory,
                                    HealthStatus_Self_Diabetes,
                                    HealthStatus_Self_LungDisease,
                                    HealthStatus_Self_CardiovascularDisease,
                                    HealthStatus_Self_Cancer,
                                    HealthStatus_Self_Obesity,
                                    HealthStatus_Self_Allergies_Skin,
                                    HealthStatus_Self_Autoimmune,
                                    HealthStatus_Self_HIV,
                                    HealthStatus_Self_OtherChronic,
                                    Medication_Immunosuppressive,
                                    Medication_Cortisone,
                                    Medication_Bloodpressure,
                                    IgG_result_new,
                                    IgA_Result_new,status1,status2,
                                    R2_Result,date_R2, timeelapsed))
koco_l2$Visit<- 2

colnames(koco_l1)[c(49,50)]<-c("Result","VisitDate")
colnames(koco_l2)[c(49,50)]<-c("Result","VisitDate")

koco_long<- rbind.fill(koco_l1,koco_l2)

koco_long$Visit<- factor(koco_long$Visit,levels=c("1","2"),
                         labels=c("V1","V2"))
koco_long$Visit1<- as.numeric(koco_long$Visit)
koco_long$Result1<- as.numeric(koco_long$Result)-1

koco_l$R2_Result_1<- as.numeric(koco_l$R2_Result)-1
koco_l$R1_Result_1<- as.numeric(koco_l$R1_Result)-1

levels(koco_long$Result)
koco_long$id <- as.numeric(as.factor(koco_long$ind_id))
koco_long$id1 <- as.numeric(as.factor(koco_long$hh_id))
koco_l$Agegroup <- factor(koco_l$Agegroup)
koco_l$Agegroup<- relevel(koco_l$Agegroup,ref="20-34")
koco_l$smokestatus<- factor(koco_l$smokestatus)
koco_l$smokestatus<- relevel(koco_l$smokestatus,ref="Non smoker")
koco_l$NetIncome_monthly_1<- factor(koco_l$NetIncome_monthly_1)
koco_l$NetIncome_monthly_1<- relevel(koco_l$NetIncome_monthly_1,ref="2500-4000")
koco_l$Housing_Type<- factor(koco_l$Housing_Type, levels = c("Buildings with 1-2 apartments",
                                                             "Buildings with 3-4 apartments",
                                                             "Buildings with >=5 apartments",
                                                             "Others"))
koco_l$HealthStatus_Self_Overall_1<- factor(koco_l$HealthStatus_Self_Overall_1,
                                            levels = c("excellent"   ,"very good"   ,"good","not good" ),
                                            labels=c("excellent","very good","good or lower","good or lower"))


pos_hh<-reshape2::dcast(koco_r2_dat, hh_id~Roche_Result_new)
pos_hh$hh_base<- ifelse(pos_hh$Positive>0,"Not null","Null")

koco_l<- merge.data.frame(koco_l,pos_hh[c(1,4)],by="hh_id",all.x = TRUE)
koco_l$hh_base<- factor(koco_l$hh_base, levels = c("Null","Not null"))


koco_l$risk_R1_c<- factor(koco_l$risk_R1_c, levels = c("Very Low","Low","Moderate","High"))
koco_l$risk_R2_c<- factor(koco_l$risk_R2_c, levels = c("Very Low","Low","Moderate","High"))
koco_l$Inf_risk_R1_c<- factor(koco_l$Inf_risk_R1_c, levels = c("Low","High"))
koco_l$Inf_risk_R2_c<- factor(koco_l$Inf_risk_R2_c, levels = c("Low","High"))
koco_l$Inf_grade_R1_c<- factor(koco_l$Inf_grade_R1_c, levels = c("Low","High"))
koco_l$Inf_grade_R2_c<- factor(koco_l$Inf_grade_R2_c, levels = c("Low","High"))



koco_r2_dat$risk_R1_c<- factor(koco_r2_dat$risk_R1_c, levels = c("Very Low","Low","Moderate","High"))
koco_r2_dat$risk_R2_c<- factor(koco_r2_dat$risk_R2_c, levels = c("Very Low","Low","Moderate","High"))
koco_r2_dat$Inf_risk_R1_c<- factor(koco_r2_dat$Inf_risk_R1_c, levels = c("Low","High"))
koco_r2_dat$Inf_risk_R2_c<- factor(koco_r2_dat$Inf_risk_R2_c, levels = c("Low","High"))
koco_r2_dat$Inf_grade_R1_c<- factor(koco_r2_dat$Inf_grade_R1_c, levels = c("Low","High"))
koco_r2_dat$Inf_grade_R2_c<- factor(koco_r2_dat$Inf_grade_R2_c, levels = c("Low","High"))




# One way analysis

model_bin <- glmmTMB(Result~Visit+ (Visit|hh_id/ind_id),
                     family=binomial,data=koco_long)
summary(model_bin)


tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

library(bayesplot)


set_cmdstan_path(path = "C:/COVID Sentinel Kohorte Munich/KoCo19 Data/KOCO19V2")

cmdstan_path()

cmdstan_version()

m1<-brm(Result1~Visit+Sex+ (Visit|hh_id/ind_id),
    data = koco_long,
    chains = 3,warmup = 2000, 
    iter = 5000,
    control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m1, 
         type = "trace")


m2<-brm(Result1~timeelapsed+Sex+(timeelapsed|hh_id/ind_id),
        data = koco_long,
        chains = 3,warmup = 2000, 
        iter = 5000,
        control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m2, 
          type = "trace")


m3 <- brm(R2_Result_1~R1_Result*timeelapsed+ (1|hh_id),
                     data=koco_l,
              chains = 3,warmup = 2000, 
              iter = 5000,
              control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m3, 
          type = "trace")

mcmc_plot(m3, 
          type = "trace")

t1<-summary(m3)

exp(t1$fixed[,c(1,3,4)])

p_significance(m3)

m4 <- brm(R2_Result_1~R1_Result*timeelapsed+Sex+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m4, 
          type = "trace")


m5 <- brm(R2_Result_1~R1_Result+timeelapsed+Sex+ (1|hh_id),
          data=koco_l,
          chains = 3,warmup = 2000, 
          iter = 5000,
          control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m5, 
          type = "trace")


mcmc_plot(m5, 
          type = "trace")

t1<-summary(m5)

exp(t1$fixed[,c(1,3,4)])

p_significance(m5)

m6 <- brm(R2_Result_1~R1_Result+timeelapsed+ (1|hh_id),
          data=koco_l,
          chains = 3,warmup = 2000, 
          iter = 5000,
          control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m6, 
          type = "trace")

mcmc_plot(m6, 
          type = "trace")

t1<-summary(m6)

exp(t1$fixed[,c(1,3,4)])

p_significance(m6)
model_bin <- glmmTMB(R2_Result~R1_Result*timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

model_bin <- glmmTMB(R2_Result~R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)

x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

model_bin <- glmmTMB(R2_Result~Sex +R1_Result+timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)

x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))



model_bin <- glmmTMB(R2_Result~Agegroup +R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)

x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_3 <- brm(R2_Result~Agegroup +R1_Result+ timeelapsed+ (1|hh_id),
          data=koco_l,
          chains = 3,warmup = 2000, 
          iter = 5000,
          control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_3, 
          type = "trace")

t1<-summary(m_3)

exp(t1$fixed[,c(1,3,4)])
    
p_significance(m_3)

model_bin <- glmmTMB(R2_Result~ Agegroup + Sex+R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_4 <- brm(R2_Result~Agegroup + Sex+R1_Result+ timeelapsed+ (1|hh_id),
           data=koco_l,
           chains = 3,warmup = 2000, 
           iter = 5000,
           control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_4, 
          type = "trace")

t1<-summary(m_4)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_4)


model_bin <- glmmTMB(R2_Result~  Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_5 <- brm(R2_Result~Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
           data=koco_l,
           chains = 3,warmup = 2000, 
           iter = 5000,
           control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_5, 
          type = "trace")

t1<-summary(m_5)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_5)

model_bin <- glmmTMB(R2_Result ~ Birth_Country + Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_6 <- brm(R2_Result ~ Birth_Country + Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
           data=koco_l,
           chains = 3,warmup = 2000, 
           iter = 5000,
           control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_6, 
          type = "trace")

t1<-summary(m_6)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_6)

model_bin <- glmmTMB(R2_Result ~ EmploymentStatus+Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_7 <- brm(R2_Result ~ EmploymentStatus+Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
           data=koco_l,
           chains = 3,warmup = 2000, 
           iter = 5000,
           control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_7, 
          type = "trace")

t1<-summary(m_7)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_7)

model_bin <- glmmTMB(R2_Result ~ Employment_ind_lastyear+Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_8 <- brm(R2_Result ~ Employment_ind_lastyear+Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
           data=koco_l,
           chains = 3,warmup = 2000, 
           iter = 5000,
           control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_8, 
          type = "trace")

t1<-summary(m_8)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_8)

model_bin <- glmmTMB(R2_Result ~ EducationLevel+Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_9 <- brm(R2_Result ~ EducationLevel+Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
           data=koco_l,
           chains = 3,warmup = 2000, 
           iter = 5000,
           control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_9, 
          type = "trace")

t1<-summary(m_9)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_9)


model_bin <- glmmTMB(R2_Result ~ smokestatus+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))
m_10 <- brm(R2_Result ~ smokestatus+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
           data=koco_l,
           chains = 3,warmup = 2000, 
           iter = 5000,
           control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_10, 
          type = "trace")

t1<-summary(m_10)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_10)

model_bin <- glmmTMB(R2_Result ~ any_risk_emp+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_11 <- brm(R2_Result ~ any_risk_emp+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_11, 
          type = "trace")

t1<-summary(m_11)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_11)

model_bin <- glmmTMB(R2_Result ~ NetIncome_monthly_1+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_12 <- brm(R2_Result ~ NetIncome_monthly_1+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_12, 
          type = "trace")

t1<-summary(m_12)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_12)

model_bin <- glmmTMB(R2_Result ~ Household_Type_1+ Age+ Sex+ R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_13 <- brm(R2_Result ~ Household_Type_1+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_13, 
          type = "trace")

t1<-summary(m_13)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_13)

model_bin <- glmmTMB(R2_Result ~ Housing_Type+ Age+ Sex+  R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l[which(koco_l$Housing_Type!="Others"),])
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_14 <- brm(R2_Result ~ Housing_Type+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l[which(koco_l$Housing_Type!="Others"),],
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_14, 
          type = "trace")

t1<-summary(m_14)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_14)

model_bin <- glmmTMB(R2_Result ~ Household_Inhabitants+ Age+ Sex+ R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_15 <- brm(R2_Result ~ Household_Inhabitants+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_15, 
          type = "trace")

t1<-summary(m_15)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_15)
model_bin <- glmmTMB(R2_Result ~ LivingArea_perInhabitant_cat+ Age+ Sex+  R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_16 <- brm(R2_Result ~ LivingArea_perInhabitant_cat+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_16, 
          type = "trace")

t1<-summary(m_16)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_16)

model_bin <- glmmTMB(R2_Result ~ hh_base+ Age+ Sex+  R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_17 <- brm(R2_Result ~ hh_base+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_17, 
          type = "trace")

t1<-summary(m_17)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_17)


model_bin <- glmmTMB(R2_Result ~ Facemask_public+ Age+ Sex+ R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_18 <- brm(R2_Result ~Facemask_public+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_18, 
          type = "trace")

t1<-summary(m_18)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_18)



model_bin <- glmmTMB(R2_Result ~ Symptoms_past2w_SenseLoss.Smell.Taste.+ Age+ Sex+ R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))


m_19 <- brm(R2_Result ~Symptoms_past2w_SenseLoss.Smell.Taste.+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_19, 
          type = "trace")

t1<-summary(m_19)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_19)



model_bin <- glmmTMB(R2_Result ~Symptoms_past2w_Drycough + Age+ Sex+ R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_20 <- brm(R2_Result ~Symptoms_past2w_Drycough+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_20, 
          type = "trace")

t1<-summary(m_20)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_20)


model_bin <- glmmTMB(R2_Result ~Symptoms_past2w_Sorethroat + Age+ Sex+  R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_21 <- brm(R2_Result ~Symptoms_past2w_Sorethroat+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_21, 
          type = "trace")

t1<-summary(m_21)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_21)



model_bin <- glmmTMB(R2_Result ~Symptoms_past2w_Fever + Age+ Sex+  R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_22 <- brm(R2_Result ~Symptoms_past2w_Fever+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_22, 
          type = "trace")

t1<-summary(m_22)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_22)

model_bin <- glmmTMB(R2_Result ~Symptoms_past2w_Chills + Age+ Sex+ R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_23 <- brm(R2_Result ~Symptoms_past2w_Chills+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_23, 
          type = "trace")

t1<-summary(m_23)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_23)


model_bin <- glmmTMB(R2_Result ~Symptoms_past2w_Runningnose + Age+ Sex+ R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_24 <- brm(R2_Result ~Symptoms_past2w_Runningnose+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_24, 
          type = "trace")

t1<-summary(m_24)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_24)

model_bin <- glmmTMB(R2_Result ~num_symptoms_pos + Age+ Sex+ R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_25 <- brm(R2_Result ~num_symptoms_pos+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_25, 
          type = "trace")

t1<-summary(m_25)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_25)


model_bin <- glmmTMB(R2_Result ~any_symptoms + Age+ Sex+ R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_26 <- brm(R2_Result ~any_symptoms+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_26, 
          type = "trace")

t1<-summary(m_26)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_26)



model_bin <- glmmTMB(R2_Result ~ HealthStatus_Self_Overall_1 + Age+ Sex+  R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

m_27 <- brm(R2_Result ~HealthStatus_Self_Overall_1+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_27, 
          type = "trace")

t1<-summary(m_27)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_27)

cbind(x[,3],round(x[c(4,8,9,7)],3))
model_bin <- glmmTMB(R2_Result ~ HealthStatus_Self_Allergies_Respiratory + Age+ Sex+ R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_28 <- brm(R2_Result ~HealthStatus_Self_Allergies_Respiratory+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_28, 
          type = "trace")

t1<-summary(m_28)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_28)

model_bin <- glmmTMB(R2_Result ~ HealthStatus_Self_Allergies_Skin + Age+ Sex+ R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_29 <- brm(R2_Result ~HealthStatus_Self_Allergies_Skin+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_29, 
          type = "trace")

t1<-summary(m_29)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_29)


model_bin <- glmmTMB(R2_Result ~ HealthStatus_Self_Autoimmune + Age+ Sex+ R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))


m_30 <- brm(R2_Result ~ HealthStatus_Self_Autoimmune + Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_30, 
          type = "trace")

t1<-summary(m_30)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_30)

model_bin <- glmmTMB(R2_Result~HealthStatus_Self_Cancer + Age+ Sex+  R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_31 <- brm(R2_Result ~ HealthStatus_Self_Cancer + Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_31, 
          type = "trace")

t1<-summary(m_31)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_31)

model_bin <- glmmTMB(R2_Result ~HealthStatus_Self_CardiovascularDisease + Age+ Sex+ R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_32 <- brm(R2_Result ~ HealthStatus_Self_CardiovascularDisease+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_32, 
          type = "trace")

t1<-summary(m_32)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_32)

model_bin <- glmmTMB(R2_Result ~HealthStatus_Self_Diabetes + Age+ Sex+  R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_33 <- brm(R2_Result ~ HealthStatus_Self_Diabetes+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_33, 
          type = "trace")

t1<-summary(m_33)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_33)

model_bin <- glmmTMB(R2_Result ~HealthStatus_Self_LungDisease + Age+ Sex+  R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))

m_34 <- brm(R2_Result ~ HealthStatus_Self_LungDisease+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_34, 
          type = "trace")

t1<-summary(m_34)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_34)

model_bin <- glmmTMB(R2_Result ~HealthStatus_Self_Obesity + Age+ Sex+ R1_Result+ timeelapsed+ (1|hh_id),
                     family=binomial,data=koco_l)
summary(model_bin)
x<-tidy(model_bin ,conf.int=TRUE,exponentiate=TRUE,effects="fixed",conf.method="profile")

cbind(x[,3],round(x[c(4,8,9,7)],3))


m_35 <- brm(R2_Result ~ HealthStatus_Self_Obesity+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_35, 
          type = "trace")

t1<-summary(m_35)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_35)


m_36 <- brm(R2_Result ~ risk_R1_1+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_36, 
          type = "trace")

t1<-summary(m_36)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_36)



m_37 <- brm(R2_Result ~ risk_R1_c+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_37, 
          type = "trace")

t1<-summary(m_37)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_37)


m_38 <- brm(R2_Result ~ Inf_risk_R1_1+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_38, 
          type = "trace")

t1<-summary(m_38)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_38)

m_39 <- brm(R2_Result ~ Inf_risk_R1_c+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_39, 
          type = "trace")

t1<-summary(m_39)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_39)


m_40 <- brm(R2_Result ~ Inf_grade_R1_1+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_40, 
          type = "trace")

t1<-summary(m_40)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_40)

m_41 <- brm(R2_Result ~ Inf_grade_R1_c+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_41, 
          type = "trace")

t1<-summary(m_41)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_41)


m_42 <- brm(R2_Result ~ risk_R2_1+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_42, 
          type = "trace")

t1<-summary(m_42)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_42)
brm_pval(m_42)

m_43 <- brm(R2_Result ~ risk_R2_c+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_43, 
          type = "trace")

t1<-summary(m_43)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_43)


m_44 <- brm(R2_Result ~ Inf_risk_R2+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_44, 
          type = "trace")

t1<-summary(m_44)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_44)

m_45 <- brm(R2_Result ~ Inf_risk_R2_c+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_45, 
          type = "trace")

t1<-summary(m_45)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_45)


m_46 <- brm(R2_Result ~ Inf_grade_R2+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_46, 
          type = "trace")

t1<-summary(m_46)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_46)

m_47 <- brm(R2_Result ~ Inf_grade_R2_c+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_47, 
          type = "trace")

t1<-summary(m_47)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_47)


m_48 <- brm(R2_Result ~ sum_contact_R1+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_48, 
          type = "trace")

t1<-summary(m_48)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_48)

m_49 <- brm(R2_Result ~ sum_contact_R2+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_49, 
          type = "trace")

t1<-summary(m_49)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_49)


m_50 <- brm(R2_Result ~ FZ_AKT+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_50, 
          type = "trace")

t1<-summary(m_50)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_50)

m_51 <- brm(R2_Result ~ FZ_FEB+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_51, 
          type = "trace")

t1<-summary(m_51)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_51)


m_52 <- brm(R2_Result ~ AKT_c+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_52, 
          type = "trace")

t1<-summary(m_52)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_52)

m_53 <- brm(R2_Result ~ FEB_c+ Age+ Sex+R1_Result+ timeelapsed+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_53, 
          type = "trace")

t1<-summary(m_53)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_53)


m_54 <- brm(R2_Result ~risk_R2+risk_R1*R1_Result+ Age+ Sex+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_54, 
          type = "trace")

t1<-summary(m_54)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_54)
brm_pval(m_54)

m_55 <- brm(R2_Result ~risk_R2_c+risk_R1*R1_Result+
             Age+ Sex+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_55, 
          type = "trace")

t1<-summary(m_55)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_55)


m_56 <- brm(R2_Result ~Inf_risk_R2+Inf_risk_R1*R1_Result+
              Age+ Sex+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_56, 
          type = "trace")

t1<-summary(m_56)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_56)

m_57 <- brm(R2_Result ~Inf_risk_R2_c+Inf_risk_R1*R1_Result+
              Age+ Sex+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_57, 
          type = "trace")

t1<-summary(m_57)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_57)


m_58 <- brm(R2_Result ~Inf_grade_R2+Inf_grade_R1*R1_Result+
              Age+ Sex+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_58, 
          type = "trace")

t1<-summary(m_58)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_58)

m_59 <- brm(R2_Result ~Inf_grade_R2_c+Inf_grade_R1*R1_Result+
              Age+ Sex+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_59, 
          type = "trace")

t1<-summary(m_59)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_59)

m_60 <- brm(R2_Result ~sum_contact_R2+sum_contact_R1*R1_Result+
              Age+ Sex+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_60, 
          type = "trace")

t1<-summary(m_60)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_60)


m_61 <- brm(R2_Result ~FZ_AKT+FZ_FEB*R1_Result+
              Age+ Sex+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_61, 
          type = "trace")

t1<-summary(m_61)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_61)

m_62 <- brm(R2_Result ~AKT_c+FZ_FEB*R1_Result+
              Age+ Sex+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_62, 
          type = "trace")

t1<-summary(m_62)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_62)




m_54 <- brm(R2_Result ~risk_R2+risk_R1+R1_Result+ Age+ Sex+ (1|hh_id),
            data=koco_l,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_54, 
          type = "trace")

t1<-summary(m_54)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_54)
brm_pval(m_54)




# Ever Positive
koco_r2_dat$Agegroup <- factor(koco_r2_dat$Agegroup)
koco_r2_dat$Agegroup<- relevel(koco_r2_dat$Agegroup,ref="20-34")
koco_r2_dat$smokestatus<- factor(koco_r2_dat$smokestatus)
koco_r2_dat$smokestatus<- relevel(koco_r2_dat$smokestatus,ref="Non smoker")
koco_r2_dat$NetIncome_monthly_1<- factor(koco_r2_dat$NetIncome_monthly_1)
koco_r2_dat$NetIncome_monthly_1<- relevel(koco_r2_dat$NetIncome_monthly_1,ref="2500-4000")
koco_r2_dat$Housing_Type<- factor(koco_r2_dat$Housing_Type, levels = c("Buildings with 1-2 apartments",
                                                             "Buildings with 3-4 apartments",
                                                             "Buildings with >=5 apartments",
                                                             "Others"))
koco_r2_dat$HealthStatus_Self_Overall_1<- factor(koco_r2_dat$HealthStatus_Self_Overall_1,
                                            levels = c("excellent"   ,"very good"   ,"good","not good" ),
                                            labels=c("excellent","very good","good or lower","good or lower"))

koco_r2_dat$status<- factor(koco_r2_dat$status)
koco_r2_dat$status<- relevel(koco_r2_dat$status, ref="Never Positive")

koco_r2_dat<- merge.data.frame(koco_r2_dat,pos_hh[c(1,4)],by="hh_id",all.x = TRUE)
koco_r2_dat$hh_base<- factor(koco_r2_dat$hh_base, levels = c("Null","Not null"))


m_1 <- brm(status~Sex+ (1|hh_id),
           data=koco_r2_dat,
           chains = 3,warmup = 2000, 
           iter = 5000,
           control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_1, 
          type = "trace")

t1<-summary(m_1)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_1)


m_2 <- brm(status~Agegroup +Sex+ (1|hh_id),
           data=koco_r2_dat,
           chains = 3,warmup = 2000, 
           iter = 5000,
           control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_2, 
          type = "trace")

t1<-summary(m_2)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_2)


m_3 <- brm(status~Agegroup + (1|hh_id),
             data=koco_r2_dat,
             chains = 3,warmup = 2000, 
             iter = 5000,
             control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_3, 
          type = "trace")

t1<-summary(m_3)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_3)



m_3.1 <- brm(status~Age +Sex+ (1|hh_id),
           data=koco_r2_dat,
           chains = 3,warmup = 2000, 
           iter = 5000,
           control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_3.1, 
          type = "trace")

t1<-summary(m_3.1)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_3.1)



m_6 <- brm(status~ Birth_Country + Age+ Sex+(1|hh_id),
           data=koco_r2_dat,
           chains = 3,warmup = 2000, 
           iter = 5000,
           control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_6, 
          type = "trace")

t1<-summary(m_6)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_6)

m_7 <- brm(status ~ EmploymentStatus+Age+ Sex +(1|hh_id),
           data=koco_r2_dat,
           chains = 3,warmup = 2000, 
           iter = 5000,
           control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_7, 
          type = "trace")

t1<-summary(m_7)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_7)

m_8 <- brm(status ~ Employment_ind_lastyear+Age+ Sex+ (1|hh_id),
           data=koco_r2_dat,
           chains = 3,warmup = 2000, 
           iter = 5000,
           control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_8, 
          type = "trace")

t1<-summary(m_8)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_8)


m_9 <- brm(status ~ EducationLevel+Age+ Sex+ (1|hh_id),
           data=koco_r2_dat,
           chains = 3,warmup = 2000, 
           iter = 5000,
           control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_9, 
          type = "trace")

t1<-summary(m_9)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_9)


m_10 <- brm(status ~ smokestatus+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_10, 
          type = "trace")

t1<-summary(m_10)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_10)


m_11 <- brm(status ~ any_risk_emp+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_11, 
          type = "trace")

t1<-summary(m_11)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_11)


m_12 <- brm(status ~ NetIncome_monthly_1+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_12, 
          type = "trace")

t1<-summary(m_12)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_12)


m_13 <- brm(status ~ Household_Type_1+ Age+ Sex+(1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_13, 
          type = "trace")

t1<-summary(m_13)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_13)



m_14 <- brm(status~ Housing_Type+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat[which(koco_r2_dat$Housing_Type!="Others"),],
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_14, 
          type = "trace")

t1<-summary(m_14)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_14)


m_15 <- brm(status ~ Household_Inhabitants + Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_15, 
          type = "trace")

t1<-summary(m_15)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_15)


m_16 <- brm(status ~ LivingArea_perInhabitant_cat+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_16, 
          type = "trace")

t1<-summary(m_16)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_16)


m_17 <- brm(status ~ hh_base+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_17, 
          type = "trace")

t1<-summary(m_17)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_17)



m_18 <- brm(status ~Facemask_public+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_18, 
          type = "trace")

t1<-summary(m_18)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_18)




m_19 <- brm(status ~Symptoms_past2w_SenseLoss.Smell.Taste.+ Age+ Sex+(1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_19, 
          type = "trace")

t1<-summary(m_19)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_19)





m_20 <- brm(status ~Symptoms_past2w_Drycough+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_20, 
          type = "trace")

t1<-summary(m_20)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_20)


m_21 <- brm(status ~Symptoms_past2w_Sorethroat+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_21, 
          type = "trace")

t1<-summary(m_21)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_21)




m_22 <- brm(status ~Symptoms_past2w_Fever+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_22, 
          type = "trace")

t1<-summary(m_22)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_22)


m_23 <- brm(status ~Symptoms_past2w_Chills+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_23, 
          type = "trace")

t1<-summary(m_23)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_23)


m_24 <- brm(status ~Symptoms_past2w_Runningnose+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_24, 
          type = "trace")

t1<-summary(m_24)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_24)



m_25 <- brm(status ~num_symptoms_pos+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_25, 
          type = "trace")

t1<-summary(m_25)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_25)



m_26 <- brm(status ~any_symptoms+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_26, 
          type = "trace")

t1<-summary(m_26)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_26)


m_27 <- brm(status ~HealthStatus_Self_Overall_1+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_27, 
          type = "trace")

t1<-summary(m_27)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_27)


m_28 <- brm(status ~HealthStatus_Self_Allergies_Respiratory+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_28, 
          type = "trace")

t1<-summary(m_28)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_28)


m_29 <- brm(status ~HealthStatus_Self_Allergies_Skin+ Age+ Sex+(1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_29, 
          type = "trace")

t1<-summary(m_29)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_29)


m_30 <- brm(status ~ HealthStatus_Self_Autoimmune + Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_30, 
          type = "trace")

t1<-summary(m_30)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_30)


m_31 <- brm(status ~ HealthStatus_Self_Cancer + Age+ Sex+(1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_31, 
          type = "trace")

t1<-summary(m_31)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_31)


m_32 <- brm(status ~ HealthStatus_Self_CardiovascularDisease+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_32, 
          type = "trace")

t1<-summary(m_32)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_32)


m_33 <- brm(status ~ HealthStatus_Self_Diabetes+ Age+ Sex+(1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_33, 
          type = "trace")

t1<-summary(m_33)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_33)


m_34 <- brm(status ~ HealthStatus_Self_LungDisease+ Age+ Sex+(1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_34, 
          type = "trace")

t1<-summary(m_34)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_34)

m_35 <- brm(status ~ HealthStatus_Self_Obesity+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_35, 
          type = "trace")

t1<-summary(m_35)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_35)


m_36 <- brm(status ~ risk_R1+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_36, 
          type = "trace")

t1<-summary(m_36)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_36)

m_37 <- brm(status ~risk_R1_c+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_37, 
          type = "trace")

t1<-summary(m_37)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_37)


m_38 <- brm(status ~Inf_risk_R1+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_38, 
          type = "trace")

t1<-summary(m_38)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_38)

m_39 <- brm(status ~Inf_risk_R1_c+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_39, 
          type = "trace")

t1<-summary(m_39)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_39)


m_40 <- brm(status ~Inf_grade_R1+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_40, 
          type = "trace")

t1<-summary(m_40)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_40)

m_41 <- brm(status ~Inf_grade_R1_c+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_41, 
          type = "trace")

t1<-summary(m_41)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_41)


m_42 <- brm(status ~risk_R2+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_42, 
          type = "trace")

t1<-summary(m_42)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_42)
brm_pval(m_42)

m_43 <- brm(status ~risk_R2_c+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_43, 
          type = "trace")

t1<-summary(m_43)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_43)


m_44 <- brm(status ~Inf_risk_R2+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_44, 
          type = "trace")

t1<-summary(m_44)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_44)

m_45 <- brm(status ~Inf_risk_R2_c+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_45, 
          type = "trace")

t1<-summary(m_45)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_45)


m_46 <- brm(status ~Inf_grade_R2+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_46, 
          type = "trace")

t1<-summary(m_46)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_46)

m_47 <- brm(status ~Inf_grade_R2_c+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_47, 
          type = "trace")

t1<-summary(m_47)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_47)


m_48 <- brm(status ~sum_contact_R1+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_48, 
          type = "trace")

t1<-summary(m_48)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_48)

m_49 <- brm(status ~sum_contact_R2+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_49, 
          type = "trace")

t1<-summary(m_49)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_49)


m_50 <- brm(status ~FZ_AKT+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_50, 
          type = "trace")

t1<-summary(m_50)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_50)

m_51 <- brm(status ~FZ_FEB+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_51, 
          type = "trace")

t1<-summary(m_51)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_51)



m_52 <- brm(status ~AKT_c+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_52, 
          type = "trace")

t1<-summary(m_52)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_52)

m_53 <- brm(status ~FEB_c+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_53, 
          type = "trace")

t1<-summary(m_53)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_53)


m_54 <- brm(status ~risk_R2+risk_R1+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_54, 
          type = "trace")

t1<-summary(m_54)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_54)
brm_pval(m_54)

m_55 <- brm(status ~risk_R2_c+risk_R1_c+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_55, 
          type = "trace")

t1<-summary(m_55)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_55)


m_56 <- brm(status ~Inf_risk_R2+Inf_risk_R1+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_56, 
          type = "trace")

t1<-summary(m_56)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_56)

m_57 <- brm(status ~Inf_risk_R2_c+ Inf_risk_R1_c+Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_57, 
          type = "trace")

t1<-summary(m_57)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_57)


m_58 <- brm(status ~Inf_grade_R2+Inf_grade_R1+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_58, 
          type = "trace")

t1<-summary(m_58)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_58)

m_59 <- brm(status ~Inf_grade_R2_c+Inf_grade_R1_c+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_59, 
          type = "trace")

t1<-summary(m_59)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_59)

m_60 <- brm(status ~sum_contact_R2+sum_contact_R1+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_60, 
          type = "trace")

t1<-summary(m_60)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_60)


m_61 <- brm(status ~FZ_AKT+FZ_FEB+ Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_61, 
          type = "trace")

t1<-summary(m_61)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_61)

m_62 <- brm(status ~AKT_c+FEB_c +Age+ Sex+ (1|hh_id),
            data=koco_r2_dat,
            chains = 3,warmup = 2000, 
            iter = 5000,
            control = list(adapt_delta = 0.95),family = bernoulli(link = "logit"),seed=1234)
mcmc_plot(m_62, 
          type = "trace")

t1<-summary(m_62)

exp(t1$fixed[,c(1,3,4)])

p_significance(m_62)


##############################################################################################
# Imputation Round 2
koco_r2_dat$Sex<- factor(koco_r2_dat$Sex)
m_1<-glmer_imp(status ~Sex,
               data=koco_r2_dat,
               random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
               thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_1)
summary(m_1)
cc_m_1<- cbind(summary(m_1)$res$status$regcoef[, c('Mean')]
               ,summary(m_1)$res$status$regcoef[, c('2.5%')],
               summary(m_1)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_1)

m_2<-glmer_imp(status~Agegroup +Sex,
               data=koco_r2_dat,
               random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
               thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_2)
summary(m_2)

cc_m_2<- cbind(summary(m_2)$res$status$regcoef[, c('Mean')]
               ,summary(m_2)$res$status$regcoef[, c('2.5%')],
               summary(m_2)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_2)

m_3<-glmer_imp(status~Agegroup,
               data=koco_r2_dat,
               random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
               thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_3)
summary(m_3)
cc_m_3<- cbind(summary(m_3)$res$status$regcoef[, c('Mean')]
               ,summary(m_3)$res$status$regcoef[, c('2.5%')],
               summary(m_3)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_3)

m_4<-glmer_imp(status~Age +Sex,
               data=koco_r2_dat,
               random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
               thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_4)
summary(m_4)
cc_m_4<- cbind(summary(m_4)$res$status$regcoef[, c('Mean')]
               ,summary(m_4)$res$status$regcoef[, c('2.5%')],
               summary(m_4)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_4)



koco_r2_dat$Birth_Country<- factor(koco_r2_dat$Birth_Country)
m_6<-glmer_imp(status~ Birth_Country + Age+ Sex,
               data=koco_r2_dat,
               random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
               thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_6)
summary(m_6)

cc_m_6<- cbind(summary(m_6)$res$status$regcoef[, c('Mean')]
               ,summary(m_6)$res$status$regcoef[, c('2.5%')],
               summary(m_6)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_6)

koco_r2_dat$EmploymentStatus<- factor(koco_r2_dat$EmploymentStatus)
m_7<-glmer_imp(status ~ EmploymentStatus+Age+ Sex ,
               data=koco_r2_dat,
               random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
               thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_7)
summary(m_7)

cc_m_7<- cbind(summary(m_7)$res$status$regcoef[, c('Mean')]
               ,summary(m_7)$res$status$regcoef[, c('2.5%')],
               summary(m_7)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_7)

koco_r2_dat$Employment_ind_lastyear<- factor(koco_r2_dat$Employment_ind_lastyear)
m_8<-glmer_imp(status ~ Employment_ind_lastyear+Age+ Sex ,
               data=koco_r2_dat,
               random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
               thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_8)
summary(m_8)

cc_m_8<- cbind(summary(m_8)$res$status$regcoef[, c('Mean')]
               ,summary(m_8)$res$status$regcoef[, c('2.5%')],
               summary(m_8)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_8)

koco_r2_dat$EducationLevel<-factor(koco_r2_dat$EducationLevel)
m_9<-glmer_imp(status ~ EducationLevel+Age+ Sex ,
               data=koco_r2_dat,
               random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
               thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_9)
summary(m_9)
cc_m_9<- cbind(summary(m_9)$res$status$regcoef[, c('Mean')]
               ,summary(m_9)$res$status$regcoef[, c('2.5%')],
               summary(m_9)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_9)


koco_r2_dat$smokestatus<- factor(koco_r2_dat$smokestatus)
m_10<-glmer_imp(status ~ smokestatus+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_10)
summary(m_10)

cc_m_10<- cbind(summary(m_10)$res$status$regcoef[, c('Mean')]
                ,summary(m_10)$res$status$regcoef[, c('2.5%')],
                summary(m_10)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_10)


koco_r2_dat$any_risk_emp<- factor(koco_r2_dat$any_risk_emp)
m_11<-glmer_imp(status ~ any_risk_emp+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_11)
summary(m_11)

cc_m_11<- cbind(summary(m_11)$res$status$regcoef[, c('Mean')]
                ,summary(m_11)$res$status$regcoef[, c('2.5%')],
                summary(m_11)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_11)

m_12<-glmer_imp(status ~ NetIncome_monthly_1+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))

traceplot(m_12)
summary(m_12)

cc_m_12<- cbind(summary(m_12)$res$status$regcoef[, c('Mean')]
                ,summary(m_12)$res$status$regcoef[, c('2.5%')],
                summary(m_12)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_12)

koco_r2_dat$Household_Type_1<- factor(koco_r2_dat$Household_Type_1)
m_13<-glmer_imp(status ~ Household_Type_1+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))

traceplot(m_13)
summary(m_13)

cc_m_13<- cbind(summary(m_13)$res$status$regcoef[, c('Mean')]
                ,summary(m_13)$res$status$regcoef[, c('2.5%')],
                summary(m_13)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_13)

koco_r2_dat$Housing_Type<- factor(koco_r2_dat$Housing_Type)
m_14<-glmer_imp(status ~ Housing_Type+Age+ Sex ,
                data=koco_r2_dat[which(koco_r2_dat$Housing_Type!="Others"),],
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))

traceplot(m_14)
summary(m_14)

cc_m_14<- cbind(summary(m_14)$res$status$regcoef[, c('Mean')]
                ,summary(m_14)$res$status$regcoef[, c('2.5%')],
                summary(m_14)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_14)


koco_r2_dat$Household_Inhabitants<- factor(koco_r2_dat$Household_Inhabitants)
m_15<-glmer_imp(status ~ Household_Inhabitants+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_15)
summary(m_15)

cc_m_15<- cbind(summary(m_15)$res$status$regcoef[, c('Mean')]
                ,summary(m_15)$res$status$regcoef[, c('2.5%')],
                summary(m_15)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_15)


koco_r2_dat$LivingArea_perInhabitant_cat<- factor(koco_r2_dat$LivingArea_perInhabitant_cat)
m_16<-glmer_imp(status ~ LivingArea_perInhabitant_cat+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_16)
summary(m_16)

cc_m_16<- cbind(summary(m_16)$res$status$regcoef[, c('Mean')]
                ,summary(m_16)$res$status$regcoef[, c('2.5%')],
                summary(m_16)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_16)




koco_r2_dat$Facemask_public<- factor(koco_r2_dat$Facemask_public)
m_18<-glmer_imp(status ~ Facemask_public+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_18)
summary(m_18)

cc_m_18<- cbind(summary(m_18)$res$status$regcoef[, c('Mean')]
                ,summary(m_18)$res$status$regcoef[, c('2.5%')],
                summary(m_18)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_18)

koco_r2_dat$Symptoms_past2w_SenseLoss.Smell.Taste.<- factor(koco_r2_dat$Symptoms_past2w_SenseLoss.Smell.Taste.)
m_19<-glmer_imp(status ~Symptoms_past2w_SenseLoss.Smell.Taste.+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_19)
summary(m_19)

cc_m_19<- cbind(summary(m_19)$res$status$regcoef[, c('Mean')]
                ,summary(m_19)$res$status$regcoef[, c('2.5%')],
                summary(m_19)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_19)



koco_r2_dat$Symptoms_past2w_Drycough<- factor(koco_r2_dat$Symptoms_past2w_Drycough)
m_20<-glmer_imp(status ~Symptoms_past2w_Drycough+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_20)
summary(m_20)

cc_m_20<- cbind(summary(m_20)$res$status$regcoef[, c('Mean')]
                ,summary(m_20)$res$status$regcoef[, c('2.5%')],
                summary(m_20)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_20)



koco_r2_dat$Symptoms_past2w_Sorethroat<- factor(koco_r2_dat$Symptoms_past2w_Sorethroat)
m_21<-glmer_imp(status ~Symptoms_past2w_Sorethroat+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_21)
summary(m_21)

cc_m_21<- cbind(summary(m_21)$res$status$regcoef[, c('Mean')]
                ,summary(m_21)$res$status$regcoef[, c('2.5%')],
                summary(m_21)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_21)


koco_r2_dat$Symptoms_past2w_Fever<- factor(koco_r2_dat$Symptoms_past2w_Fever)
m_22<-glmer_imp(status ~Symptoms_past2w_Fever+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_22)
summary(m_22)

cc_m_22<- cbind(summary(m_22)$res$status$regcoef[, c('Mean')]
                ,summary(m_22)$res$status$regcoef[, c('2.5%')],
                summary(m_22)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_22)


koco_r2_dat$Symptoms_past2w_Chills<- factor(koco_r2_dat$Symptoms_past2w_Chills)
m_23<-glmer_imp(status ~Symptoms_past2w_Chills+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_23)
summary(m_23)

cc_m_23<- cbind(summary(m_23)$res$status$regcoef[, c('Mean')]
                ,summary(m_23)$res$status$regcoef[, c('2.5%')],
                summary(m_23)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_23)

koco_r2_dat$Symptoms_past2w_Runningnose<- factor(koco_r2_dat$Symptoms_past2w_Runningnose)
m_24<-glmer_imp(status ~Symptoms_past2w_Runningnose+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_24)
summary(m_24)

cc_m_24<- cbind(summary(m_24)$res$status$regcoef[, c('Mean')]
                ,summary(m_24)$res$status$regcoef[, c('2.5%')],
                summary(m_24)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_24)


m_25<-glmer_imp(status ~num_symptoms_pos+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_25)
summary(m_25)

cc_m_25<- cbind(summary(m_25)$res$status$regcoef[, c('Mean')]
                ,summary(m_25)$res$status$regcoef[, c('2.5%')],
                summary(m_25)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_25)

koco_r2_dat$any_symptoms<- factor(koco_r2_dat$any_symptoms)
m_26<-glmer_imp(status ~any_symptoms+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_26)
summary(m_26)

cc_m_26<- cbind(summary(m_26)$res$status$regcoef[, c('Mean')]
                ,summary(m_26)$res$status$regcoef[, c('2.5%')],
                summary(m_26)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_26)


m_27<-glmer_imp(status ~HealthStatus_Self_Overall_1+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_27)
summary(m_27)

cc_m_27<- cbind(summary(m_27)$res$status$regcoef[, c('Mean')]
                ,summary(m_27)$res$status$regcoef[, c('2.5%')],
                summary(m_27)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_27)

koco_r2_dat$HealthStatus_Self_Allergies_Respiratory<- factor(koco_r2_dat$HealthStatus_Self_Allergies_Respiratory)
m_28<-glmer_imp(status ~HealthStatus_Self_Allergies_Respiratory+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_28)
summary(m_28)

cc_m_28<- cbind(summary(m_28)$res$status$regcoef[, c('Mean')]
                ,summary(m_28)$res$status$regcoef[, c('2.5%')],
                summary(m_28)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_28)


koco_r2_dat$HealthStatus_Self_Allergies_Skin<- factor(koco_r2_dat$HealthStatus_Self_Allergies_Skin)
m_29<-glmer_imp(status ~HealthStatus_Self_Allergies_Skin+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_29)
summary(m_29)

cc_m_29<- cbind(summary(m_29)$res$status$regcoef[, c('Mean')]
                ,summary(m_29)$res$status$regcoef[, c('2.5%')],
                summary(m_29)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_29)

koco_r2_dat$HealthStatus_Self_Autoimmune<- factor(koco_r2_dat$HealthStatus_Self_Autoimmune)
m_30<-glmer_imp(status ~ HealthStatus_Self_Autoimmune+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_30)
summary(m_30)

cc_m_30<- cbind(summary(m_30)$res$status$regcoef[, c('Mean')]
                ,summary(m_30)$res$status$regcoef[, c('2.5%')],
                summary(m_30)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_30)

koco_r2_dat$HealthStatus_Self_Cancer<- factor(koco_r2_dat$HealthStatus_Self_Cancer)
m_31<-glmer_imp(status ~ HealthStatus_Self_Cancer+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_31)
summary(m_31)

cc_m_31<- cbind(summary(m_31)$res$status$regcoef[, c('Mean')]
                ,summary(m_31)$res$status$regcoef[, c('2.5%')],
                summary(m_31)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_31)

koco_r2_dat$HealthStatus_Self_CardiovascularDisease<- factor(koco_r2_dat$HealthStatus_Self_CardiovascularDisease)
m_32<-glmer_imp(status ~ HealthStatus_Self_CardiovascularDisease+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_32)
summary(m_32)

cc_m_32<- cbind(summary(m_32)$res$status$regcoef[, c('Mean')]
                ,summary(m_32)$res$status$regcoef[, c('2.5%')],
                summary(m_32)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_32)

koco_r2_dat$HealthStatus_Self_Diabetes<- factor(koco_r2_dat$HealthStatus_Self_Diabetes)
m_33<-glmer_imp(status ~ HealthStatus_Self_Diabetes+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_33)
summary(m_33)

cc_m_33<- cbind(summary(m_33)$res$status$regcoef[, c('Mean')]
                ,summary(m_33)$res$status$regcoef[, c('2.5%')],
                summary(m_33)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_33)


koco_r2_dat$HealthStatus_Self_LungDisease<- factor(koco_r2_dat$HealthStatus_Self_LungDisease)
m_34<-glmer_imp(status ~ HealthStatus_Self_LungDisease+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_34)
summary(m_34)

cc_m_34<- cbind(summary(m_34)$res$status$regcoef[, c('Mean')]
                ,summary(m_34)$res$status$regcoef[, c('2.5%')],
                summary(m_34)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_34)

koco_r2_dat$HealthStatus_Self_Obesity<- factor(koco_r2_dat$HealthStatus_Self_Obesity)
m_35<-glmer_imp(status ~ HealthStatus_Self_Obesity+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_35)
summary(m_35)

cc_m_35<- cbind(summary(m_35)$res$status$regcoef[, c('Mean')]
                ,summary(m_35)$res$status$regcoef[, c('2.5%')],
                summary(m_35)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_35)

m_36<-glmer_imp(status ~ risk_R1+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_36)
summary(m_36)

cc_m_36<- cbind(summary(m_36)$res$status$regcoef[, c('Mean')]
                ,summary(m_36)$res$status$regcoef[, c('2.5%')],
                summary(m_36)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_36)

m_37<-glmer_imp(status ~ risk_R1_c+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_37)
summary(m_37)

cc_m_37<- cbind(summary(m_37)$res$status$regcoef[, c('Mean')]
                ,summary(m_37)$res$status$regcoef[, c('2.5%')],
                summary(m_37)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_37)


m_38<-glmer_imp(status ~ Inf_risk_R1+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_38)
summary(m_38)

cc_m_38<- cbind(summary(m_38)$res$status$regcoef[, c('Mean')]
                ,summary(m_38)$res$status$regcoef[, c('2.5%')],
                summary(m_38)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_38)

m_39<-glmer_imp(status ~ Inf_risk_R1_c+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_39)
summary(m_39)

cc_m_39<- cbind(summary(m_39)$res$status$regcoef[, c('Mean')]
                ,summary(m_39)$res$status$regcoef[, c('2.5%')],
                summary(m_39)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_39)


m_40<-glmer_imp(status ~ Inf_grade_R1+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_40)
summary(m_40)

cc_m_40<- cbind(summary(m_40)$res$status$regcoef[, c('Mean')]
                ,summary(m_40)$res$status$regcoef[, c('2.5%')],
                summary(m_40)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_40)

m_41<-glmer_imp(status ~ Inf_grade_R1_c+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_41)
summary(m_41)

cc_m_41<- cbind(summary(m_41)$res$status$regcoef[, c('Mean')]
                ,summary(m_41)$res$status$regcoef[, c('2.5%')],
                summary(m_41)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_41)



m_42<-glmer_imp(status ~ risk_R2+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_42)
summary(m_42)

cc_m_42<- cbind(summary(m_42)$res$status$regcoef[, c('Mean')]
                ,summary(m_42)$res$status$regcoef[, c('2.5%')],
                summary(m_42)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_42)

m_43<-glmer_imp(status ~ risk_R2_c+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_43)
summary(m_43)

cc_m_43<- cbind(summary(m_43)$res$status$regcoef[, c('Mean')]
                ,summary(m_43)$res$status$regcoef[, c('2.5%')],
                summary(m_43)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_43)


m_44<-glmer_imp(status ~ Inf_risk_R2+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_44)
summary(m_44)

cc_m_44<- cbind(summary(m_44)$res$status$regcoef[, c('Mean')]
                ,summary(m_44)$res$status$regcoef[, c('2.5%')],
                summary(m_44)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_44)

m_45<-glmer_imp(status ~ Inf_risk_R2_c+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_45)
summary(m_45)

cc_m_45<- cbind(summary(m_45)$res$status$regcoef[, c('Mean')]
                ,summary(m_45)$res$status$regcoef[, c('2.5%')],
                summary(m_45)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_45)


m_46<-glmer_imp(status ~ Inf_grade_R2+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_46)
summary(m_46)

cc_m_46<- cbind(summary(m_46)$res$status$regcoef[, c('Mean')]
                ,summary(m_46)$res$status$regcoef[, c('2.5%')],
                summary(m_46)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_46)

m_47<-glmer_imp(status ~ Inf_grade_R2_c+Age+ Sex ,
                data=koco_r2_dat,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_47)
summary(m_47)

cc_m_47<- cbind(summary(m_47)$res$status$regcoef[, c('Mean')]
                ,summary(m_47)$res$status$regcoef[, c('2.5%')],
                summary(m_47)$res$status$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_47)



##############################################################################################
# Imputation
koco_l$Sex<- factor(koco_l$Sex)
koco_l$R1_Result<- factor(koco_l$R1_Result)
m_1<-glmer_imp(R2_Result~ Sex+R1_Result+ timeelapsed ,
               data=koco_l,
                 random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                 thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_1)
summary(m_1)
cc_m_1<- cbind(summary(m_1)$res$R2_Result$regcoef[, c('Mean')]
                       ,summary(m_1)$res$R2_Result$regcoef[, c('2.5%')],
                       summary(m_1)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_1)

m_2<-glmer_imp(R2_Result~ Agegroup+ Sex+R1_Result+ timeelapsed ,
               data=koco_l,
              random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
              thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_2)
summary(m_2)

cc_m_2<- cbind(summary(m_2)$res$R2_Result$regcoef[, c('Mean')]
               ,summary(m_2)$res$R2_Result$regcoef[, c('2.5%')],
               summary(m_2)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_2)

m_3<-glmer_imp(R2_Result~ Agegroup +R1_Result+ timeelapsed ,
               data=koco_l,
               random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
               thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_3)
summary(m_3)
cc_m_3<- cbind(summary(m_3)$res$R2_Result$regcoef[, c('Mean')]
               ,summary(m_3)$res$R2_Result$regcoef[, c('2.5%')],
               summary(m_3)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_3)

m_4<-glmer_imp(
               R2_Result~ Age+ Sex+R1_Result+ timeelapsed ,
               data=koco_l,
               random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
               thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_4)
summary(m_4)
cc_m_4<- cbind(summary(m_4)$res$R2_Result$regcoef[, c('Mean')]
               ,summary(m_4)$res$R2_Result$regcoef[, c('2.5%')],
               summary(m_4)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_4)



koco_l$Birth_Country<- factor(koco_l$Birth_Country)
m_6<-glmer_imp(R2_Result~ Birth_Country+ Sex+R1_Result+ timeelapsed ,
               data=koco_l,
               random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
               thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_6)
summary(m_6)

cc_m_6<- cbind(summary(m_6)$res$R2_Result$regcoef[, c('Mean')]
               ,summary(m_6)$res$R2_Result$regcoef[, c('2.5%')],
               summary(m_6)$res$R2_Result$regcoef[, c('97.5%')]) ## slow (~ 11 seconds)

exp(cc_m_6)

koco_l$EmploymentStatus<- factor(koco_l$EmploymentStatus)
m_7<-glmer_imp(R2_Result~ EmploymentStatus+ Sex+R1_Result+ timeelapsed ,
               data=koco_l,
               random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
               thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_7)
summary(m_7)

cc_m_7<- cbind(summary(m_7)$res$R2_Result$regcoef[, c('Mean')]
               ,summary(m_7)$res$R2_Result$regcoef[, c('2.5%')],
               summary(m_7)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_7)

koco_l$Employment_ind_lastyear<- factor(koco_l$Employment_ind_lastyear)
m_8<-glmer_imp(R2_Result~ Employment_ind_lastyear+ Sex+R1_Result+ timeelapsed ,
               data=koco_l,
               random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
               thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_8)
summary(m_8)

cc_m_8<- cbind(summary(m_8)$res$R2_Result$regcoef[, c('Mean')]
               ,summary(m_8)$res$R2_Result$regcoef[, c('2.5%')],
               summary(m_8)$res$R2_Result$regcoef[, c('97.5%')]) ## slow (~ 11 seconds)

exp(cc_m_8)

koco_l$EducationLevel<-factor(koco_l$EducationLevel)
m_9<-glmer_imp(R2_Result~ EducationLevel+ Sex+R1_Result+ timeelapsed ,
               data=koco_l,
               random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
               thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_9)
summary(m_9)
cc_m_9<- cbind(summary(m_9)$res$R2_Result$regcoef[, c('Mean')]
               ,summary(m_9)$res$R2_Result$regcoef[, c('2.5%')],
               summary(m_9)$res$R2_Result$regcoef[, c('97.5%')]) ## slow (~ 11 seconds)

exp(cc_m_9)


koco_l$smokestatus<- factor(koco_l$smokestatus)
m_10<-glmer_imp(R2_Result~ smokestatus+ Sex+R1_Result+ timeelapsed ,
                data=koco_l,
               random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
               thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_10)
summary(m_10)

cc_m_10<- cbind(summary(m_10)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_10)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_10)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_10)


koco_l$any_risk_emp<- factor(koco_l$any_risk_emp)
m_11<-glmer_imp(R2_Result~ any_risk_emp+ Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_11)
summary(m_11)

cc_m_11<- cbind(summary(m_11)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_11)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_11)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_11)

m_12<-glmer_imp(
                R2_Result~ NetIncome_monthly_1+ Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))

traceplot(m_12)
summary(m_12)

cc_m_12<- cbind(summary(m_12)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_12)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_12)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_12)

koco_l$Household_Type_1<- factor(koco_l$Household_Type_1)
m_13<-glmer_imp(R2_Result~ Household_Type_1+ Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))

traceplot(m_13)
summary(m_13)

cc_m_13<- cbind(summary(m_13)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_13)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_13)$res$R2_Result$regcoef[, c('97.5%')])   ## slow (~ 11 seconds)

exp(cc_m_13)

koco_l$Housing_Type<- factor(koco_l$Housing_Type)
m_14<-glmer_imp(R2_Result~ Housing_Type + Sex+R1_Result+ timeelapsed ,
                data=koco_l[which(koco_l$Housing_Type!="Others"),],
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))

traceplot(m_14)
summary(m_14)

cc_m_14<- cbind(summary(m_14)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_14)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_14)$res$R2_Result$regcoef[, c('97.5%')])   ## slow (~ 11 seconds)

exp(cc_m_14)


koco_l$Household_Inhabitants<- factor(koco_l$Household_Inhabitants)
m_15<-glmer_imp(R2_Result~Household_Inhabitants + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_15)
summary(m_15)

cc_m_15<- cbind(summary(m_15)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_15)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_15)$res$R2_Result$regcoef[, c('97.5%')])   ## slow (~ 11 seconds)

exp(cc_m_15)


koco_l$LivingArea_perInhabitant_cat<- factor(koco_l$LivingArea_perInhabitant_cat)
m_16<-glmer_imp(R2_Result~LivingArea_perInhabitant_cat + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_16)
summary(m_16)

cc_m_16<- cbind(summary(m_16)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_16)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_16)$res$R2_Result$regcoef[, c('97.5%')])   ## slow (~ 11 seconds)

exp(cc_m_16)




koco_l$Facemask_public<- factor(koco_l$Facemask_public)
m_18<-glmer_imp(R2_Result~Facemask_public + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_18)
summary(m_18)

cc_m_18<- cbind(summary(m_18)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_18)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_18)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_18)

koco_l$Symptoms_past2w_SenseLoss.Smell.Taste.<- factor(koco_l$Symptoms_past2w_SenseLoss.Smell.Taste.)
m_19<-glmer_imp(R2_Result~Symptoms_past2w_SenseLoss.Smell.Taste. + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_19)
summary(m_19)

cc_m_19<- cbind(summary(m_19)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_19)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_19)$res$R2_Result$regcoef[, c('97.5%')])   ## slow (~ 11 seconds)

exp(cc_m_19)



koco_l$Symptoms_past2w_Drycough<- factor(koco_l$Symptoms_past2w_Drycough)
m_20<-glmer_imp(R2_Result~Symptoms_past2w_Drycough + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_20)
summary(m_20)

cc_m_20<- cbind(summary(m_20)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_20)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_20)$res$R2_Result$regcoef[, c('97.5%')])   ## slow (~ 11 seconds)

exp(cc_m_20)



koco_l$Symptoms_past2w_Sorethroat<- factor(koco_l$Symptoms_past2w_Sorethroat)
m_21<-glmer_imp(R2_Result~Symptoms_past2w_Sorethroat + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_21)
summary(m_21)

cc_m_21<- cbind(summary(m_21)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_21)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_21)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_21)


koco_l$Symptoms_past2w_Fever<- factor(koco_l$Symptoms_past2w_Fever)
m_22<-glmer_imp(R2_Result~Symptoms_past2w_Fever + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_22)
summary(m_22)

cc_m_22<- cbind(summary(m_22)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_22)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_22)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_22)


koco_l$Symptoms_past2w_Chills<- factor(koco_l$Symptoms_past2w_Chills)
m_23<-glmer_imp(R2_Result~Symptoms_past2w_Chills + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_23)
summary(m_23)

cc_m_23<- cbind(summary(m_23)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_23)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_23)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_23)

koco_l$Symptoms_past2w_Runningnose<- factor(koco_l$Symptoms_past2w_Runningnose)
m_24<-glmer_imp(R2_Result~Symptoms_past2w_Runningnose + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_24)
summary(m_24)

cc_m_24<- cbind(summary(m_24)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_24)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_24)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_24)


m_25<-glmer_imp(R2_Result~num_symptoms_pos + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_25)
summary(m_25)

cc_m_25<- cbind(summary(m_25)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_25)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_25)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_25)

koco_l$any_symptoms<- factor(koco_l$any_symptoms)
m_26<-glmer_imp(R2_Result~any_symptoms+Age + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_26)
summary(m_26)

cc_m_26<- cbind(summary(m_26)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_26)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_26)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_26)


m_27<-glmer_imp(R2_Result~HealthStatus_Self_Overall_1+Age + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_27)
summary(m_27)

cc_m_27<- cbind(summary(m_27)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_27)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_27)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_27)

koco_l$HealthStatus_Self_Allergies_Respiratory<- factor(koco_l$HealthStatus_Self_Allergies_Respiratory)
m_28<-glmer_imp(R2_Result~HealthStatus_Self_Allergies_Respiratory+Age + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_28)
summary(m_28)

cc_m_28<- cbind(summary(m_28)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_28)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_28)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_28)


koco_l$HealthStatus_Self_Allergies_Skin<- factor(koco_l$HealthStatus_Self_Allergies_Skin)
m_29<-glmer_imp(R2_Result~HealthStatus_Self_Allergies_Skin+Age + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_29)
summary(m_29)

cc_m_29<- cbind(summary(m_29)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_29)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_29)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_29)

koco_l$HealthStatus_Self_Autoimmune<- factor(koco_l$HealthStatus_Self_Autoimmune)
m_30<-glmer_imp(R2_Result~HealthStatus_Self_Autoimmune+Age + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_30)
summary(m_30)

cc_m_30<- cbind(summary(m_30)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_30)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_30)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_30)

koco_l$HealthStatus_Self_Cancer<- factor(koco_l$HealthStatus_Self_Cancer)
m_31<-glmer_imp(R2_Result~HealthStatus_Self_Cancer+Age + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_31)
summary(m_31)

cc_m_31<- cbind(summary(m_31)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_31)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_31)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_31)

koco_l$HealthStatus_Self_CardiovascularDisease<- factor(koco_l$HealthStatus_Self_CardiovascularDisease)
m_32<-glmer_imp(R2_Result~HealthStatus_Self_CardiovascularDisease+Age + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_32)
summary(m_32)

cc_m_32<- cbind(summary(m_32)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_32)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_32)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_32)

koco_l$HealthStatus_Self_Diabetes<- factor(koco_l$HealthStatus_Self_Diabetes)
m_33<-glmer_imp(R2_Result~HealthStatus_Self_Diabetes+Age + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_33)
summary(m_33)

cc_m_33<- cbind(summary(m_33)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_33)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_33)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_33)


koco_l$HealthStatus_Self_LungDisease<- factor(koco_l$HealthStatus_Self_LungDisease)
m_34<-glmer_imp(R2_Result~HealthStatus_Self_LungDisease+Age + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_34)
summary(m_34)

cc_m_34<- cbind(summary(m_34)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_34)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_34)$res$R2_Result$regcoef[, c('97.5%')]) ## slow (~ 11 seconds)

exp(cc_m_34)

koco_l$HealthStatus_Self_Obesity<- factor(koco_l$HealthStatus_Self_Obesity)
m_35<-glmer_imp(R2_Result~HealthStatus_Self_Obesity+Age + Sex+R1_Result+ timeelapsed ,
                data=koco_l,
                random=~(1|hh_id), n.chains = 3, n.adapt = 2000, n.iter = 3000,
                thin = 1,seed=1234,family =binomial(link = "logit"))
traceplot(m_35)
summary(m_35)

cc_m_35<- cbind(summary(m_35)$res$R2_Result$regcoef[, c('Mean')]
                ,summary(m_35)$res$R2_Result$regcoef[, c('2.5%')],
                summary(m_35)$res$R2_Result$regcoef[, c('97.5%')])  ## slow (~ 11 seconds)

exp(cc_m_35)


# Tables for behavioural data
bdat_l<- subset(koco_l, 
                select = c(ind_id,Sex,Age,R2_Result,
                           risk_R1,risk_R1_c,                              
                           Inf_risk_R1,Inf_risk_R1_c,                          
                           Inf_grade_R1,Inf_grade_R1_c,                         
                           risk_R2,risk_R2_c,                              
                           Inf_risk_R2,Inf_risk_R2_c,                          
                           Inf_grade_R2,Inf_grade_R2_c,                         
                           sum_contact_R1,sum_contact_R2,                         
                           FZ_AKT,FZ_FEB,                                 
                           AKT_c,FEB_c))
bdat_l$outcome<- paste0(bdat_l$R2_Result,bdat_l$Sex)

bdat_l$Agegroup<-ifelse(bdat_l$Age<35,"<35",
                        ifelse(bdat_l$Age>=65,">=65","35-65"))

bdat_l$Agegroup<- factor(bdat_l$Agegroup,levels = c("<35","35-65",">=65"))
bdat_l$outcome<- paste0(bdat_l$R2_Result,bdat_l$Sex)
bdat_l$outcome1<- paste0(bdat_l$R2_Result,bdat_l$Agegroup)

bdat_l$outcome1<- factor(bdat_l$outcome1,
                         levels = c("Negative<35","Negative35-65","Negative>=65",
                                    "Positive<35","Positive35-65","Positive>=65"))

t7<-tableby(outcome~ 
              risk_R1+chisq(risk_R1_c)+
              Inf_risk_R1+ chisq(Inf_risk_R1_c)+
              Inf_grade_R1+chisq(Inf_grade_R1_c)+
              risk_R2+chisq(risk_R2_c)+
              Inf_risk_R2+chisq(Inf_risk_R2_c)+
              Inf_grade_R2+chisq(Inf_grade_R2_c)+
              sum_contact_R1+chisq(sum_contact_R2)+
              FZ_AKT+FZ_FEB+
              chisq(AKT_c)+chisq(FEB_c)
            ,data =bdat_l[which(bdat_l$Sex=="Male"),],
            control = my_controls , simulate.p.value = TRUE , B=10000)

write2html(summary(t7, pfootnote=TRUE),"BehaviourPosR2_Full_Male.html")

t7<-tableby(outcome~ 
              risk_R1+chisq(risk_R1_c)+
              Inf_risk_R1+ chisq(Inf_risk_R1_c)+
              Inf_grade_R1+chisq(Inf_grade_R1_c)+
              risk_R2+chisq(risk_R2_c)+
              Inf_risk_R2+chisq(Inf_risk_R2_c)+
              Inf_grade_R2+chisq(Inf_grade_R2_c)+
              sum_contact_R1+chisq(sum_contact_R2)+
              FZ_AKT+FZ_FEB+
              chisq(AKT_c)+chisq(FEB_c)
            ,data =bdat_l[which(bdat_l$Sex=="Female"),],
            control = my_controls , simulate.p.value = TRUE , B=10000)

write2html(summary(t7, pfootnote=TRUE),"BehaviourPosR2_Full_Female.html")


t7<-tableby(outcome1~ 
              risk_R1+chisq(risk_R1_c)+
              Inf_risk_R1+ chisq(Inf_risk_R1_c)+
              Inf_grade_R1+chisq(Inf_grade_R1_c)+
              risk_R2+chisq(risk_R2_c)+
              Inf_risk_R2+chisq(Inf_risk_R2_c)+
              Inf_grade_R2+chisq(Inf_grade_R2_c)+
              sum_contact_R1+chisq(sum_contact_R2)+
              FZ_AKT+FZ_FEB+
              chisq(AKT_c)+chisq(FEB_c)
            ,data =bdat_l[which(bdat_l$Agegroup=="<35"),],
            control = my_controls , simulate.p.value = TRUE , B=10000)

write2html(summary(t7, pfootnote=TRUE),"BehaviourPosR2_Full_Agegroup_U35.html")

t7<-tableby(outcome1~ 
              risk_R1+chisq(risk_R1_c)+
              Inf_risk_R1+ chisq(Inf_risk_R1_c)+
              Inf_grade_R1+chisq(Inf_grade_R1_c)+
              risk_R2+chisq(risk_R2_c)+
              Inf_risk_R2+chisq(Inf_risk_R2_c)+
              Inf_grade_R2+chisq(Inf_grade_R2_c)+
              sum_contact_R1+chisq(sum_contact_R2)+
              FZ_AKT+FZ_FEB+
              chisq(AKT_c)+chisq(FEB_c)
            ,data =bdat_l[which(bdat_l$Agegroup=="35-65"),],
            control = my_controls , simulate.p.value = TRUE , B=10000)

write2html(summary(t7, pfootnote=TRUE),"BehaviourPosR2_Full_Agegroup_35_65.html")

t7<-tableby(outcome1~ 
              risk_R1+chisq(risk_R1_c)+
              Inf_risk_R1+ chisq(Inf_risk_R1_c)+
              Inf_grade_R1+chisq(Inf_grade_R1_c)+
              risk_R2+chisq(risk_R2_c)+
              Inf_risk_R2+chisq(Inf_risk_R2_c)+
              Inf_grade_R2+chisq(Inf_grade_R2_c)+
              sum_contact_R1+chisq(sum_contact_R2)+
              FZ_AKT+FZ_FEB+
              chisq(AKT_c)+chisq(FEB_c)
            ,data =bdat_l[which(bdat_l$Agegroup==">=65"),],
            control = my_controls , simulate.p.value = TRUE , B=10000)

write2html(summary(t7, pfootnote=TRUE),"BehaviourPosR2_Full_Agegroup_Over65.html")

bdat<- subset(koco_r2_dat, 
                select = c(ind_id,Sex,Age,status,
                           risk_R1,risk_R1_c,                              
                           Inf_risk_R1,Inf_risk_R1_c,                          
                           Inf_grade_R1,Inf_grade_R1_c,                         
                           risk_R2,risk_R2_c,                              
                           Inf_risk_R2,Inf_risk_R2_c,                          
                           Inf_grade_R2,Inf_grade_R2_c,                         
                           sum_contact_R1,sum_contact_R2,                         
                           FZ_AKT,FZ_FEB,                                 
                           AKT_c,FEB_c))
bdat$outcome<- paste0(bdat$status,bdat$Sex)
bdat$Agegroup<-ifelse(bdat$Age<35,"<35",
                        ifelse(bdat$Age>=65,">=65","35-65"))

bdat$Agegroup<- factor(bdat$Agegroup,levels = c("<35","35-65",">=65"))

bdat$outcome1<- paste0(bdat$status,bdat$Agegroup)
bdat$outcome1<- factor(bdat$outcome1,
                         levels = c("Never Positive<35","Never Positive35-65","Never Positive>=65",
                                    "Atleast once Positive<35","Atleast once Positive35-65","Atleast once Positive>=65"))

write.csv(bdat,"bdat_total.csv")
write.csv(bdat_l,"bdat_R2.csv")
t8<-tableby(outcome~ 
              risk_R1+chisq(risk_R1_c)+
              Inf_risk_R1+ chisq(Inf_risk_R1_c)+
              Inf_grade_R1+chisq(Inf_grade_R1_c)+
              risk_R2+chisq(risk_R2_c)+
              Inf_risk_R2+chisq(Inf_risk_R2_c)+
              Inf_grade_R2+chisq(Inf_grade_R2_c)+
              sum_contact_R1+chisq(sum_contact_R2)+
              FZ_AKT+FZ_FEB+
              chisq(AKT_c)+chisq(FEB_c)
            ,data =bdat[which(bdat$Sex=="Male"),],
            control = my_controls , simulate.p.value = TRUE , B=10000)

write2html(summary(t8, pfootnote=TRUE),"BehaviourEverPos_Full_Male.html")

t8<-tableby(outcome~ 
              risk_R1+chisq(risk_R1_c)+
              Inf_risk_R1+ chisq(Inf_risk_R1_c)+
              Inf_grade_R1+chisq(Inf_grade_R1_c)+
              risk_R2+chisq(risk_R2_c)+
              Inf_risk_R2+chisq(Inf_risk_R2_c)+
              Inf_grade_R2+chisq(Inf_grade_R2_c)+
              sum_contact_R1+chisq(sum_contact_R2)+
              FZ_AKT+FZ_FEB+
              chisq(AKT_c)+chisq(FEB_c)
            ,data =bdat[which(bdat$Sex=="Female"),],
            control = my_controls , simulate.p.value = TRUE , B=10000)

write2html(summary(t8, pfootnote=TRUE),"BehaviourEverPos_Full_Female.html")

t8<-tableby(outcome1~ 
              risk_R1+chisq(risk_R1_c)+
              Inf_risk_R1+ chisq(Inf_risk_R1_c)+
              Inf_grade_R1+chisq(Inf_grade_R1_c)+
              risk_R2+chisq(risk_R2_c)+
              Inf_risk_R2+chisq(Inf_risk_R2_c)+
              Inf_grade_R2+chisq(Inf_grade_R2_c)+
              sum_contact_R1+chisq(sum_contact_R2)+
              FZ_AKT+FZ_FEB+
              chisq(AKT_c)+chisq(FEB_c)
            ,data =bdat[which(bdat$Agegroup=="<35"),],
            control = my_controls , simulate.p.value = TRUE , B=10000)

write2html(summary(t8, pfootnote=TRUE),"BehaviourEverPos_Full_U35.html")

t8<-tableby(outcome1~ 
              risk_R1+chisq(risk_R1_c)+
              Inf_risk_R1+ chisq(Inf_risk_R1_c)+
              Inf_grade_R1+chisq(Inf_grade_R1_c)+
              risk_R2+chisq(risk_R2_c)+
              Inf_risk_R2+chisq(Inf_risk_R2_c)+
              Inf_grade_R2+chisq(Inf_grade_R2_c)+
              sum_contact_R1+chisq(sum_contact_R2)+
              FZ_AKT+FZ_FEB+
              chisq(AKT_c)+chisq(FEB_c)
            ,data =bdat[which(bdat$Agegroup=="35-65"),],
            control = my_controls , simulate.p.value = TRUE , B=10000)

write2html(summary(t8, pfootnote=TRUE),"BehaviourEverPos_Full_35_65.html")

t8<-tableby(outcome1~ 
              risk_R1+chisq(risk_R1_c)+
              Inf_risk_R1+ chisq(Inf_risk_R1_c)+
              Inf_grade_R1+chisq(Inf_grade_R1_c)+
              risk_R2+chisq(risk_R2_c)+
              Inf_risk_R2+chisq(Inf_risk_R2_c)+
              Inf_grade_R2+chisq(Inf_grade_R2_c)+
              sum_contact_R1+chisq(sum_contact_R2)+
              FZ_AKT+FZ_FEB+
              chisq(AKT_c)+chisq(FEB_c)
            ,data =bdat[which(bdat$Agegroup==">=65"),],
            control = my_controls , simulate.p.value = TRUE , B=10000)

write2html(summary(t8, pfootnote=TRUE),"BehaviourEverPos_Full_Over65.html")


