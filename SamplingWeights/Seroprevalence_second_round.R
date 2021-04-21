#' This script calibrates the sampling weights (from the script Sampling_weights.R)
#' for the second round of visit and calculates the variance associated with the estimation
#' of the proportion of ppl. infected by the virus (positive tested).
#' 

rm(list = ls())

here_weights = function(...) here::here("SamplingWeights", ...)


source(here_weights("Sampling_weights.R"))


here_koco_data = function (...) here::here("Data", ...)
here_weights = function(...) here::here("SamplingWeights", ...)

# CAUTION: For Calibration we need KoCo_BLab, therefore you have to run the script Sampling_weights.R first

library(descr)
library(reshape2)
library(sampling)
# library(ggplot2)


#############################
# Load data sets
#############################


###
# Second round of KoCo19 
###

KoCo19_sr <- read.csv(here_koco_data("Second_visit_Koco19/koco_r2_matched.csv"), stringsAsFactors = F)


###
# Complete sample of KoCo19 (with missing values)
###

# We need this data set to identify households with or without children
KoCo19 <- read.csv(here_koco_data("KoCo19_Datasets/ind_lab_baseline_new.csv"), stringsAsFactors = F)


# We have 134 missing ages in this data set. We consider that:
# - in households with 2 members, the missing age is another adult
# - in households with 3 or more members, the missing age(s) is/are a child(ren)


# In households with missing ages and with more than 2 members, we impute the missing ages with
# the value 10. It does not matter here, we will not use directly this value. We just want to say
# that these missing values are children.
KoCo19$Age[is.na(KoCo19$Age) & KoCo19$obs_hh_members > 2] <- 10

###
# Necessary datasets for the treatment of non response
###

# Collect information on respondents and non-respondents at the individual and at the household levels
info_ind <- read.csv(here_koco_data("KoCo19_Datasets/Koco_baseline.csv"), stringsAsFactors = TRUE)
info_hh <- read.csv(here_koco_data("KoCo19_Datasets/hh_characteristics_new.csv"), stringsAsFactors = TRUE)



###
# Data set with the constituencies
###

# reading the data set with with the consituencies info for each household
Const <- read.csv(here_koco_data("KoCo19_Datasets/KoCo19_Haushalte4Modeler_wRRstartConstituency_20200910.csv"), stringsAsFactors = TRUE)
# 3007 hh

# Remove some hh that are not in the final study population
Const <- Const[Const$hht_ID %in% KoCo_BLab$hh_id, ]
# 2994 hh

# Recoding of the IDs of the starting constituencies
Const$const_start <- paste0(substr(Const$constituency_RRstart, 1, 2), "0",
                            substr(Const$constituency_RRstart, 3, 4))

Const$const_start[Const$constituency_RRstart < 1000] <- 
  paste0("0", substr(Const$constituency_RRstart[Const$constituency_RRstart < 1000], 1, 1),
         "0", substr(Const$constituency_RRstart[Const$constituency_RRstart < 1000], 2, 3))

# Keep only the variables of interest
Const <- Const[, c("hht_ID", "const_start")]


###
# Data sets from the Statistisches Amt
###


### Sex/Age distribution

margins_age_sex <- read.csv(here_koco_data("KoCo19_Datasets/muc_age_structure_KoCo_bis_study_pop.csv"), stringsAsFactors = F)

# Rename variables to match the included information 
names(margins_age_sex)[names(margins_age_sex)=="Group.1"] <- "age_cat"
names(margins_age_sex)[names(margins_age_sex)=="Group.2"] <- "Sex"
names(margins_age_sex)[names(margins_age_sex)=="x"] <- "Freq"

# Changing the first letters of sex to match the information in the study data set
margins_age_sex$Sex[margins_age_sex$Sex=="male"] <- "Male"
margins_age_sex$Sex[margins_age_sex$Sex=="female"] <- "Female"

# Creating a variable including the information about the sex/age categories
margins_age_sex$sex_age <- paste(margins_age_sex$Sex, margins_age_sex$age_cat,
                                 sep = "_")

### Number of households per constituency

# Information about number of hh in Munich depending on constituency based on stat.Amt
Munich_hh <- read.csv(here_koco_data("KoCo19_Datasets/muc_hh_structure_KoCo_bis.csv"), stringsAsFactors = TRUE)

# Recoding of the IDs of the constituencies
Munich_hh$const <- as.character(Munich_hh$Const_ID)

Munich_hh$const[Munich_hh$Const_ID < 10000] <- paste0("0", Munich_hh$const[Munich_hh$Const_ID < 10000])

# Keep only the variables of interest
Munich_hh <- Munich_hh[, c("const", "Nb_hh")]



#############################
# Define sample for the second round
#############################

###
# Seroprevalence based on Roche and DBS values
# We shrink the sample to ppl for whom we have test results. 
# Recalculate the sampling weights within the households
##




###
# Shrink the sample to ppl for whom we have test results. Recalculate the sampling weights within the households
###

KoCo_BLab_sr <- merge(KoCo19_sr, KoCo_BLab, by = c("ind_id", "hh_id"), all.x =  TRUE)

### Recalculate the weights within the households (non response treatment, since some participants are now missing)

# Creating a variable with the number of people (ppl) in the households that are taking part in the study
KoCo_BLab_sr <- merge(KoCo_BLab_sr, as.data.frame(table(KoCo_BLab_sr$hh_id)), by.x = "hh_id", by.y = "Var1", all.x=T)
names(KoCo_BLab_sr)[names(KoCo_BLab_sr)=="Freq"] <- "n_ppl_part_hh"


# Calculating the weight for ppl in the household depending on other mm of the hh (not taking part in the study but in the study pop)
KoCo_BLab_sr$w_ind <- KoCo_BLab_sr$obs_hh_members_study / KoCo_BLab_sr$n_ppl_part_hh 

# Remove unnecessary variable
KoCo_BLab_sr$n_ppl_part_hh <- NULL


# Change name for first round
KoCo_BLab_fr <- KoCo_BLab

# Change name for second round
KoCo_BLab <- KoCo_BLab_sr


#############################
# Collecting auxiliary data for calibration (definition of the margins)
#############################

###
# Age/Sex structure and country of birth in Munich (based on KoCo19, stat. Amt)
###

# Selecting the margins for sex/age categories 
totals <- margins_age_sex$Freq
names(totals) <- margins_age_sex$sex_age

# Selecting the margins for country of origin 
country_germany_unknow <- 1085145
country_other <- sum(totals) - country_germany_unknow

###
# Nb of single/multi-ppl hh and hh with/without children in Munich (based on KoCo19, stat. Amt)
###

# Selecting the margins for household size (single vs. multi) 
# only main residences
n_single_house <- 449561
n_multi_house <- 828672 - n_single_house

# Selecting the margins for household with or without children in Munich
n_house_0_child <- 682708
n_house_1plus_child <- 828672 - n_house_0_child


###
# Summarizing the information about Munich
###
totals <- c(totals, country_germany_unknow = country_germany_unknow, country_other = country_other,
            n_single_house = n_single_house, n_multi_house = n_multi_house,
            n_house_0_child = n_house_0_child, n_house_1plus_child = n_house_1plus_child)

# Removing temporary data
rm(KoCo_BLab_sr, margins_age_sex, country_germany_unknow, country_other, n_single_house, n_multi_house, n_house_0_child, n_house_1plus_child)


#############################
# Create in the sample the variables needed for the calibration
#############################

###
# Age/Sex structure, country of origin 
###


# Creating age categories matching the information about Munich 
KoCo_BLab$age_cat <- cut(KoCo_BLab$Age, c(14, 19, 34, 49, 64, 79, 150),
                            include.lowest = TRUE,
                            right = TRUE,
                            labels = c("14-19", "20-34", "35-49", "50-64",
                                       "65-79", ">=80"))

# Creating a variable including the information about the sex/age categories
KoCo_BLab$sex_age <- paste(KoCo_BLab$Sex, KoCo_BLab$age_cat, sep = "_")

# Checking the nb of ppl in the sex/age categories
table(KoCo_BLab$sex_age, useNA = "always")

# Creating a variable including the information about hh_id and sex/age categories
KoCo_BLab$hh_sex_age <- paste(KoCo_BLab$hh_id, KoCo_BLab$sex_age, sep = "_")

table(KoCo_BLab$Birth_Country, useNA = "always")

# Creating a variable including the information about the country of origin
KoCo_BLab$Birth_Country[is.na(KoCo_BLab$Birth_Country)] <- "Germany"

# Checking the nb of ppl by country of origin
table(KoCo_BLab$Birth_Country, useNA = "always")

# Creating a variable including the information about hh_id and country categories
KoCo_BLab$hh_country <- paste(KoCo_BLab$hh_id, KoCo_BLab$Birth_Country, sep = "_")


### Sex/Age



# Calculating the number of males, females in each age group in each household
n_sex_age <- freq(KoCo_BLab$hh_sex_age, w = KoCo_BLab$w_ind, plot=F)

# Removing the line for totals
n_sex_age <- n_sex_age[-nrow(n_sex_age), ]

# Splitting the information saved in the variable about hh_id and sex/age categories
n_sex_age <- data.frame(hh_id = substr(rownames(n_sex_age), 1, 7),
                              sex_age = substr(rownames(n_sex_age), 9, 30),
                              freq = n_sex_age[, "Frequency"])

# Reshaping the dataframe to have for each hh the number of ppl based on sex and age
data_house <- reshape(n_sex_age, v.names = "freq",
                      timevar = "sex_age", idvar = "hh_id", direction = "wide")


### Country of origin

# Calculating the number Germans/other in each household
n_country <- freq(KoCo_BLab$hh_country, w = KoCo_BLab$w_ind, plot=F)

# Removing the line for totals
n_country <- n_country[-nrow(n_country), ]

# Splitting the information saved in the variable about hh_id and country of origin
n_country <- data.frame(hh_id = substr(rownames(n_country), 1, 7),
                              country = substr(rownames(n_country), 9, 15),
                              freq = n_country[, "Frequency"])


# Reshaping the dataframe to have for each hh the number of ppl based on contry of origin
data_house2 <- reshape(n_country, v.names = "freq",
                       timevar = "country", idvar = "hh_id", direction = "wide")

# Merging of the two information
data_house <- merge(data_house, data_house2, by = "hh_id")

# Setting empty cells to 0
data_house[is.na(data_house)] <- 0

# Renaming columns and rows 
colnames(data_house) <- gsub(x = colnames(data_house), pattern = "freq.", replacement = "")
rownames(data_house) <- NULL


# Reordering the columns
data_house <- data_house[, c("hh_id", "Male_14-19", "Male_20-34", "Male_35-49", "Male_50-64", "Male_65-79", "Male_>=80",
                             "Female_14-19", "Female_20-34", "Female_35-49", "Female_50-64", "Female_65-79", "Female_>=80",
                             "Germany", "Other")]



###
# Nb of members in the hh, children/no children in the hh 
###

# Adding information about number of observed individuals in hh
data_house <- merge(data_house, unique(KoCo_BLab[, c("hh_id", "obs_hh_members")]), by="hh_id", all.x = T)


# Adding information about type of housing (single vs. multi)
data_house$house <- cut(data_house$obs_hh_members, breaks = c(0, 1, 100),
                        labels = c("single", "multi"))

# Creating dummy variable(s) for type of housing
data_house <- cbind(data_house, model.matrix(~ house - 1, data = data_house))


# Adding information about children in the household (children vs. no children)
# Nb of children in each household
nb_child <- tapply(KoCo19$Age < 18, KoCo19$hh_id, function(x) sum(x, na.rm = T))
nb_child <- data.frame(hh_id = names(nb_child), nb_child = nb_child)

# Dummy variable: Child vs. no child
nb_child$child <- ifelse(nb_child$nb_child == 0, "no", "yes")

data_house <- merge(data_house, nb_child, by="hh_id", all.x = T)

# Creating dummy variable(s) for Children
data_house <- cbind(data_house, model.matrix(~ child - 1, data = data_house))

# Adding rownames
rownames(data_house) <- data_house$hh_id


# Removing unnecessary variables
data_house[, c("hh_id", "house", "obs_hh_members", "nb_child", "child")] <- NULL


# removing temporary data
rm(n_sex_age, n_country, KoCo19, data_house2, nb_child)

#############################
# Non response treatment
#############################

#############
# Uniform non response
#############

# p_rep <- length(unique((KoCo_BLab$hh_id))) / length(unique((KoCo_BLab_fr$hh_id)))

#############
# Estimated probability of NR
#############

# Link article Juillard, Chauvet:
# https://www150.statcan.gc.ca/n1/pub/12-001-x/2018002/article/54952-eng.htm

# Dummy variable of response at the household level
ind_nr <- as.numeric(unique(KoCo_BLab_fr$hh_id) %in% unique(KoCo_BLab$hh_id))

# Create a data frame data_nr to anakyse the NR behaviour
data_nr <- data.frame(hh_id = unique(KoCo_BLab_fr$hh_id), ind_nr = ind_nr)

# Add information at the household level
# Note: More information was included at the beginning, e.g. smoker in the household, living area, migration background, constituency, etc...
# However, these variables did not have an effect on the NR at the household level. Therefore, we do not include them here.
data_nr <- merge(data_nr, info_hh[, c("hh_id", "HousingType", "HouseholdSize",
                                      "NetIncome_month_hh", "HighestEducation_hh")], by = "hh_id", all.x = TRUE)

# Add information at the individual level (aggregated at the household level)
nb_pos_age <- by(info_ind, info_ind$hh_id, function(x){
  nb_pos <- sum(x[, "Roche_Result_new"] == "Positive", na.rm = TRUE)
  mean_age <- mean(x[, "Age"])
  return(c(nb_pos = nb_pos, mean_age = mean_age))
})

nb_pos_age <- as.data.frame(do.call(rbind, nb_pos_age))
nb_pos_age$hh_id <- rownames(nb_pos_age)

# Combine information at the individual level and at the household level
data_nr <- merge(data_nr, nb_pos_age, by = "hh_id")

###
# Process information at the individual level
###

# If at least 1 member was positive, household positive
data_nr$pos <- as.factor(ifelse(data_nr$nb_pos == 0, 0, 1))

# Categorize average age of the household
data_nr$mean_age <- cut(data_nr$mean_age, breaks = c(0, 34, 49, 64, 79, Inf),
                        labels = c("0-34", "35-49", "50-64", "65-79", "80+"))



###
# Process information at the household level
###

# Income
data_nr$income <- "No Answer"
data_nr$income[data_nr$NetIncome_month_hh %in% c("<500", "500-750", "750-1000", "1000-1500", "1500-2000", "2000-2500")] <- "<=2500"
data_nr$income[data_nr$NetIncome_month_hh %in% c("2500-3000", "3000-4000")] <- "2501-4000"
data_nr$income[data_nr$NetIncome_month_hh %in% c("4000-5000", "5000-6000")] <- "4001-6000"
data_nr$income[data_nr$NetIncome_month_hh == ">6000"] <- "6001+"
data_nr$income <- factor(data_nr$income, levels = c("No Answer", "<=2500", "2501-4000", "4001-6000", "6001+"))

# House type
data_nr$house_type <- "5+ apt"
data_nr$house_type[data_nr$HousingType %in% c("Type1", "Type2")] <- "1-2 apt"
data_nr$house_type[data_nr$HousingType == "Type3"] <- "3-4 apt"
data_nr$house_type[data_nr$HousingType == "Type7"] <- "Others"
data_nr$house_type <- as.factor(data_nr$house_type)

# Household size
data_nr$hh_members <- cut(data_nr$HouseholdSize, breaks = c(0, 1, 2, Inf),
                          labels = c("1", "2", "3+"))

# Size is missing for one household. We checked and household size = 1
data_nr[is.na(data_nr$hh_members), "hh_members"] <- 1


# Highest education
data_nr$education <- "Unknow"
data_nr$education[data_nr$HighestEducation_hh %in% c("Type1", "Type2", "Type3", "Type4", "Type7", "Type8", "Type9")] <- "<= 12 yrs"
data_nr$education[data_nr$HighestEducation_hh %in% c("Type5", "Type6")] <- "> 12 yrs"
data_nr$education <- factor(data_nr$education, levels = c("Unknow", "<= 12 yrs", "> 12 yrs"))

# No response
data_nr$ind_nr <- as.factor(data_nr$ind_nr)


###
# Model the non response to predict the probability of response
###

model <- glm(ind_nr ~ income + mean_age + house_type + hh_members + education + pos, data = data_nr, family = "binomial")
summary(model)

# Predicted probability of response
phat <- predict(model, type = "response")

###
# Create homogeneous group of response (HRG)
###

# We use here a very simple approach to create these groups by cutting the probabilities of response at each decile
# A more sophisticated approach was also applied (using k-means to cluster the probabilities), providing similar results.

hrg <- cut(phat, breaks = quantile(phat, probs = seq(0,1, 0.1)), include.lowest = TRUE)

data_nr$hrg <- hrg

# Estimate the probabilities in each HRG by number of respondents / number of participant
phat <- by(data_nr, data_nr$hrg, function(x){
  prop.table(table(x[, "ind_nr"]))[2]
})


# Add these probabilities to the KoCo_BLab table (participants to the follow-up)
phat <- do.call(rbind, as.list(phat))

phat <- data.frame(hrg = rownames(phat), phat = phat)

data_nr <- merge(data_nr, phat, by = "hrg")

KoCo_BLab <- merge(KoCo_BLab, data_nr[, c("hh_id", "phat")], by = "hh_id", all.x = TRUE)



#############################
# Calibration Roche DBS
#############################

# Checking if the order of age/sex categories are the same
colnames(data_house)
names(totals)

# Weights at the household level (corrected for the NR)
d_house <- KoCo_BLab$w_constituency * KoCo_BLab$w_household / KoCo_BLab$phat
names(d_house) <- KoCo_BLab$hh_id
d_house <- d_house[!duplicated(names(d_house))]

# Ordering weights in the same order than the household data
d_house <- d_house[rownames(data_house)]

###
# Calculating the calibrated weights by means of different calibration methods.
###


### method "logit"
g_logit <- calib(data_house, d = d_house, totals, method="logit",
           bounds=c(0.3,3), max_iter = 2000)  
w_hh_cal_logit <- g_logit * d_house

# Checking if calibrated weights sum up to auxiliary data
sum(t(w_hh_cal_logit) %*% as.matrix(data_house) - totals[])
summary(w_hh_cal_logit)
plot(density(g_logit))
# Ok


calib_weights <- data.frame(hh_id = names(d_house), w_hh_cal = w_hh_cal_logit,
                            g = g_logit)

KoCo_BLab <- merge(KoCo_BLab, calib_weights, by = "hh_id")

KoCo_BLab$w_ind_cal <- KoCo_BLab$w_hh_cal * KoCo_BLab$w_ind




#############################
# Calibrated estimators
#############################

# Estimation of the proportion with/without the calibrated weights

freq(x = KoCo_BLab$R2_Result, plot=F)
freq(x = KoCo_BLab$R2_Result, w = KoCo_BLab$w_ind_cal, plot=F)

sum(KoCo_BLab$w_ind_cal)

# Estimated number of ppl tested positive during the first round
freq(x = KoCo_BLab$Roche_Result_new, w = KoCo_BLab$w_ind_cal, plot=F)
# 22064

