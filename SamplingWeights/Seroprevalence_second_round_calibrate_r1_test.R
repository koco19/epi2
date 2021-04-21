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
# Nb ppl. tested positive during the first round
###

nb_pos_r1 <- 25900
nb_neg_r1 <- sum(totals) - nb_pos_r1


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
            Negativ = nb_neg_r1, Positiv = nb_pos_r1,
            n_single_house = n_single_house, n_multi_house = n_multi_house,
            n_house_0_child = n_house_0_child, n_house_1plus_child = n_house_1plus_child)

# Removing temporary data
rm(KoCo_BLab_sr, margins_age_sex, country_germany_unknow, country_other, n_single_house, n_multi_house, n_house_0_child, n_house_1plus_child, nb_pos_r1)


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

# Checking the nb of ppl who were tested positive in R1
table(KoCo_BLab$Roche_Result_new, useNA = "always")

# Creating a variable including the information about hh_id and positive R1
KoCo_BLab$hh_pos_r1 <- paste(KoCo_BLab$hh_id, KoCo_BLab$Roche_Result_new, sep = "_")



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

### Tested positive Round 1

# Calculating the number pos R1 / neg R1 in each household
n_pos <- freq(KoCo_BLab$hh_pos_r1, w = KoCo_BLab$w_ind, plot=F)

# Removing the line for totals
n_pos <- n_pos[-nrow(n_pos), ]

# Splitting the information saved in the variable about hh_id and country of origin
n_pos <- data.frame(hh_id = substr(rownames(n_pos), 1, 7),
                        pos_r1 = substr(rownames(n_pos), 9, 15),
                        freq = n_pos[, "Frequency"])


# Reshaping the dataframe to have for each hh the number of ppl based on contry of origin
data_house2 <- reshape(n_country, v.names = "freq",
                       timevar = "country", idvar = "hh_id", direction = "wide")

# Merging of the two information
data_house <- merge(data_house, data_house2, by = "hh_id")

# Reshaping the dataframe to have for each hh the number of ppl tested positive
data_house2 <- reshape(n_pos, v.names = "freq",
                       timevar = "pos_r1", idvar = "hh_id", direction = "wide")

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
                             "Germany", "Other", "Negativ", "Positiv")]



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
rm(n_sex_age, n_country, n_pos, KoCo19, data_house2, nb_child)



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

freq(x = KoCo_BLab$R2_Result[KoCo_BLab$Roche_Result_new == "Negative"], w = KoCo_BLab$w_ind_cal[KoCo_BLab$Roche_Result_new == "Negative"], plot=F)


#############################
# Exporting file with the weights
#############################

### Adding the weights that sum up to the size of the sample

KoCo_Weights <- KoCo_BLab

# At the household level

# Roche
tab_hh <- KoCo_Weights[!duplicated(KoCo_Weights$hh_id), "w_hh_cal" ]

KoCo_Weights$w_hh_samp <- KoCo_Weights$w_hh_cal * length(tab_hh) / sum(tab_hh)


# Check if the sum is correct
sum(KoCo_Weights[!duplicated(KoCo_Weights$hh_id), "w_hh_samp" ])

# At the individual level
KoCo_Weights$w_ind_samp <- KoCo_Weights$w_ind_cal * nrow(KoCo_Weights) / sum(KoCo_Weights$w_ind_cal)

# Check if the sum is correct
sum(KoCo_Weights$w_ind_samp)


### Keep the variables of interest
KoCo_Weights <- KoCo_Weights[, c("hh_id", "ind_id", "Age", "Sex", "Birth_Country",
                                 "obs_hh_members", "Roche_Result_new", "R2_Result",
                                 "w_hh_cal", "w_ind_cal", "w_hh_samp", "w_ind_samp")]

### Export table
write.csv(KoCo_Weights, here_koco_data("Second_visit_Koco19/KoCo_weights.csv"), row.names = FALSE)

# Removing temporary data
rm(list = setdiff(ls(), c("here_weights", "here_koco_data", "KoCo_weights", "KoCo_BLab", "data_house", "Const", "Munich_hh",
                          "d_house")))


#############################
# Variance of the calibrated estimator
#############################


###
# calculate the linearized variable for a proportion
###

### Linearization at the individual level

# Dummy variables for DBS
table(KoCo_BLab$R2_Result)

KoCo_BLab$DBS_results <- ifelse(KoCo_BLab$R2_Result == "Positive", 1, 0)


# Weighted proportion for DBS
p_dbs <- sum(KoCo_BLab$DBS_results * KoCo_BLab$w_ind_cal)/sum(KoCo_BLab$w_ind_cal)


# Linearized variables
KoCo_BLab$dbs_u_k_0 <- (KoCo_BLab$DBS_results - p_dbs)/sum(KoCo_BLab$w_ind_cal)


### Linearization for the two groups: positive at the first round and the negative at the first round

group <- split(KoCo_BLab, KoCo_BLab$Roche_Result_new)

# Weighted proportions for DBS in each group:
# - Negative at baseline: proportion of positive in the second round
# - Positive at baseline: Proportions of positive AND negative in the second round

p_dbs <- sapply(group, function(z){
  pos_r2 <- sum(z[, "DBS_results"] * z[, "w_ind_cal"])/sum(z[, "w_ind_cal"])
  neg_r2 <- sum((1-z[, "DBS_results"]) * z[, "w_ind_cal"])/sum(z[, "w_ind_cal"])
  return(c(pos_r2 = pos_r2, neg_r2 = neg_r2))
})

p_dbs <- c(Negative_Positive = p_dbs[1], Positive_Positive = p_dbs[3], Positive_Negative = p_dbs[4])

group[[3]] <- group[[2]]
group[[3]][, "DBS_results"] <- 1 - group[[3]][, "DBS_results"]

# Linearized variables

KoCo_BLab_group <- lapply(1:length(group), function(x){
  group[[x]][, paste0("dbs_u_k_", x)] <- (group[[x]][, "DBS_results"] - p_dbs[x])/sum(group[[x]][, "w_ind_cal"])
  return(group[[x]])
})

# Add linearized variables Positive_Negative to Positive_Positive
KoCo_BLab_group[[3]][, "DBS_results"] <- 1 - KoCo_BLab_group[[3]][, "DBS_results"]
KoCo_BLab_group[[2]] <- merge(KoCo_BLab_group[[2]], KoCo_BLab_group[[3]])

# Combine datasets
KoCo_BLab_group <- do.call(dplyr::bind_rows, KoCo_BLab_group[1:2])


# Change NAs to zero
KoCo_BLab_group[, paste0("dbs_u_k_", 1:length(group))][is.na(KoCo_BLab_group[, paste0("dbs_u_k_", 1:length(group))])] <- 0


### Linearization at the household level

split_hh <- split(KoCo_BLab_group[, c("hh_id", "w_ind", paste0("dbs_u_k_", 0:length(group)))], KoCo_BLab_group$hh_id)

hh_u_k <- sapply(split_hh, function(z){
  res <- sapply(0:length(group), function(x){
    sum(z[, paste0("dbs_u_k_", x)] * z[, "w_ind"])
  })
  return(res)
})


hh_u_k <- as.data.frame(t(hh_u_k))

hh_u_k$hh_id <- rownames(hh_u_k)



###
# Create a data frame at the hh level with the linearized variable,
# the auxiliary information, the calibrated weights and the ratio
# of the calibrated weights and the sampling weights g
###

data_house$hh_id <- rownames(data_house)

data_house <- merge(data_house, KoCo_BLab[!duplicated(KoCo_BLab$hh_id),
                                          c("hh_id", "w_hh_noshare", "w_household", "w_hh_cal", "g", "phat")],
                    by = "hh_id")


data_house <- merge(data_house, hh_u_k, by = "hh_id")

identical(data_house$hh_id, names(d_house))

###
# Run a linear model with the linearized variable as response variable and the
# auxiliary information as covariates including the sampling weights
###

res.lm <- lapply(0:length(group), function(x){
  # x <- 1
  res.lm <- lm(as.formula(paste(paste0("V",x+1), paste(paste0("`", names(data_house[, 2:19]), "`"), collapse=" + "), sep=" ~ ")), 
               weights = d_house,
               data = data_house[, c(2:21, ncol(data_house) - length(group) + x)])
})

###############
# Variance estimation
###############

# For the variance estimation, we will use the residuals e_k -> ge


data_res <- cbind(data_house, sapply(res.lm, function(x){
  x$residuals
}))

ge <- sapply(0:length(group), function(x){
  data_res[, paste0("ge_dbs_", x)] <- data_res[, ncol(data_res) - length(group) + x] 
})


colnames(ge) <- paste0("ge_dbs_", 0:length(group))

# Keep useful variables
data_res <- cbind(data_res[, c("hh_id", "w_hh_noshare", "w_household" , "phat")], ge)


# Add the constituency ID
data_res <- merge(data_res, Const, by.x = "hh_id", by.y = "hht_ID")

########
### First term of the variance estimation
########

# Total number of constituencies
N_const <- length(unique(Munich_hh$const))

# number of constituencies surveyed
n_const <- length(unique(Const$const_start))


# Calculate the estimated totals of ge at the constituency level
split_const <- split(data_res, data_res$const_start)

tot_e_const <- sapply(split_const, function(z){
  res <- sapply(0:length(group), function(x){
    sum(z[, paste0("ge_dbs_", x)] * z[, "w_household"] / z[, "phat"] )
  })
  return(res)
})


# Calculate the dispersion of ge
s_e_dbs <- apply(tot_e_const, 1, var)


# First term of the variance estimation
V1_dbs <- N_const^2 * (1/n_const - 1/N_const) * s_e_dbs


########
### Second term of the variance estimation including the non response term
########


# Own function for the variance including the probability of response
# (inspired by varHT from the library sampling)
varHT_home <- function (y, pikl, p_rep, method = 1) 
{
  if (any(is.na(pikl))) 
    stop("there are missing values in pikl")
  if (any(is.na(p_rep))) 
    stop("there are missing values in p_rep")
  if (any(is.na(y))) 
    stop("there are missing values in y")
  if (!(is.data.frame(pikl) | is.matrix(pikl))) 
    stop("pikl should be a matrix or a data frame")
  if (!(is.data.frame(p_rep) | is.matrix(p_rep))) 
    stop("p_rep should be a matrix or a data frame")
  if (is.data.frame(pikl) | is.matrix(pikl)) 
    if (nrow(pikl) != ncol(pikl)) 
      stop("pikl is not a square matrix")
  if (is.data.frame(p_rep) | is.matrix(p_rep)) 
    if (nrow(p_rep) != ncol(p_rep)) 
      stop("p_rep is not a square matrix")
  if (length(y) != nrow(pikl)) 
    stop("y and pik have different sizes")
  if (length(y) != nrow(p_rep)) 
    stop("y and p_rep have different sizes")
  if (!missing(method) & !(method %in% c(1, 2))) 
    stop("the method should be 1 or 2")
  if (is.data.frame(pikl)) 
    pikl = as.matrix(pikl)
  if (is.data.frame(p_rep)) 
    p_rep = as.matrix(p_rep)
  pik = diag(pikl)
  pik1 = outer(pik, pik, "*")
  delta = pikl - pik1
  diag(delta) = pik * (1 - pik)
  y1 = outer(y, y, "*")
  if (method == 1) 
    return(sum(y1 * delta/(pik1 * pikl * p_rep)))
  if (method == 2) {
    y2 = outer(y/pik, y/pik, "-")^2
    return(0.5 * sum(y2 * (pik1 - pikl)/(pikl * p_rep)))
  }
}


# Nb hh surveyed
nb_hh_s <- as.data.frame(table(Const$const_start))

# Add to the data the nb of hh in each constituency and the number of hh surveyed
data_res <- merge(data_res, nb_hh_s, by.x = "const_start", by.y = "Var1")
data_res <- merge(data_res, Munich_hh, by.x = "const_start", by.y = "const", all.x = TRUE)

# Split data based on the constituencies (in order to calculate the second order inclusion probabilities)
data_res_split <- split(data_res, data_res$const_start)

### Calculate the second term of the variance associated to the sampling design and to the non response
V2_dbs <- sapply(data_res_split, function(x){
  # Matrix containing the second order inclusion probabilities of each hh in consituency x
  pi2 <- matrix(unique(x[, "Freq"] * (x[, "Freq"] - 1)) / (x[, "Nb_hh"] * (x[, "Nb_hh"] - 1)), nrow = nrow(x), ncol = nrow(x))
  diag(pi2) <- unique(x[, "Freq"] / x[, "Nb_hh"])
  # Matrix containig the single and joint probabilities of response of each household in constituency x
  mat_p_rep <- outer(x[, "phat"], x[, "phat"], "*")
  diag(mat_p_rep) <- x[, "phat"]
  # For each variable of interest, we calculate the second term of the variance in constituency x
  V2_dbs <- sapply(0:length(group), function(y){
    # Variance associated to the sampling design for variable y
    var_ht <- varHT_home(x[, paste0("ge_dbs_", y)], pi2, mat_p_rep, method = 1)
    # Variance associated to the NR for variable y
    var_nr <- sum((x[, paste0("ge_dbs_", y)] / diag(pi2))^2 * (1 - diag(mat_p_rep))/diag(mat_p_rep)^2)
    return(sum(var_ht, var_nr))
  }
  )
  
  return(V2_dbs)
})

# Second term of the variance

V2_dbs <- N_const / n_const * rowSums(V2_dbs)


###
# Variance estimation and confidence intervals
###

### DBS


w_p_dbs <- sapply(group, function(z){
  freq(z[, "DBS_results"], w = z[, "w_ind_cal"], plot = F)[1:2, 2]
})

w_p_dbs <- c(All = freq(x = KoCo_BLab$DBS_results, w = KoCo_BLab$w_ind_cal, plot=F)[2, 2],
             Negative_Positive = w_p_dbs[2], Positive_Positive = w_p_dbs[4], Positive_Negative = w_p_dbs[3])

V_dbs <- V1_dbs + V2_dbs

# Upper bound CI
w_u_ci_dbs <- w_p_dbs + qnorm(0.975) * sqrt(V_dbs) * 100

# Lower bound CI
w_l_ci_dbs <- w_p_dbs - qnorm(0.975) * sqrt(V_dbs) * 100


res_w <- data.frame(Adjustment = "Unadjusted", Calculation = "Weighted", 
                    Categories = c("Prevalence R2", "Negative R1 Positive R2", "Positive R1 Positive R2", "Positive R1 Negative R2"), 
                  Estimates = w_p_dbs, Lower_95_CB = w_l_ci_dbs, Upper_95_CB = w_u_ci_dbs,
                  row.names = NULL)


###
# Adjust for specificity and sensitivity
###

spec_sens_classifier <- readRDS(file = here_koco_data("KoCo19_Datasets/specificity_sensitivity_classifier.RData"))

res_adj <- res_w

res_adj$Adjustment <- "Adjusted"

res_adj[res_adj$Categories != "Positive R1 Negative R2", c("Estimates", "Lower_95_CB", "Upper_95_CB")] <- 
  (res_adj[res_adj$Categories != "Positive R1 Negative R2", c("Estimates", "Lower_95_CB", "Upper_95_CB")]/100 + spec_sens_classifier["Roche N pan-Ig optimized cut-off", "Specificity"] - 1)/(spec_sens_classifier["Roche N pan-Ig optimized cut-off", "Sensitivity"] + spec_sens_classifier["Roche N pan-Ig optimized cut-off", "Specificity"] - 1)*100

res_adj[res_adj$Categories == "Positive R1 Negative R2", c("Estimates", "Lower_95_CB", "Upper_95_CB")] <- 
  (res_adj[res_adj$Categories == "Positive R1 Negative R2", c("Estimates", "Lower_95_CB", "Upper_95_CB")]/100 + spec_sens_classifier["Roche N pan-Ig optimized cut-off", "Sensitivity"] - 1)/(spec_sens_classifier["Roche N pan-Ig optimized cut-off", "Sensitivity"] + spec_sens_classifier["Roche N pan-Ig optimized cut-off", "Specificity"] - 1)*100



###
# All results
###

res_all <- rbind(res_w, res_adj)

write.csv(res_all, here_weights("Estimates_Cal_R1_Test.csv"), row.names = FALSE)

