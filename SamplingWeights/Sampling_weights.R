#' This script calculates the sampling weights for confounding effects based 
#' on the sampling design (selection of the constituencies, the households,
#' and the individuals) during the first round of interviews.
#' The weights of the households are shared to account
#' for the random routes crossing the boundaries of the constituencies.



rm(list = ls())

here_koco_data = function (...) here::here("Data", ...)
here_weights = function(...) here::here("SamplingWeights", ...)


#############################
# Load data sets
#############################

###
# Individual data sets
###

# Reading the lab data set with all observations from the first round
KoCo_BLab <- read.csv(here_koco_data("KoCo19_Datasets/Koco_baseline.csv"), stringsAsFactors = TRUE)

table(KoCo_BLab$Roche_Result_old)

# Keep only the variables of interest
KoCo_BLab <- KoCo_BLab[, c("ind_id", "hh_id", "Age", "obs_hh_members")]

# We need this other data set to calculate the number of ppl in the study population (age >= 14) in each hh
KoCo19 <- read.csv(here_koco_data("KoCo19_Datasets/ind_lab_baseline_new.csv"), stringsAsFactors = F)

# We have 134 missing ages in this data set. We consider that:
# - in households with 2 members, the missing age is another adult (age >= 14)
# - in households with 3 or more members, the missing age(s) is/are a child(ren) under 14 years old

 
# In each household, calculate the number of members aged 14 or more
member_study_hh <- as.data.frame(table(KoCo19$hh_id[KoCo19$Age >= 14]))

colnames(member_study_hh) <- c("hh_id", "obs_hh_members_study")

KoCo19 <- merge(KoCo19, member_study_hh)

# Identify households with missing ages
missing_age <- KoCo19[is.na(KoCo19$Age), ]

# In the households of size 2 with missing ages, change the number of members in 
# the study population to 2
KoCo19$obs_hh_members_study[KoCo19$obs_hh_members == 2 & KoCo19$hh_id %in% missing_age$hh_id] <- 2

# Add the variable obs_hh_members_study to the data set KoCo_BLab
obs_hh_members_study <- KoCo19[!duplicated(KoCo19$hh_id), c("hh_id", "obs_hh_members_study")]

KoCo_BLab <- merge(KoCo_BLab, obs_hh_members_study, all.x = TRUE)

###
# Data set with the constituencies
###

# reading the data set with with the consituencies info for each household
Const <- read.csv(here_koco_data("KoCo19_Datasets/KoCo19_Haushalte4Modeler_wRRstartConstituency_20200910.csv"))
# 3007 hh

# Remove some hh that are not in the final study population
Const <- Const[Const$hht_ID %in% KoCo_BLab$hh_id, ]
# 2994 hh

# Recoding of the IDs of the ending constituencies
Const$const_end <- paste0(substr(Const$constituency, 1, 2), "0",
                          substr(Const$constituency, 3, 4))

Const$const_end[Const$constituency < 1000] <- 
  paste0("0", substr(Const$constituency[Const$constituency < 1000], 1, 1),
         "0", substr(Const$constituency[Const$constituency < 1000], 2, 3))

# Recoding of the IDs of the starting constituencies
Const$const_start <- paste0(substr(Const$constituency_RRstart, 1, 2), "0",
                            substr(Const$constituency_RRstart, 3, 4))

Const$const_start[Const$constituency_RRstart < 1000] <- 
  paste0("0", substr(Const$constituency_RRstart[Const$constituency_RRstart < 1000], 1, 1),
         "0", substr(Const$constituency_RRstart[Const$constituency_RRstart < 1000], 2, 3))

# Keep only the variables of interest
Const <- Const[, c("hht_ID", "const_end", "const_start")]

###
# Data set from the Statistisches Amt
###

### Number of households per constituency

# Information about number of hh in Munich depending on constituency based on stat.Amt
Munich_hh <- read.csv(here_koco_data("KoCo19_Datasets/muc_hh_structure_KoCo_bis.csv"))

# Recoding of the IDs of the constituencies
Munich_hh$const <- as.character(Munich_hh$Const_ID)

Munich_hh$const[Munich_hh$Const_ID < 10000] <- paste0("0", Munich_hh$const[Munich_hh$Const_ID < 10000])

# Keep only the variables of interest
Munich_hh <- Munich_hh[, c("const", "Nb_hh")]


###
# Contiguity matrices
###

# First order neighborhood of the constituencies
load(here_weights("contiguity_matrices-munich-constituencies_200914.RData"))

# First and second order neighborhood of the constituencies (neighbors of the neighbors)
contiguity_Mat <- matrix(0, nrow = 755, ncol = 755)

for(i in 1:755){
  for(j in 1:755){
    if(M.binary[i, j] == 1){
      for(k in 1:755){
        if(M.binary[k, j] == 1){
          contiguity_Mat[i, k] <- 1
        }
      }
    }
  }
}

# Changing the rownames
rownames(contiguity_Mat) <- paste0(substr(rownames(M.binary), 1, 2), "0",
                                   substr(rownames(M.binary), 3, 4)) 

rownames(contiguity_Mat)[as.numeric(rownames(M.binary)) < 1000] <- 
  paste0("0", substr(rownames(M.binary)[as.numeric(rownames(M.binary)) < 1000], 1, 1),
         "0", substr(rownames(M.binary)[as.numeric(rownames(M.binary)) < 1000], 2, 3))

# Changing the colnames
colnames(contiguity_Mat) <- paste0(substr(colnames(M.binary), 1, 2), "0",
                                   substr(colnames(M.binary), 3, 4)) 

colnames(contiguity_Mat)[as.numeric(colnames(M.binary)) < 1000] <- 
  paste0("0", substr(colnames(M.binary)[as.numeric(colnames(M.binary)) < 1000], 1, 1),
         "0", substr(colnames(M.binary)[as.numeric(colnames(M.binary)) < 1000], 2, 3))

# Reorder the rows and the columns of the matrix
contiguity_Mat <- contiguity_Mat[order(rownames(contiguity_Mat)), order(colnames(contiguity_Mat))]


#############################
# Sampling weights of the constituencies
#############################

# Weight of constituency based an random sampling (nb of all constituencies / nb of sampled constituencies)
w_constituency <- length(unique(Munich_hh$const)) / length(unique(Const$const_start))

# Add this weight to the individual data
KoCo_BLab$w_constituency <- w_constituency

#############################
# Conditional sampling weights of the households
#############################

# Creating a table with, for each constituency, the number of households (hh) that have been surveyed
# and the number of hh that would have been surveyed is the constituency was selected (30)
w_hh <- merge(Munich_hh, as.data.frame(table(Const$const_start)),
              by.x = "const", by.y = "Var1", all.x=T)
names(w_hh)[names(w_hh)=="Freq"] <- "n_samp_hh"

# Indicator that the constituency was selected
w_hh$samp <- 0
w_hh$samp[!is.na(w_hh$n_samp_hh)] <- 1

# For the constituencies that were not drawn, we would have selected 30 households
w_hh$n_samp_hh[is.na(w_hh$n_samp_hh)] <- 30

# Weights of the households
w_hh$w_hh <- w_hh$Nb_hh/w_hh$n_samp_hh


#############################
# Sharing the weights of the households
#############################


# Nb of neighbors
nb <- rowSums(contiguity_Mat)

# Normalized matrix
contiguity_Mat <- contiguity_Mat / nb

# Reorder the table with the weights of each hh in each constituency
w_hh <- w_hh[order(w_hh$const), ]

# For each selected consituency (samp = 1), we will share the weights with the 
# neighboring ones (first and second order neighbors)
shared_w_hh <- contiguity_Mat %*% w_hh$w_hh 


# Add the initial weights of the households and shared weights to the hh data
shared_w_hh <- data.frame(const = rownames(shared_w_hh), w_household = shared_w_hh)
Const <- merge(Const, w_hh[, c("const", "w_hh")], by.x = "const_start", by.y = "const", all.x = TRUE)
Const <- merge(Const, shared_w_hh, by.x = "const_start", by.y = "const", all.x = TRUE)

# Add the shared weights to the individual data
KoCo_BLab <- merge(KoCo_BLab, Const[, c("hht_ID", "w_hh", "w_household")], by.x = "hh_id", by.y = "hht_ID", all.x = TRUE)

# Rename variable
names(KoCo_BLab)[names(KoCo_BLab) == "w_hh"] <- "w_hh_noshare"

#############################
# Conditional sampling weights at the individual level (some observations were removed due to missingness in the data)
#############################

# Creating a variable with the number of people (ppl) in the households that are taking part in the study
KoCo_BLab <- merge(KoCo_BLab, as.data.frame(table(KoCo_BLab$hh_id)), by.x = "hh_id", by.y = "Var1", all.x=T)
names(KoCo_BLab)[names(KoCo_BLab)=="Freq"] <- "n_ppl_part_hh"


# Calculating the weight for ppl in the household depending on other mm of the hh (not taking part in the study but in the study pop)
KoCo_BLab$w_ind <- KoCo_BLab$obs_hh_members_study / KoCo_BLab$n_ppl_part_hh 
# w_ind_study is the weight ppl have for their household (because of other mm not taking part in the study) 
# e.g the ppl with weight 6 represent additionally 5 others of their hh.

# Remove unnecessary variable
KoCo_BLab$n_ppl_part_hh <- NULL

#############################
# Final(non conditional) sampling weights at the individual level
#############################

KoCo_BLab$w_ind_samp <- KoCo_BLab$w_constituency * KoCo_BLab$w_household * KoCo_BLab$w_ind




#############################
# Tidy
#############################

KoCo_BLab <- KoCo_BLab[setdiff(names(KoCo_BLab), "Age")]

rm(list = setdiff(ls(), "KoCo_BLab"))
