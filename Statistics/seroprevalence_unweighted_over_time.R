#' This script calculates the unweighted point estimates and associated confidence intervals
#' (unadjusted and adjusted for test performance) for overall and weekly prevalences 
#' for the second round of visit of the KoCo19 project. Numerical results are stored in the "Output" folder 
#' as .xlsx files and the pictures as .pdf/.png files.

rm(list = ls())

# location of data
here_koco_data = function (...)
  here::here("Data", ...)

# location of output files
here_output = function (...)
  here::here("Output", ...)

# function to install and load packages
Install_And_Load <- function(Required_Packages)
{
  Remaining_Packages <-
    Required_Packages[!(Required_Packages %in% installed.packages()[, "Package"])]
  
  if (length(Remaining_Packages)) {
    install.packages(Remaining_Packages)
  }
  for (package_name in Required_Packages) {
    library(package_name,
            character.only = TRUE,
            quietly = TRUE)
  }
}

# packages needed
Required_Packages <-
  c("tidyverse",
    "magrittr",
    "mgcv",
    "xlsx",
    "data.table",
    "writexl")

# install (if necessary) and load packages
Install_And_Load(Required_Packages)

# number of bootstrap samples for cluster bootstrap
boot <- 10000

# should intermediates be dropped? If FALSE, results of the main body of the manuscript are created,
# if TRUE, appendix figure is created
drop_intermediates <- FALSE

# test performance (from Olbricht et al, 2020)
specificity <- 0.9972041
sensitivity <- 0.8860104

# Load data
d_r2 <- read.csv(here_koco_data("Second_visit_Koco19/koco_r2_dat.csv"),
                stringsAsFactors = TRUE) %>%
  as_tibble()

# select relevant variables
d_r2 %<>%
  dplyr::select(
    ind_id,
    hh_id,
    contains("BloodID"),
    VisitDate,
    contains("esult"),
    Roche_R2_new,
    contains("quant"),
    date_R2,
    starts_with("Date"),
    -IgG_result,
    -R_Result
  )

# create data set without missing values for second round
d_r2_full <- d_r2 %>%
  filter(!is.na(R2_Result))

# unweighted seroprevalence -----------------------------------------------------------

# crude number of positives and intermediates at second round 
# according to DBS or venous blood sampling if DBS intermediate
d_r2_full %>% dplyr::select(R2_Result) %>% table()

# seroprevalence
sp_crude <- prop.table(table(d_r2$R2_Result))[2]
sp_crude

## cluster bootstrap for confidence intervals

seroprevalence_data <- d_r2_full %>%
  # select only household id and test result variable
  dplyr::select(hh_id, R2_Result) %>%
  # recode test result variable
  mutate(R2_Result = ifelse(R2_Result == "Positive", 1, 0))  %>%
  as.data.table()

# placeholder for bootstrap estimates
pos_rate_clu <- numeric(boot)
# cluster by household id
setkey(seroprevalence_data, "hh_id")
# store number of households
n_hh <- length(unique(seroprevalence_data$hh_id))

# set seed for replicability
set.seed(22)
# start loop
for (i in 1:boot)
{
  # create the bootstrap dataset by resampling clusters
  samp <- sample(unique(seroprevalence_data$hh_id), n_hh, replace = TRUE)
  bootdat <- seroprevalence_data[J(samp), allow.cartesian = TRUE]
  # calculate bootstrap estimate
  pos_rate_clu[i] <- mean(bootdat$R2_Result == 1)
}
# results: point estimate and 95% percentile intervals
c(mean(pos_rate_clu), quantile(pos_rate_clu, probs = c(0.025, 0.975)))

# save results: crude...
prev <- tibble(
  Prevalence = c("Mean", "Lower Bound", "Upper Bound"),
  crude = c(sp_crude, quantile(pos_rate_clu, probs = c(0.025, 0.975)))
) %>%
  #...and adjusted prevalence
  mutate(adjusted = (crude + specificity - 1) / (specificity + sensitivity - 1))
prev

write_xlsx(prev, here_output("prevalence_unweighted_bootstrap.xlsx"))


# Seroprevalence over time -------------------------------------------------
# look at DBS results over time
d_r2_full %>%
  # reformat date variable
  mutate(
    week = as.Date(Date_DBS_prick),
    # group by week
    week = cut.Date(week, breaks = "1 week", labels = FALSE)
  ) %>%
  group_by(week) %>%
  # count tests per week
  mutate(obs_per_week = n(),
         # count number of households tested per week
         hh_per_week = length(unique(hh_id))) %>%
  # show tests, number of households, number of positive test results and 
  # positive rate per week as well as
  # number and rate of intermediates
  summarise(
    obs = mean(obs_per_week),
    n_hh = mean(hh_per_week),
    pos_obs = sum(R2_Result == "Positive"),
    pos_rate = mean(R2_Result == "Positive"),
    int_obs = sum(DBS_Result == "Intermediate"),
    int_rate = mean(DBS_Result == "Intermediate")
  )

# use DBS date or venous blood sampling date for the intermediates
# recode two entries that reported finger prick before it was actually possible
d_r2_full %<>%
  mutate(
    week = as.Date(date_R2),
    # recode two entries that reported finger prick before it was 
    # actually possible by shifting weeks
    # and recoding the one that were to early ("in week 0")
    week = cut.Date(week, breaks = "1 week", labels = FALSE) - 1,
    week = ifelse(week == 0, 1, week)
  ) %>%
  group_by(week) %>%
  # count tests and households per week
  mutate(obs_per_week = n(),
         hh_per_week = length(unique(hh_id))) %>%
  ungroup()

# for appendix: drop intermediates as they would introduce a "bias" since these persons are 
# retested at the end of the sampling period with venous
# blood sampling and are more likely to be positive
if (drop_intermediates)
prev_time <- d_r2_full  %>%
  filter(DBS_Result != "Intermediate") else
    prev_time <- d_r2_full

# show tests and number of positives per week;
# see whether share of intermediates changes over time
prev_time %>%
  group_by(week) %>%
  summarise(
    obs = mean(obs_per_week),
    n_hh = mean(hh_per_week),
    pos_obs = sum(R2_Result == "Positive"),
    pos_rate = mean(R2_Result == "Positive")
  )

# group weeks: more than five weeks as the most recent category
prev_time %<>%
  mutate(
    week = ifelse(week >= 5, 5, week),
    week = as.character(week),
    week = ifelse(week == "5", "5+", week)
  ) %>%
  group_by(week) %>%
  # create new variables for number of tests and households per time period
  mutate(obs = n(),
         hh_per_week = length(unique(hh_id))) %>%
  ungroup()

# show tests and number of positives and negatives per time frame
prev_time_results <- prev_time %>%
  group_by(week) %>%
  summarise(
    n_hh = mean(hh_per_week),
    n_ind = mean(n()),
    pos_obs = sum(R2_Result == "Positive"),
    neg_obs = sum(R2_Result == "Negative"),
    pos_rate = mean(R2_Result == "Positive")
  ) %>%
  # crude standard error using normal distribution as approximation and
  # once considering too conservatively household clustering by using number of households in
  # standard deviation formula...
  mutate(
    se_hh = sqrt(pos_rate * (1 - pos_rate) / n_hh),
    lb_hh = pos_rate - 1.96 * se_hh,
    ub_hh = pos_rate + 1.96 * se_hh,
    # and once by ignoring household clustering
    se_ind = sqrt(pos_rate * (1 - pos_rate) / n_ind),
    lb_ind = pos_rate - 1.96 * se_ind,
    ub_ind = pos_rate + 1.96 * se_ind
  )
prev_time_results

# using the cluster bootstrap
prev_time <- prev_time %>%
  # select relevant variables
  dplyr::select(ind_id, hh_id, week, R2_Result) %>%
  # recode test result
  mutate(R2_Result = ifelse(R2_Result == "Positive", 1, 0)) %>%
  as.data.table()

# placeholder for results, one column for each time frame
pred_boot_clu <- matrix(0, boot, length(unique(prev_time$week)))
# cluster by household id
setkey(prev_time, "hh_id")
# store number of households
n_hh <- length(unique(prev_time$hh_id))
# set seed for replicability
set.seed(22)

# start loop
for (i in 1:boot)
{
  # create the bootstrap data set
  samp <- sample(unique(prev_time$hh_id), n_hh, replace = TRUE)
  bootdat <- prev_time[J(samp), allow.cartesian = TRUE]
  
  # calculate number of observations and prevalence for each time frame
  d_boot_grouped <- bootdat %>%
    group_by(week) %>%
    summarise(
      n_ind = mean(n()),
      pos_obs = sum(R2_Result == 1),
      pos_rate = pos_obs / n_ind
    )
  # save results
  pred_boot_clu[i, ] <- d_boot_grouped$pos_rate
}

# show results: point estimates for each time frame and 95% percentile intervals
apply(pred_boot_clu, 2, mean)
apply(pred_boot_clu, 2, function (x)
  quantile(x, probs = c(0.025, 0.975)))

# store results for unadjusted estimates
prev_time_results %<>%
  add_column(
    lb_boot = apply(pred_boot_clu, 2, function (x)
      quantile(x, probs = c(0.025))),
    ub_boot = apply(pred_boot_clu, 2, function (x)
      quantile(x, probs = c(0.975)))
  )  %>%
  # select relevant variables
  dplyr::select(week, n_ind, pos_obs, pos_rate, lb_boot, ub_boot) %>%
  # store results for adjusted estimates
  mutate(
    pos_rate_adj = (pos_rate + specificity - 1) / (specificity + sensitivity - 1),
    lb_boot_adj = (lb_boot + specificity - 1) / (specificity + sensitivity - 1),
    ub_boot_adj = (ub_boot + specificity - 1) / (specificity + sensitivity - 1)
  )
prev_time_results

# save results in an excel sheet
write_xlsx(
  prev_time_results,
  here_output("prevalence_over_time_unweighted_bootstrap.xlsx")
)

# plot results (not further commented)
plot <- prev_time_results %>%
  mutate(Week = as.numeric(str_replace(week, "\\+", ""))) %>%
  dplyr::select(-week) %>%
  gather(variable, value,-Week,-n_ind,-pos_obs) %>%
  mutate(value = value * 100) %>%
  mutate(Adjustment = case_when(
    variable %in% c("pos_rate", "lb_boot", "ub_boot") ~ "Unadjusted",
    variable %in% c("pos_rate_adj", "lb_boot_adj", "ub_boot_adj") ~ "Adjusted"
  )) %>%
  mutate(
    variable = case_when(
      variable %in% c("pos_rate", "pos_rate_adj") ~ "pos_rate",
      variable %in% c("lb_boot", "lb_boot_adj") ~ "lb_boot",
      variable %in% c("ub_boot", "ub_boot_adj") ~ "ub_boot"
    )
  ) %>%
  spread(variable, value) %>%
  ggplot(aes(x = Week)) +
  geom_ribbon(aes(ymin = lb_boot, ymax = ub_boot),
              alpha = 0.15,
              fill = "#0272b0") +
  geom_line(aes(y = pos_rate), col = "#0272b0") +
  geom_point(aes(y = pos_rate), col = "#0272b0") +
  facet_wrap( ~ Adjustment) +
  theme_bw() +
  ylab("Seroprevalence (%)") +
  scale_x_continuous(
    breaks = 1:5,
    labels = paste(
      c("1", "2", "3", "4", "5-11"),
      "\n(n=",
      prev_time_results$n_ind,
      ")",
      sep = ""
    )
  )
plot

# save plots
ggsave(
  here_output(file = "2nd-Round-Prevalence-Over-Time.pdf"),
  device = cairo_pdf,
  width = 10,
  height = 4
)

ggsave(
  here_output(file = "2nd-Round-Prevalence-Over-Time.png"),
  width = 10,
  height = 4
)
