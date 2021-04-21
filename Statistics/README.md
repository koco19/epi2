# Responsible persons

* Peter Pütz
* Mercè Garí


# Content

The folder consists of the following scripts:

* seroprevalence_unweighted_over_time.R: This script calculates the unweighted point estimates and associated confidence intervals (unadjusted and adjusted for specificity and sensitivity for overall and weekly prevalences for the second round of visit of the KoCo19 project.
* 2nd-Round-Study-Framework.R: This script takes the LMU nowcast data (available at https://corona.stat.uni-muenchen.de/nowcast/) to plot the 2nd Round Study Framework (Figure 1).
* Plot-Prevalence-Incidence.R: This scripts takes the prevalence/incidence data calculated in the Sampling Weights to plot Figure 4.
* Plot-OR.R: This script takes the model output data to plot Figure 5.
* epi_maps-prevalence.R: This script takes the Density population and seroprevalence estimates of the Munich districts to plot the Prevalence Map.


# Requirements

Necessary data sets: 

* koco_r2_dat.csv: Test Results for Koco19 participants for the first two rounds.
* results_nowcast_2021-03-16.csv: Number of Cases per day/month in Munich.
* Estimates_Cal_R1_Test.csv: Estimates for Prevalence and Incidence for the 2nd round.
* OR_R2Positive.xlsx: OR estimates (and 95% credible intervals).
* Muc_popn_density_Wiki.xlsx and District_Seroprev.xlsx (Density population and seroprevalence of the Munich districts).


# Run

All scripts mentioned in the README are written with global paths.


# Results
Using the according script you get the following output which can be found in the folder "Output".

* Sampling_weights.R -> 
  - prevalence_unweighted_bootstrap.xlsx: Contains the unweighted prevalence estimates (unadjusted and adjusted for   specificity and sensitivity) with associated 95% cluster bootstrap percentile intervals.
  - prevalence_over_time_unweighted_bootstrap.xlsx: Contains the weekly unweighted prevalence estimates (unadjusted and adjusted for specificity and sensitivity) with associated 95% cluster bootstrap percentile intervals.
  - 2nd-Round-Prevalence-Over-Time.pdf / 2nd-Round-Prevalence-Over-Time.png: Plots the weekly unweighted prevalence estimates (unadjusted and adjusted for specificity and sensitivity) with associated 95% cluster bootstrap percentile intervals.
  - Figure_1.pdf: Figure 1 of the manuscript regarding the Study Framework of the KoCo19 2nd round.
  - Figure_4.pdf: Figure 4 of the manuscript regarding the Prevalence and Incidence of the KoCo19 2nd round.
  - Figure_5.pdf: Figure 5 of the manuscript regarding the OR estimates (and 95% credible intervals).
  - Epi2-PrevMap-Bin-Global-NA-Weighted.pdf: Map of prevalences in the city of Munich by districts.
  
