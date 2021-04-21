# Responsible persons

* Ronan Le Gleut
* Turid Frahnow


# Content

The folder consists of the following scripts:

* Sampling_weights.R: The script computes the sampling weights for each household and individual according to the sampling design (all participants to the first round).
* Seroprevalence_second_round_calibrate_r1_test.R: The script computes the calibrated weights and calculate the variance associated to the sero-prevalence estimates for the second round: overall sero-prevalence at follow-up, sero-incidence (negative at baseline, positive at follow-up), sero-remission (positive at baseline, negative at follow-up), sero-persistence (positice at baseline, positive at follow-up). These estimates are adjusted for specificity and sensitivity of the test.
* Seroprevalence_second_round.R: Same as Seroprevalence_second_round_calibrate_r1_test.R but without the number of positive tests at baseline as a margin in the calibration process. This script is just used as a sensitivity analysis, therefore, no final estimate is computed.

# Requirements

Necessary data sets: 

* contiguity_matrices-munich-constituencies_200914.RData (for Sampling_weights.R): contiguity matrices of the constituencies
* Koco_baseline.csv (for Sampling_weights.R, Seroprevalence_second_round.R, Seroprevalence_second_round_calibrate_r1_test.R): all participants at baseline
* hh_characteristics_new.csv (for Seroprevalence_second_round.R, Seroprevalence_second_round_calibrate_r1_test.R): information on households at baseline
* ind_lab_baseline_new.csv (for Sampling_weights.R, Seroprevalence_second_round.R, Seroprevalence_second_round_calibrate_r1_test.R): all members in all households participating in the study (to identify households with children)
* KoCo19_Haushalte4Modeler_wRRstartConstituency_20200910.csv (for Sampling_weights.R, Seroprevalence_second_round.R, Seroprevalence_second_round_calibrate_r1_test.R): for each random route, starting and ending constituencies
* muc_hh_structure_KoCo_bis.csv (for Sampling_weights.R, Seroprevalence_second_round.R, Seroprevalence_second_round_calibrate_r1_test.R): total number of households per constituency in Munich
* muc_age_structure_KoCo_bis_study_pop.csv (for Seroprevalence_second_round.R, Seroprevalence_second_round_calibrate_r1_test.R): age and sex structure in Munich at 31.12.2020
* koco_r2_matched.csv (for Seroprevalence_second_round.R, Seroprevalence_second_round_calibrate_r1_test.R): participants to the follow-up (with lab results)
* specificity_sensitivity_classifier.RData (for Seroprevalence_second_round_calibrate_r1_test.R): sensitivity and specificity for Roche test


# Run

All scripts mentioned in the README are written with global paths.


# Results
Using the according scripts you get the following output.

* Sampling_weights.R -> data frame "KoCo_BLab" containing the sampling weights of the participants to the first round.
* Seroprevalence_second_round_calibrate_r1_test.R -> csv file "Estimates_Cal_R1_Test" containing the weighted seroprevalence estimates adjusted or not for specificity and sensitivity with the corresponding 95% CI. 
* Seroprevalence_second_round.R -> no real output, sensitivity analysis without including the number of tests at baseline as a margin in the calibration process.
