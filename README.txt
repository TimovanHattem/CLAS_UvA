FILE DESCRIPTIONS (by Timo van Hattem, 5-8-2021)

- 'old-irrelevant' stores MATLAB scripts that are old or irrelevant. That is, these scripts contain exploratory/interim/try-out code that was eventually not used anymore for the final analysis (e.g., comparisons of roseplots between the green and blue room; corrections in conversion of ECG electrodes; scripts for ERP/TFA analysis without loops, etc.). These scripts might still contain usable and valuable code for other projects. 

- 'Accuracy_averageplots_finalfigures_timo.m' is a MATLAB script to get the final average roseplots (with final lay-out etc.). This script uses the ACCURACY_INFORMATION_ALLPARPTICIPANTS_CLEAN.mat file as input.

- 'Accuracy_individualroseplots_finalfigures_timo.m' is a MATLAB script to get the individual roseplots of each participant (with final lay-out etc.). This script uses the ACCURACY_INFORMATION_ALLPARTICIPANTS_CLEAN.mat file as input.

- 'AccuracyScript_mar2021_LOOP2021_timo.m' is a MATLAB script to compute the algorithm performance (i.e., roseplots showing accuracy and precision) for each participant using three different filters. The output of this script is the file ACCURACY_INFORMATION_ALLPARTICIPANTS_CLEAN.mat. 

- 'ArousalScript_jul2021_timo.m' is a MATLAB script to get information on arousals during REM sleep. It uses sleep scoring file as input and gives the file AROUSAL_INFORMATION_ALLPARTICIPANTS.mat as output.

- 'CFWAnalysis_joao.m' is a MATLAB script to get ERP information of markers that are on-target or off-target. 

- 'ConversionScript_feb2021.m' is a MATLAB script that converts the raw EEG recording data to a .SET file.

- 'EpochScript_mar2021_LOOP_timo.m' is a MATLAB script that epochs the data for each participant around stimulus onset. Then this script checks if the markers are properly seperated (that is, more than two seconds apart). The output file is used for ERP and TFA analysis. 

- 'ERPAnalysis_mar2021_LOOP_timo.m' is a MATLAB script that performs ERP analysis for each participant and the total average. The output is the file ERP_INFORMATION_ALLPARTICIPANTS.mat. The script also computes the grand average ERP and gives the final figure as output. 

- 'permutest.m' is a MATLAB function that performs a cluster-based permutation test. This script is used for ERP and TFA statistical analysis.

- 'PreprocessingScript_mar2021_LOOP2_timo.m' is a MATLAB script that preprocesses the initial EEG data (e.g., filering, rereferencing). This script uses the sleep scoring file as input to extract only the REM signal of the EEG data. This scripts gives various files as output.

- 'stdshade.m' is a MATLAB function that adds the standard deviation or standard error as shade in a figure.

- 'TFAPSDAnalysis_mar2021_LOOP_timo.m' is a MATLAB script that performs power analysis for each participant and the total average. The output is the file TFA_INFORMATION_ALLPARTICIPANTS.mat.  

- 'TotalAverage_AccuracyScript_LOOP2_mar2021_timo.m' is a MATLAB script that computes the total average roseplot across all participants for three different filters. This script is relatively old - use 'Accuracy_averageroseplots_finalfigures_timo.m' for the final average roseplots instead. 

(For further information, please contact João Patriota or Timo van Hattem)