% ANALYSIS ACCURACY THETA TARGETING
%
% Timo van Hattem & Joao Patriota
% Experiment Theta_Targeting
% Updated: 12-3-2021

%% Settings
% Clear workspace
clear all
close all
clc

% Set path
if isempty(strfind(path, ['/home/Public/Projects/sleep-and-cognition/generic/matlab-scripts/toolboxes/eeglab_v2020', pathsep]))
    addpath('/home/Public/Projects/sleep-and-cognition/generic/matlab-scripts/toolboxes/eeglab_v2020');
end
addpath(genpath('/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data'));
addpath('/home/Public/Projects/sleep-and-cognition/generic/matlab-scripts/toolboxes/CircStat2012a');
addpath('/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/matlab-scripts/2. ZMax Analysis/Algorithm Accuracy'); % to get e.g. 'bandpass.m'
addpath(genpath('/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results'));
fprintf('Paths added!\n')

% Setting output file path
DATAOUT = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/Algorithm performance/Fz/Clean data roseplots/TOTAL AVERAGE'];

%% Load data participants
load('ACCURACY_INFORMATION_ALLPARTICIPANTS_CLEAN.mat');

for i = 1:3
    if i == 1
        filt = 'FILT1';
    elseif i == 2
        filt = 'FILT2';
    elseif i == 3
        filt = 'FILT3';
    end
    
    total_phasevec_stim = [];
    total_phasevec_sham = [];
    for j = 1:size(accuracy_info_allpp,2)
        if strcmpi(accuracy_info_allpp(j).Filter, filt)
            total_phasevec_stim = [total_phasevec_stim, accuracy_info_allpp(j).PhasesStimMarkers];
            total_phasevec_sham = [total_phasevec_sham, accuracy_info_allpp(j).PhasesShamMarkers];
        end
    end
    
    %% Statistics & Analysis
    cmean_stim = circ_mean(total_phasevec_stim', [], 1);
    if cmean_stim <0
        cmean_stim = (2*pi + cmean_stim);
    end
    
    cmedian_stim = circ_median(total_phasevec_stim');
    if cmedian_stim<0
        cmedian_stim= (2*pi + cmedian_stim);
    end
    
    [P_value_stim, Z_value_stim] = circ_rtest(total_phasevec_stim');
    ave_stimmarkers = mean([accuracy_info_allpp.TotalStimMarkers],2);
    ave_predfreqstim = mean([accuracy_info_allpp.PredictionFrequencyStimMarkers],2);
    
    cmean_sham = circ_mean(total_phasevec_sham', [], 1);
    if cmean_sham <0
        cmean_sham = (2*pi + cmean_sham);
    end
    
    cmedian_sham = circ_median(total_phasevec_sham');
    if cmedian_sham<0
        cmedian_sham= (2*pi + cmedian_sham);
    end
    
    [P_value_sham, Z_value_sham] = circ_rtest(total_phasevec_sham');
    ave_shammarkers = mean([accuracy_info_allpp.TotalShamMarkers],2);
    ave_predfreqsham = mean([accuracy_info_allpp.PredictionFrequencyShamMarkers],2);
    
    
    %% Plot Figure Stimulation condition
    s=sprintf('Average All Participants: Stimulation condition');
    s0=sprintf('Average number of markers = \t%.0f',ave_stimmarkers);
    s1=sprintf(['Mean phase = \t%.2f ' char(176) 'C'],circ_rad2ang(cmean_stim));
    s2=sprintf(['Median phase = \t%.2f ' char(176) 'C'],circ_rad2ang(cmedian_stim));
    s3=sprintf('SD = \t%.2f',circ_rad2ang(circ_std(total_phasevec_stim')));
    s4=sprintf('SE = \t%.2f',circ_rad2ang(circ_std(total_phasevec_stim')/(sqrt(length(total_phasevec_stim')))));
    %s5=sprintf('Skewness (0=symm) = \t%.2f',circ_skewness(total_phasevec'));
    %s6=sprintf('Kurtosis (0=norm.d.) = \t%.2f',circ_kurtosis(total_phasevec'));
    %s7 = sprintf('Z-value = \t%.5f,\t P-value = \t%.5f', Z_value, P_value);
    s8 = sprintf('Average prediction frequency (total REM) = \t%.3f Hz', ave_predfreqstim);
    
    figure
    circ_plot(total_phasevec_stim','hist',[], 20,false,true,'linewidth',2,'color','r');
    title({s,s0,s8,s1,s2,s3,s4});
    set(gcf, 'MenuBar', 'none')
    set(gcf, 'Toolbar', 'none')
    set(gcf, 'Position', get(0, 'Screensize'));
    
    figname = [DATAOUT, '/accuracy_roseplot_AVERAGE_STIM_', filt, '.jpg'];
    saveas (gcf,figname);
    
    %% Plot figure Sham condition
    s=sprintf('Average All Participants: Sham condition');
    s0=sprintf('Average number of markers = \t%.0f',ave_shammarkers);
    s1=sprintf(['Mean phase = \t%.2f ' char(176) 'C'],circ_rad2ang(cmean_sham));
    s2=sprintf(['Median phase = \t%.2f ' char(176) 'C'],circ_rad2ang(cmedian_sham));
    s3=sprintf('SD = \t%.2f',circ_rad2ang(circ_std(total_phasevec_sham')));
    s4=sprintf('SE = \t%.2f',circ_rad2ang(circ_std(total_phasevec_sham')/(sqrt(length(total_phasevec_sham')))));
    %s5=sprintf('Skewness (0=symm) = \t%.2f',circ_skewness(total_phasevec'));
    %s6=sprintf('Kurtosis (0=norm.d.) = \t%.2f',circ_kurtosis(total_phasevec'));
    %s7 = sprintf('Z-value = \t%.5f,\t P-value = \t%.5f', Z_value, P_value);
    s8 = sprintf('Average prediction frequency (total REM) = \t%.3f Hz', ave_predfreqsham);
    
    figure
    circ_plot(total_phasevec_sham','hist',[], 20,false,true,'linewidth',2,'color','r');
    title({s,s0,s8,s1,s2,s3,s4});
    set(gcf, 'MenuBar', 'none')
    set(gcf, 'Toolbar', 'none')
    set(gcf, 'Position', get(0, 'Screensize'));
    
    figname = [DATAOUT, '/accuracy_roseplot_AVERAGE_SHAM_', filt, '.jpg'];
    saveas (gcf,figname);
    
    close all
end
