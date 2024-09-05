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

i = 2;
filt = 'FILT2';

total_phasevec_stim = [];
total_phasevec_sham = [];
for j = 1:size(accuracy_info_allpp,2)
    if strcmpi(accuracy_info_allpp(j).Filter, filt)
        total_phasevec_stim = [total_phasevec_stim, accuracy_info_allpp(j).PhasesStimMarkers];
        total_phasevec_sham = [total_phasevec_sham, accuracy_info_allpp(j).PhasesShamMarkers];
    end
end
totaltotal_phasevec = [total_phasevec_stim, total_phasevec_sham];

%% Statistics & Analysis
cmean_stim = circ_mean(total_phasevec_stim', [], 1);
if cmean_stim <0
    cmean_stim = (2*pi + cmean_stim);
end

[P_value_stim, Z_value_stim] = circ_rtest(total_phasevec_stim');
ave_stimmarkers = mean([accuracy_info_allpp.TotalStimMarkers],2);

cmean_sham = circ_mean(total_phasevec_sham', [], 1);
if cmean_sham <0
    cmean_sham = (2*pi + cmean_sham);
end

[P_value_sham, Z_value_sham] = circ_rtest(total_phasevec_sham');
ave_shammarkers = mean([accuracy_info_allpp.TotalShamMarkers],2);

cmean_total = circ_mean(totaltotal_phasevec', [], 1);
if cmean_total <0
    cmean_total = (2*pi + cmean_total);
end

[P_value_total, Z_value_total] = circ_rtest(totaltotal_phasevec');

[pval, table] = circ_wwtest(total_phasevec_stim', total_phasevec_sham');

%% Plot Figure Stimulation condition
s=sprintf(['Stimulation condition']);
s1=sprintf(['(M = %.2f' char(176), ', SD = %.2f)'],circ_rad2ang(cmean_stim),circ_rad2ang(circ_std(total_phasevec_stim')));
s0=sprintf('');

figure(1)
polarhistogram(total_phasevec_stim',20,'FaceColor','blue','FaceAlpha',.3);
title({s,s1,s0}, 'FontSize', 20, 'FontWeight', 'bold');
%circ_plot(total_phasevec_stim','hist',[], 20,false,true,'linewidth',2,'color','r');
    
%% Plot figure Sham condition
s=sprintf('Sham condition');
s1=sprintf(['(M = %.2f' char(176), ', SD = %.2f)'],circ_rad2ang(cmean_sham),circ_rad2ang(circ_std(total_phasevec_sham')));
s0=sprintf('');

figure(2)
polarhistogram(total_phasevec_sham', 20,'FaceColor','red','FaceAlpha',.3);
title({s,s1,s0},'FontSize', 20, 'FontWeight', 'bold');
%circ_plot(total_phasevec_sham','hist',[], 20,false,true,'linewidth',2,'color','r');
    
%% Plot figure total average
s=sprintf('TOTAL (Stimulation + Sham)');
s1=sprintf(['(M = %.2f' char(176), ', SD = %.2f)'],circ_rad2ang(cmean_total),circ_rad2ang(circ_std(totaltotal_phasevec')));
s2=sprintf('');

figure(3)
polarhistogram(totaltotal_phasevec', 20,'FaceColor','black','FaceAlpha',.3);
title({s,s1,s2},'FontSize', 20, 'FontWeight', 'bold');
%circ_plot(totaltotal_phasevec','hist',[], 20,false,true,'linewidth',2,'color','r');
    

circ_rad2ang(cmean_total)
circ_rad2ang(circ_std(totaltotal_phasevec'))

%%

% hoi(1,:) = unique([accuracy_info_allpp.TotalStimMarkers]);
% hoi(2,:) = unique([accuracy_info_allpp.TotalShamMarkers]);

[a,p1] = vartest2(hoi(1,:),hoi(2,:));
[b, p2] = jbtest(hoi(1,:));
[c, p3] = jbtest(hoi(2,:));

[d, p4, ci, stats] = ttest(hoi(1,:),hoi(2,:))