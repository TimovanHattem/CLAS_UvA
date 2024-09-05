% ANALYSIS ACCURACY THETA TARGETING
%
% Timo van Hattem & Joao Patriota
% Experiment Theta_Targeting
% Updated: 12-3-2021

%% LOOP OVER ALL PARTICIPANTS TO GET ACCURACY ROSEPLOTS

% Set path
if isempty(strfind(path, ['/home/Public/Projects/sleep-and-cognition/generic/matlab-scripts/toolboxes/eeglab_v2020', pathsep]))
    addpath('/home/Public/Projects/sleep-and-cognition/generic/matlab-scripts/toolboxes/eeglab_v2020');
end
addpath(genpath('/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data'));
addpath('/home/Public/Projects/sleep-and-cognition/generic/matlab-scripts/toolboxes/CircStat2012a');
addpath('/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/matlab-scripts/2. ZMax Analysis/Algorithm Accuracy'); % to get e.g. 'bandpass.m'
fprintf('Paths added!\n')

%%

load('ACCURACY_INFORMATION_ALLPARTICIPANTS_CLEAN.mat')

pp = ['P' num2str(24)];

phase_total = [accuracy_info_allpp(71).PhasesStimMarkers accuracy_info_allpp(71).PhasesShamMarkers];
total_marker = [accuracy_info_allpp(71).TotalMarkers];

cmean = circ_mean(phase_total', [], 1);
if cmean <0
    cmean = (2*pi + cmean);
end

s=sprintf('%s',pp);
s0=sprintf('N (markers) = %.0f',length(phase_total));
s1=sprintf(['Mean = %.2f' char(176), ', SD = %.2f', char(176)],circ_rad2ang(cmean),circ_rad2ang(circ_std(phase_total')));

figure(1)
h = circ_plot(phase_total','hist',[], 20,false,true,'linewidth',2,'color','r');
title({s,s0,s1});

