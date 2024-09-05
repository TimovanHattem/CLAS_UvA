% ANALYSIS FOR AROUSALS THETA TARGETING
%
% Timo van Hattem & Joao Patriota
% Experiment Theta_Targeting
% Updated: 22-7-2021

%% SETTINGS

clear all
close all
clc

% Set paths
if isempty(strfind(path, ['/home/Public/Projects/sleep-and-cognition/generic/matlab-scripts/toolboxes/eeglab_v2020', pathsep]))
    addpath('/home/Public/Projects/sleep-and-cognition/generic/matlab-scripts/toolboxes/eeglab_v2020');
end
addpath(genpath('/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data'));
addpath('/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/matlab-scripts/2. ZMax Analysis/Algorithm Accuracy'); % to get e.g. 'bandpass.m'
addpath(genpath('/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results'));
fprintf('Paths added!\n')

z = 1; 

for i = 12:35
    if i == 16
        continue
    end 
    pp = ['P', num2str(i)];
    
    %% Load EEG dataset in eeglab
    eeglab;
    eegfilename = [pp, '_filtered.set'];
    eegpath = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data/ThetaTargeting_2021_Preprocessing/preproc/', pp];
    fprintf('Loading the .set file containing EEG recording for this participant...\n')
    EEG = pop_loadset('filename',eegfilename,'filepath',eegpath);
    fprintf('EEG file loaded!\n')
    
    %% Load arousal file
    arousalfile = ['thetatargeting', num2str(i), '_toscore_scores_th_arousals.csv'];
    arousalpath = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data/ThetaTargeting_2021_Preprocessing/sleep scoring/', pp];
    fprintf('Loading the .csv file containing sleep scoring for this participant...\n')
    arousal_data = readtable(arousalfile);
    arousal_data = table2struct(arousal_data);
    fprintf('Sleep scoring file loaded!\n')
        
    %% Determine arousals and blocks
    arousal_nonblock = 0;
    arousal_stim = 0;
    arousal_sham = 0;
    arousal_unidentified = 0;
    for j = 1:size(arousal_data,1)
        arousalstamp = arousal_data(j).StartTime*512;
        [a, index]  = min(abs([EEG.event.latency]-arousalstamp));
        
        if index == 1 || index == size(EEG.event,2)
            arousal_nonblock = arousal_nonblock + 1;
        elseif (strcmpi(EEG.event(index).type, '101') && strcmpi(EEG.event(index-1).type, '100')) || (strcmpi(EEG.event(index).type, '100') && strcmpi(EEG.event(index+1).type, '101'))
            arousal_nonblock = arousal_nonblock + 1;
            continue
        elseif strcmpi(EEG.event(index).type, '102') || (strcmpi(EEG.event(index).type, '101') && strcmpi(EEG.event(index+1).type, '102')) || (strcmpi(EEG.event(index).type, '101') && strcmpi(EEG.event(index+1).type, '104') && strcmpi(EEG.event(index-1).type, '102')) || (strcmpi(EEG.event(index).type, '101') && strcmpi(EEG.event(index+1).type, '105') && strcmpi(EEG.event(index-1).type, '102')) || (strcmpi(EEG.event(index).type, '105') && strcmpi(EEG.event(index-2).type, '102')) || (strcmpi(EEG.event(index).type, '105') && strcmpi(EEG.event(index-2).type, '104') && strcmpi(EEG.event(index-4).type, '102'))
            arousal_stim = arousal_stim+1;
        elseif strcmpi(EEG.event(index).type, '103') || (strcmpi(EEG.event(index).type, '101') && strcmpi(EEG.event(index+1).type, '103')) || (strcmpi(EEG.event(index).type, '101') && strcmpi(EEG.event(index+1).type, '104') && strcmpi(EEG.event(index-1).type, '103')) || (strcmpi(EEG.event(index).type, '101') && strcmpi(EEG.event(index+1).type, '105') && strcmpi(EEG.event(index-1).type, '103')) || (strcmpi(EEG.event(index).type, '105') && strcmpi(EEG.event(index-2).type, '103')) || (strcmpi(EEG.event(index).type, '105') && strcmpi(EEG.event(index-2).type, '104') && strcmpi(EEG.event(index-4).type, '103'))
            arousal_sham = arousal_sham+1;
        else
            arousal_unidentified = arousal_unidentified+1;
        end
    end
    
    arousal_info_allpp(z).Participant = pp;
    arousal_info_allpp(z).TotalArousals = size(arousal_data,1);
    arousal_info_allpp(z).ArousalStimulation = arousal_stim;
    arousal_info_allpp(z).ArousalSham = arousal_sham;
    arousal_info_allpp(z).ArousalNonblock = arousal_nonblock;
    arousal_info_allpp(z).ArousalUnidentified = arousal_unidentified;
    
    z = z + 1;
    
end
%save(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/Arousals/AROUSAL_INFORMATION_ALLPARTICIPANTS.mat'], 'arousal_info_allpp');

%% Statistics

close all

load('AROUSAL_INFORMATION_ALLPARTICIPANTS.mat');

ar_stim = [arousal_info_allpp.ArousalStimulation];
ar_sham = [arousal_info_allpp.ArousalSham];

[a,b] = vartest2(ar_stim,ar_sham);
[c,d] = jbtest(ar_stim);
[f,g] = jbtest(ar_sham);

[pt,ht,statst] = signrank(ar_stim, ar_sham);