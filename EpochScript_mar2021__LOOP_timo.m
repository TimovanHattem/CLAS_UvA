% ANALYSIS EPOCH FOR STA/TFA/PSD THETA TARGETING
%
% Timo van Hattem & Joao Patriota
% Experiment Theta_Targeting
% Updated: 1-6-2021

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPROCESS DATA FOR ERP/TFA ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set pathS
if isempty(strfind(path, ['/home/Public/Projects/sleep-and-cognition/generic/matlab-scripts/toolboxes/eeglab_v2020', pathsep]))
    addpath('/home/Public/Projects/sleep-and-cognition/generic/matlab-scripts/toolboxes/eeglab_v2020');
end
addpath(genpath('/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data'));
addpath('/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/matlab-scripts/2. ZMax Analysis/Algorithm Accuracy'); % to get e.g. 'bandpass.m'
addpath(genpath('/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results'));
fprintf('Paths added!\n')

%% Loop over all participants

for i = [18,23,32]

    % Clear workspace
    close all
    clc

    % Load EEG dataset in eeglab
    eeglab;
    eegfilename = ['P', num2str(i), '_continuousonlyREMonlyFz.set'];
    eegpath = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data/ThetaTargeting_2021_Preprocessing/preproc/P', num2str(i)];
    EEG = pop_loadset('filename',eegfilename,'filepath',eegpath);
    fprintf('EEG file loaded!\n')

    % Setting output file path
    pp = EEG.part_num;
    if strcmpi(pp,'P66')
        pp = 'P16';
    end 
    mkdir(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/'], pp);
    DATAOUT = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/',pp];

    % Filtering
    %theta_rang = [ 3 30 ];
    %EEG.data = bandpass(EEG.data,1,theta_rang(1), theta_rang(2),EEG.srate);
    %EEG = pop_eegfiltnew(EEG, 4, 29, [], 0);
    
    % Removing irrelevant events
    index_nonevent = [];
    for i = 1:size(EEG.event,2)
        if strcmpi(EEG.event(i).type, '101') || strcmpi(EEG.event(i).type,'104') || strcmpi(EEG.event(i).type,'105') || strcmpi(EEG.event(i).type, 'X') || strcmpi(EEG.event(i).type, 'boundary') || strcmpi(EEG.event(i).type, '100') || strcmpi(EEG.event(i).type, '69')
            index_nonevent = [index_nonevent i];
        end
    end
    EEG.event(index_nonevent) = [];

    % Detect properly spaced stimuli
    remov = [];
    for i = 1:size(EEG.event,2)
        for z = i+1:size(EEG.event,2)
            if EEG.event(z).latency-EEG.event(i).latency <= 768
                remov = [remov z];
                if EEG.event(z).latency-EEG.event(i).latency <= 512
                    remov = [remov i];
                end
            end
        end
    end
    remov = unique(remov);

    % Epoching
    [EEG, indx2] = pop_epoch(EEG, {'102', '103'}, [-0.5,1]);

    % Baseline correction
    EEG = pop_rmbase(EEG, [-500 0]);

    % Sorting trials per condition
    EEG = pop_select(EEG, 'notrial', remov);

    stim_epochs = [];
    sham_epochs = [];
    for i = 1:size(EEG.event,2)
        if strcmpi(EEG.event(i).type, '102')
            stim_epochs = [stim_epochs i];
        else
            sham_epochs = [sham_epochs i];
        end
    end

    EEG_stim = pop_select(EEG, 'notrial', sham_epochs);
    EEG_sham = pop_select(EEG, 'notrial', stim_epochs);

    % Save trials for STA/TFA/PSD
    pop_saveset(EEG_stim, 'filename', [EEG.part_num, '_epoched_STIM.set'], 'filepath', DATAOUT);
    pop_saveset(EEG_sham, 'filename', [EEG.part_num, '_epoched_SHAM.set'], 'filepath', DATAOUT);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Artefact rejection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eeglab; % Use GUI to reject trials and save file

