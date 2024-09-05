% PREPROCESSING EEG DATA
%
% This script uses already re-referenced data (see
% ConversionScript_feb2021.m). This script filters data, removes unused
% channels, selects REM epochs and does trial rejection on bad epochs based
% on visual inspection.
% The output is a set file with clean REM data, ready for further analyses.
%
% This script is partly an adaptation of a script bij VPathak & EJuan
%
% Timo van Hattem & Joao Patriota
% Experiment ThetaTargeting
% Updated: 4-2-2021

%% Set path
if isempty(strfind(path, ['/home/Public/Projects/sleep-and-cognition/generic/matlab-scripts/toolboxes/eeglab_v2020', pathsep]))
    addpath('/home/Public/Projects/sleep-and-cognition/generic/matlab-scripts/toolboxes/eeglab_v2020');
end
addpath(genpath('/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data'));
addpath('/home/Public/Projects/sleep-and-cognition/generic/matlab-scripts/toolboxes/CircStat2012a');
fprintf('Paths added!\n')

for i = 12:35
    
    if i == 16
        continue
    end
    
    %     %% Setting output file path
    pp = ['P', num2str(i)];
    %     mkdir(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data/ThetaTargeting_2021_Preprocessing/preproc'], pp);
    DATAOUT = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data/ThetaTargeting_2021_Preprocessing/preproc/', pp];
    
    %% Load EEG dataset in eeglab
    eeglab;
    eegfilename = [pp, '_filtered.set'];
    eegpath = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data/ThetaTargeting_2021_Preprocessing/preproc/', pp];
    fprintf('Loading the .set file containing EEG recording for this participant...\n')
    EEG = pop_loadset('filename',eegfilename,'filepath',eegpath);
    fprintf('EEG file loaded!\n')
    
    %% Load sleep stage file
    sleepscorefile = ['thetatargeting', num2str(i), '_toscore_scores_th.csv'];
    sleepscorepath = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data/ThetaTargeting_2021_Preprocessing/sleep scoring/', pp];
    fprintf('Loading the .csv file containing sleep scoring for this participant...\n')
    sleep_data = readtable(sleepscorefile);
    sleep_data = table2struct(sleep_data);
    fprintf('Sleep scoring file loaded!\n')
    %
    %     %% Get participant info
    %     EEG.part_num = pp;
    %
    %     %% Input check datasets
    %     %Get sampling rate data check if it is 512
    %     Fs = EEG.srate;
    %     if Fs ~= 512
    %         fprintf('Please check dataset, sampling rate is not 512Hz')
    %         fprintf('Press any button to continue processing')
    %         pause
    %     end
    %
    %     % Epoch length that will be used in seconds
    %     Epochlength = 30;
    %
    %     % Determine the length of the data in seconds
    %     LengthDataSeconds = length(EEG.data)/Fs;
    %
    %     % Determine amount of samples in 1 epoch
    %     EpochSrate = Fs*Epochlength;
    %
    %     % Determine number of epochs based on length of data and amount of samples
    %     % in epoch
    %     xepochsEEG = floor(length(EEG.data)/EpochSrate);
    %
    %     %Determine amount of epochs csv file
    %     sleep_data(size(sleep_data,1)) = [];
    %     xepochsscore= size(sleep_data,1);
    %
    %     % Check if number of epochs matches the data and print result of check.
    %     if xepochsEEG ~= xepochsscore
    %         difference = xepochsEEG - xepochsscore;
    %         fprintf('Difference between EEG epochs and sleep score epochs is: %f\n', difference)
    %         error('Number of epochs in EEG data does not match number of epochs in the sleep score file')
    %     end
    %     fprintf('Equal amount of epochs in EEG and sleep score file :) \n')
    %
    %     %% Remove non EEG channels
    %     EEG = pop_select(EEG, 'nochannel', [65:68]);
    %
    %     %% Filter data
    %     EEG_filtered = pop_eegfiltnew(EEG, 0.15, 100, [], 0); %bandpass filter 0.15-100
    %     EEG_filtered = pop_eegfiltnew(EEG_filtered, 48, 52, [], 1); %notch filter 48-52
    %     pop_saveset(EEG_filtered,'filename',[EEG_filtered.part_num,'_filtered.set'],'filepath',DATAOUT);
    %     %pop_eegplot(EEG_filtered,1,1,1);
    %
    %% Rereferencing
    EEG_filtered = EEG;
    EEG_filtered = pop_reref(EEG_filtered,[13 19],'keepref', 'on'); %to have the same chanloc struct.
    pop_saveset(EEG_filtered,'filename',[EEG_filtered.part_num,'_rereferenced.set'],'filepath',DATAOUT);
    
    %% Divide EEG data into 30 seconds epochs
    EEG_epoched =  eeg_regepochs(EEG_filtered,'recurrence',[30],'limits',[0 30]);
    pop_saveset(EEG_epoched,'filename',[EEG_filtered.part_num,'_epoched.set'],'filepath',DATAOUT);
    
    %% Select REM epochs
    REM_epochs = [];
    for i = 1:size(sleep_data, 1)
        if strcmpi(sleep_data(i).stage, 'REM')
            REM_epochs = [REM_epochs i]; % Vector with indices of REM epochs
        end
    end
    EEG_REMepochs = pop_select(EEG_epoched, 'trial', REM_epochs);
    %pop_eegplot(EEG_REMepochs,1,1,1);
    
    epvector = 1:size(sleep_data,1);
    NREM_epochs = setdiff(epvector, REM_epochs); % Vector with indices of NREM epochs
    EEG_REMepochs.REMepochs = REM_epochs;
    EEG_REMepochs.NREMepochs = NREM_epochs;
    EEG_REMepochs.orgevents = EEG_epoched.event;
    pop_saveset(EEG_REMepochs,'filename',[EEG_filtered.part_num,'_epochedonlyREMallchannel.set'],'filepath',DATAOUT);
    
    EEGout = epoch2continuous(EEG_REMepochs);
    pop_saveset(EEGout,'filename',[EEG_filtered.part_num,'_continuousonlyREMallchannel.set'],'filepath',DATAOUT);
    
    
    %     %% Select and extract channel of interest to see REM epochs
    %     coi = 6; %%%%%%%% INPUT channel of interest, Fz = 6, Fpz = 2 %%%%%%%%
    %     EEG_newchannels = pop_select(EEG_REMepochs, 'channel', coi);
    %     %pop_eegplot(EEG_newchannels,1,1,1);
    %
    %     EEG_newchannels.REMepochs = REM_epochs;
    %     EEG_newchannels.NREMepochs = NREM_epochs;
    %     EEG_newchannels.orgevents = EEG_epoched.event;
    %     pop_saveset(EEG_newchannels,'filename',[EEG_newchannels.part_num,'_tobeinspected_', EEG_newchannels.chanlocs.labels, '.set'],'filepath',DATAOUT);
    %
    %
    close all
    clc
    clear all
    %
end

% %% MANUALLY REJECTING TRIALS AFTER LOOP
%
% %% Load data
% pp = input('Write the participant number (e.g. P01): ', 's');
% DATAOUT = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data/ThetaTargeting_2021_Preprocessing/preproc/', pp];
%
% eeglab;
% [eegfilename eegpath]= uigetfile('*.set'); % LOAD TOBEINSPECTED.SET
% fprintf('Loading the .set file containing EEG recording for this participant...\n')
% EEG_newchannels = pop_loadset('filename',eegfilename,'filepath',eegpath);
% fprintf('EEG file loaded!\n')
%
% eeglab;
% [eegfilename eegpath]= uigetfile('*.set'); % LOAD FILTERED.SET
% fprintf('Loading the .set file containing EEG recording for this participant...\n')
% EEG_filtered = pop_loadset('filename',eegfilename,'filepath',eegpath);
% fprintf('EEG file loaded!\n')
%
% [sleepscorefile, sleepscorepath] = uigetfile('*.csv');
% fprintf('Loading the .csv file containing sleep scoring for this participant...\n')
% sleep_data = readtable(sleepscorefile);
% sleep_data = table2struct(sleep_data);
% fprintf('Sleep scoring file loaded!\n')
% sleep_data(size(sleep_data,1)) = [];
%
% %% Inspect trials
% pop_eegplot(EEG_newchannels,1,1,1); % Settings > Time range = 1 epoch, range = 200
%
% %% Reject trials after visual inspection
% rejtrial = [13,18,25,34,53,59,91,103]; %%%%%%%% INPUT give vector with numbers for trials to reject %%%%%%%%%%
% save([DATAOUT, '/reject', EEG_newchannels.part_num, '.mat'], 'rejtrial')
%
% indexes = 1:size(EEG_newchannels.data,3);
% idx_keep = ismember(indexes, rejtrial); %indexes of the trials we want to reject
% EEG_reject = pop_rejepoch(EEG_newchannels, idx_keep);
%
% %% Concatenate all epochs again for further analysis
% EEGout = epoch2continuous(EEG_reject);
%
% %% Save preprocessed EEG data
% EEGout.rejecttrial = rejtrial;
% pop_saveset(EEGout,'filename',[EEGout.part_num,'_preprocessed_', EEGout.chanlocs.labels, '.set'],'filepath',DATAOUT);
%
% %% Create and save reference EEG dataset with all markers %%
% REM_epochs = EEG_newchannels.REMepochs;
% NREM_epochs = EEG_newchannels.NREMepochs;
%
% badtrials = REM_epochs(rejtrial);
% delmat = sort(cat(2,NREM_epochs,badtrials));
%
% EEG_ref = EEG_filtered;
% delete_epochs = [];
% for i = 1:length(delmat)
%     delete_epochs(i,1) = sleep_data(delmat(i)).start*512+1;
%     delete_epochs(i,2) = sleep_data(delmat(i)).xEnd*512+1;
%     EEG_ref.data(:,delete_epochs(i,1):delete_epochs(i,2)) = [-300];
% end
%
% coi = 6; %%%%%%%% INPUT channel of interest, Fz = 6, Fpz = 2 %%%%%%%%
% EEG_ref = pop_select(EEG_ref, 'channel', coi); %using this variable to save filtered, but with important epochs only
% pop_saveset(EEG_ref,'filename',[EEG_ref.part_num,'_refEEG_', EEG_ref.chanlocs.labels, '.set'],'filepath',DATAOUT);
%
% %%
%
% x = 0;
% for i = 1:size(EEG.event,2)
%     if strcmpi(EEG.event(i).type,'102')|| strcmpi(EEG.event(i).type,'103')
%         x = x + 1;
%     end
% end


