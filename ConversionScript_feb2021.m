% CONVERTING EEG DATA
%
% Complete script to open  tmsi poly5 file and save as a .set file AND a .set
% file to be sleep scored.
%
% Timo van Hattem & Joao Patriota
% Experiment ThetaTargeting
% Updated: 6-5-2021

%%

clear all
close all
clc
%function [] = tmsi2setjp(filename, folder)
filename = input('Write the filename (e.g. thetatargeting01): ', 's');
pp = input('Write the participant number (e.g. P01): ', 's'); 

addpath('/home/Public/Projects/sleep-and-cognition/generic/matlab-scripts/toolboxes/eeglab_v2020');
addpath('/home/Public/Projects/sleep-and-cognition/experiments/0010.dreams-SO-memory/matlab-scripts/1-data_format/');
addpath(genpath('/home/Public/Projects/sleep-and-cognition/experiments/0010.dreams-SO-memory/raw-data/'));
load('chanlocs_68ch.mat');
eeglab;

%DATAIN = ['/home/Public/Projects/sleep-and-cognition/experiments/009.simulation-theta/raw-data/EEG/', folder, '/tmsi/'];
%DATAOUT = ['/home/Public/Projects/sleep-and-cognition/exDATAIN = ['/home/Public/Projects/sleep-and-cognition/experiments/009.simulation-theta/raw-data/EEG/', folder, '/tmsi/'];

mkdir(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data/ThetaTargeting_2021_Preprocessing/set'], pp);
DATAIN = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data/ThetaTargeting_2021_Raw/'];
DATAOUT = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data/ThetaTargeting_2021_Preprocessing/set/', pp];

%cd(DATAIN)
[file, path] = uigetfile('*.Poly5'); % select TMSi file to convert

%  files
TMSi = tms_read(file, path); % read TMSi file

% Store Channels into dataframe
Channels = 68; % 68 channels to include ECG
for i=1:Channels
    EEG.data(i,:) = TMSi.data{i,1};
end

%% Import events
events = strcat(TMSi.path,TMSi.fname,'.events.csv');
triggers = csvimport(events);

for k=2:1:size(triggers,1)
  time=strsplit(triggers{k,1},':');
  sample=str2double(time{3})*TMSi.fs + 60*str2double(time{2})*TMSi.fs + 3600*str2double(time{1})*TMSi.fs;
  EEG.event(1,k-1).type = strrep(triggers{k,3},'"','');
  EEG.event(1,k-1).latency = sample;
  EEG.event(1,k-1).decomposed_time = time;
  EEG.event(1,k-1).time = strjoin(time,':'); % time; EJ changed to see the actual time
end
%% Convert to eeglab
% Add additional info to set file
EEG.srate = TMSi.fs;
EEG.setname = [filename, '.set']; % TMSi.fname;
EEG.filename = [DATAOUT, EEG.setname];% TMSi.filename;
EEG.filepath = DATAOUT;
EEG.chanlocs = chanlocs_68ch; % or chanlocs_68ch
EEG.recording_time = {'date', 'start_time', 'duration'; TMSi.measurementdate, TMSi.measurementtime, TMSi.measurementduration};

eeglab checkset % correct for consistency within eeglab

%% Taking channels out and preparing to export for sleepscoring
%ChannelsToKeep={'Fpz' 'F4' 'Cz' 'C4' 'Pz' 'Oz' 'O2' 'HEOG' 'VEOG' 'EMG'};
% On our current setup (feb 2021) we know that those channels belong to

indexes = 1:68; %a vector with all possible indexes
Referenced_EEG = pop_reref(EEG,[13 19],'keepref', 'on', 'exclude',[65 66 67 68]); %re-referencing data
EEG_tosleep = Referenced_EEG; % this script will emit 2 files: one for scoring and another for regular analysis

channel_number = [2 7 16 17 26 31 32 65 66 67]; %channels we want to keep
idx_tokeep = ismember(indexes, channel_number); %indexes of the channels we want to keep
EEG_datatoosleep = EEG_tosleep.data(idx_tokeep',:); %extracting relevant channels
EEG_channelstoosleep = EEG_tosleep.chanlocs(idx_tokeep); %extracting relevant channels

EEG_tosleep.data = EEG_datatoosleep;
EEG_tosleep.chanlocs = EEG_channelstoosleep;
%% Filtering data,
EEG_tosleep.nbchan = 10; %this is needed to inform that there are only 10 channels on the current structure. 
EEG_toscore_filtered = pop_eegfiltnew(EEG_tosleep, 0.5, 30); %
%% Resampling
EEG_toscore_filtered = pop_resample(EEG_toscore_filtered,256); %re-sampling the file to be scored
%% Saving both .set files (complete and to sleepscore)
filename_toscore = [filename, '_toscore']; %making the filename
pop_saveset(Referenced_EEG,filename, EEG_toscore_filtered.filepath) %saving re-referenced full .set file
pop_saveset(EEG_toscore_filtered,filename_toscore, EEG_toscore_filtered.filepath) %saving compact .set file (to sleepscore)