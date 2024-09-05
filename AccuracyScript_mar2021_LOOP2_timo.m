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

k = 1;

for i = 12:35
    
    %% Settings
    
    % Clear workspace
    close all
    clc
    
    %% Load EEG dataset in eeglab
    eeglab;
    eegfilename = ['P', num2str(i), '_preprocessed.set'];
    eegpath = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data/ThetaTargeting_2021_Preprocessing/preproc/P', num2str(i)];
    fprintf('Loading the .set file containing EEG recording for this participant...\n')
    EEG = pop_loadset('filename',eegfilename,'filepath',eegpath);
    fprintf('EEG file loaded!\n')
    
    %% Setting output file path
    pp = EEG.part_num;
    mkdir(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/Algorithm performance/Fz/Clean data roseplots'], pp);
    DATAOUT = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/Algorithm performance/Fz/Clean data roseplots/',pp];
    
    for j = 1:3
        
        %% Apply filter on theta range [4-8 Hz]
        if j == 1
            filt = 'FILT1';
            low_theta = [5];
            high_theta = [7];
            EEG_thetafilt = pop_eegfiltnew(EEG,low_theta,high_theta,[],0);
            %         elseif j == 2
            %             filt = 'FILT2';
            %             [B,A] = butter(1,[(4/(512/2)), (8/(512/2))], 'bandpass');
            %             EEG_thetafilt = filtfilt(B,A,EEG.data(6,:));
        elseif j == 2
            filt = 'FILT2';
            ord = 1; %Filter order
            Fs = EEG.srate; % 512; %Sampling rate in Hz 512 is the refa amplifier
            theta_rang = [4 8];
            EEG_thetafilt = bandpass(EEG.data,ord,theta_rang(1), theta_rang(2),Fs);
        elseif j == 3
            filt = 'FILT3';
            ord = 2; %Filter order
            Fs = EEG.srate; % 512; %Sampling rate in Hz 512 is the refa amplifier
            theta_rang = [4 8];
            EEG_thetafilt = bandpass(EEG.data,ord,theta_rang(1), theta_rang(2),Fs);
        end
        
        %% Calculate Hilbert transform
        if j == 1
            hil_eeg = hilbert(EEG_thetafilt.data);
        else
            hil_eeg = hilbert(EEG_thetafilt);
        end
        
        %% STIMULATION condition %%
        
        %% Find markers STIMULATION condition
        stimtrials = [];
        for i = 1:size(EEG.event,2)
            if strcmpi(EEG.event(i).type,'102')
                stimtrials = [stimtrials, EEG.event(i).latency];
            end
        end
        
        %% Find prediction phases STIMULATION condition
        hil_stim = hil_eeg(stimtrials);
        phase_stim = angle(hil_stim)+(pi/2);
        
        %% Statistics & Analysis STIMULATION condition
        cmean_stim = circ_mean(phase_stim', [], 1);
        if cmean_stim <0
            cmean_stim = (2*pi + cmean_stim);
        end
        
        cmedian_stim = circ_median(phase_stim');
        if cmedian_stim<0
            cmedian_stim= (2*pi + cmedian_stim);
        end
        
        [P_value_stim, Z_value_stim] = circ_rtest(phase_stim');
        
        totalremtime = length(hil_eeg)/512;
        predfreq_stim = length(stimtrials)/totalremtime;
        
        %% Plot roseplot STIMULATION condition
        s=sprintf('%s: Stimulation condition',EEG.part_num);
        s0=sprintf('Total number of markers = \t%.0f',length(phase_stim));
        s1=sprintf(['Mean phase = \t%.2f ' char(176) 'C'],circ_rad2ang(cmean_stim));
        s2=sprintf(['Median phase = \t%.2f ' char(176) 'C'],circ_rad2ang(cmedian_stim));
        s3=sprintf('SD = \t%.2f',circ_rad2ang(circ_std(phase_stim')));
        s4=sprintf('SE = \t%.2f',circ_rad2ang(circ_std(phase_stim')/(sqrt(length(phase_stim')))));
        %s5=sprintf('Skewness (0=symm) = \t%.2f',circ_skewness(phase_stim'));
        %s6=sprintf('Kurtosis (0=norm.d.) = \t%.2f',circ_kurtosis(phase_stim'));
        %s7 = sprintf('Z-value = \t%.5f,\t P-value = \t%.5f', Z_value_stim, P_value_stim);
        s8 = sprintf('Prediction frequency (total REM) = \t%.3f Hz', predfreq_stim);
        
        figure
        circ_plot(phase_stim','hist',[], 20,false,true,'linewidth',2,'color','r');
        title({s,s0,s8,s1,s2,s3,s4});
        set(gcf, 'MenuBar', 'none')
        set(gcf, 'Toolbar', 'none')
        set(gcf, 'Position', get(0, 'Screensize'));
        
        %% Save roseplot STIMULATION condition
        figname = [DATAOUT, '/', EEG.part_num, '_accuracy_roseplot_STIM_' EEG.chanlocs.labels, '_', filt, '.jpg'];
        saveas (gcf,figname);
        
        %% SHAM CONDITION  %%
        
        %% Find markers SHAM condition
        shamtrials = [];
        for i = 1:size(EEG.event,2)
            if strcmpi(EEG.event(i).type,'103')
                shamtrials = [shamtrials, EEG.event(i).latency];
            end
        end
        
        %% Find prediction phases SHAM condition
        hil_sham = hil_eeg(shamtrials);
        phase_sham = angle(hil_sham)+(pi/2);
        
        %% Statistics & Analysis SHAM condition
        cmean_sham = circ_mean(phase_sham', [], 1);
        if cmean_sham <0
            cmean_sham = (2*pi + cmean_sham);
        end
        
        cmedian_sham = circ_median(phase_sham');
        if cmedian_sham<0
            cmedian_sham = (2*pi + cmedian_sham);
        end
        
        [P_value_sham, Z_value_sham] = circ_rtest(phase_sham');
        
        totalremtime = length(hil_eeg)/512;
        predfreq_sham = length(shamtrials)/totalremtime;
        
        %% Plot figure SHAM condition
        s=sprintf('%s: Sham condition',EEG.part_num);
        s0=sprintf('Total number of markers = \t%.0f',length(phase_sham));
        s1=sprintf(['Mean phase = \t%.2f ' char(176) 'C'],circ_rad2ang(cmean_sham));
        s2=sprintf(['Median phase = \t%.2f ' char(176) 'C'],circ_rad2ang(cmedian_sham));
        s3=sprintf('SD = \t%.2f',circ_rad2ang(circ_std(phase_sham')));
        s4=sprintf('SE = \t%.2f',circ_rad2ang(circ_std(phase_sham')/(sqrt(length(phase_sham')))));
        %s5=sprintf('Skewness (0=symm) = \t%.2f',circ_skewness(phase_sham'));
        %s6=sprintf('Kurtosis (0=norm.d.) = \t%.2f',circ_kurtosis(phase_sham'));
        %s7 = sprintf('Z-value = \t%.5f,\t P-value = \t%.5f', Z_value_sham, P_value_sham);
        s8 = sprintf('Prediction frequency (total REM) = \t%.3f Hz', predfreq_sham);
        
        figure
        circ_plot(phase_sham','hist',[], 20,false,true,'linewidth',2,'color','r');
        title({s,s0,s8,s1,s2,s3,s4});
        set(gcf, 'MenuBar', 'none')
        set(gcf, 'Toolbar', 'none')
        set(gcf, 'Position', get(0, 'Screensize'));
        
        %% Save roseplot SHAM condition
        figname = [DATAOUT, '/', EEG.part_num, '_accuracy_roseplot_SHAM_' EEG.chanlocs.labels, '_', filt, '.jpg'];
        saveas (gcf,figname);
        
        %% CREATE STRUCTURE WITH ACCURACY INFORMATION ALL PARTICIPANTS
        
        accuracy_info_allpp(k).Participant = pp;
        accuracy_info_allpp(k).Filter = filt;
        accuracy_info_allpp(k).TotalREMTime = totalremtime;
        accuracy_info_allpp(k).TotalMarkers = length(phase_stim)+length(phase_sham);
        
        accuracy_info_allpp(k).TotalStimMarkers = length(phase_stim);
        accuracy_info_allpp(k).PredictionFrequencyStimMarkers = predfreq_stim;
        accuracy_info_allpp(k).LatencyStimMarkers = stimtrials;
        accuracy_info_allpp(k).PhasesStimMarkers = double(phase_stim);
        accuracy_info_allpp(k).MeanPhaseStimMarkers = circ_rad2ang(cmean_stim);
        accuracy_info_allpp(k).MedianPhaseStimMarkers = circ_rad2ang(cmedian_stim);
        accuracy_info_allpp(k).SD = circ_rad2ang(circ_std(phase_stim'));
        accuracy_info_allpp(k).SE = circ_rad2ang(circ_std(phase_stim')/(sqrt(length(phase_stim'))))
        accuracy_info_allpp(k).PvalueStim = P_value_stim;
        accuracy_info_allpp(k).ZvalueStim = Z_value_stim;
        
        accuracy_info_allpp(k).TotalShamMarkers = length(phase_sham);
        accuracy_info_allpp(k).PredictionFrequencyShamMarkers = predfreq_sham;
        accuracy_info_allpp(k).LatencyShamMarkers = shamtrials;
        accuracy_info_allpp(k).PhasesShamMarkers = double(phase_sham);
        accuracy_info_allpp(k).MeanPhaseShamMarkers = circ_rad2ang(cmean_sham);
        accuracy_info_allpp(k).MedianPhaseShamMarkers = circ_rad2ang(cmedian_sham);
        accuracy_info_allpp(k).SD = circ_rad2ang(circ_std(phase_sham'));
        accuracy_info_allpp(k).SE = circ_rad2ang(circ_std(phase_sham')/(sqrt(length(phase_sham'))))
        accuracy_info_allpp(k).PvalueSham = P_value_sham;
        accuracy_info_allpp(k).ZvalueSham = Z_value_sham;
        
        k = k + 1;
        
    end
end

%% SAVE STRUCTURE WITH ACCURACY INFORMATION ALL PARTICIPANTS
save(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/Algorithm performance/Fz/Clean data roseplots/ACCURACY_INFORMATION_ALLPARTICIPANTS_CLEAN.mat'], 'accuracy_info_allpp');
