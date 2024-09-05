% ANALYSIS FOR TFA/PSD ANALYSIS THETA TARGETING
%
% Timo van Hattem & Joao Patriota
% Experiment Theta_Targeting
% Updated: 1-6-2021

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

%% LOOP OVER ALL PARTICIPANTS FOR INDIVIDUAL TFAs
% 
% z = 1;
% 
% for i = 12:35
%     
%     % Clear workspace
%     close all
%     clc
%     
%     % Load EEG dataset in eeglab
%     eeglab;
%     eegfilename = ['P', num2str(i), '_epoched_clean_STIM.set'];
%     eegpath = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/P', num2str(i)];
%     EEG_stim = pop_loadset('filename',eegfilename,'filepath',eegpath);
%     eegfilename = ['P', num2str(i), '_epoched_clean_SHAM.set'];
%     eegpath = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/P', num2str(i)];
%     EEG_sham = pop_loadset('filename',eegfilename,'filepath', eegpath);
%     fprintf('EEG file loaded!\n')
%     
%     % Setting output file path
%     pp = EEG_stim.part_num;
%     if strcmpi(pp,'P66')
%         pp = 'P16';
%     end
%     DATAOUT = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/',pp];
%     
%     % Reshape data
%     for j = 1:size(EEG_stim.data,3)
%         data_stim(j,:) = EEG_stim.data(:,:,j);
%         EEG_data_stim = double(data_stim);
%         %save(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/',pp, '/', pp,'_epochs_stim_nonset.mat'], 'EEG_data_stim');
%     end
%     for k = 1:size(EEG_sham.data,3)
%         data_sham(k,:) = EEG_sham.data(:,:,k);
%         EEG_data_sham = double(data_sham);
%         %save(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/',pp, '/', pp,'_epochs_sham_nonset.mat'], 'EEG_data_sham');
%     end
%     
%     ERP_average_stim_evoked = mean(EEG_data_stim,1);
%     ERP_average_sham_evoked = mean(EEG_data_sham,1);
%     
%     figure(1)
%     window_spectrogram = 128;
%     nooverlap_spectrogram = 120;
%     [s1,F1,T1,P1]=spectrogram(ERP_average_stim_evoked,window_spectrogram,nooverlap_spectrogram,1024,512); % Computes the spectrogram
%     plot3([0 0],[-10 40],[0 600],'w','LineWidth',3)
%     hold on
%     plot3([-0.375 0.875],[3 3],[0 200],'--w','LineWidth',1)
%     plot3([-0.375 0.875],[8 8],[0 200],'--w','LineWidth',1)
%     hold on
%     T1 = T1-0.5;
%     surf(T1,F1,10*log10(abs(P1)),'EdgeColor','none');
%     axis xy; axis tight; view(0,90);
%     xlabel('Time', 'FontSize', 18);
%     ylabel('Frequency (Hz)', 'FontSize', 18);
%     title('Spectrogram Stimulation')
%     colorbar
%     ylim([2 30])
%     caxis([-25 12])
%     
%     % Save figure
%     figname = [DATAOUT, '/', EEG_stim.part_num, '_TFA_STIM.jpg'];
%     saveas (gcf,figname);
%         
%     figure(2)
%     [s2,F2,T2,P2]=spectrogram(ERP_average_sham_evoked,window_spectrogram,nooverlap_spectrogram,1024,512); % Computes the spectrogra
%     T2 = T2-0.5;
%     surf(T2,F2,10*log10(abs(P2)),'EdgeColor','none');
%     hold on
%     plot3([0 0],[-10 40],[0 600],'w','LineWidth',3)
%     plot3([-0.375 0.875],[3 3],[0 200],'--w','LineWidth',1)
%     plot3([-0.375 0.875],[8 8],[0 200],'--w','LineWidth',1)
%     axis xy; axis tight; view(0,90);
%     xlabel('Time', 'FontSize', 18);
%     ylabel('Frequency (Hz)', 'FontSize', 18);
%     title('Spectogram Sham')
%     colorbar
%     ylim([2 30])
%     caxis([-25 12])
%     
%     % Save figure
%     figname = [DATAOUT, '/', EEG_stim.part_num, '_TFA_SHAM.jpg'];
%     saveas (gcf,figname);
%         
%     figure(3)
%     surf(T2,F2,10*log10(abs(P1-P2)),'EdgeColor','none');
%     hold on
%     plot3([0 0],[-10 40],[0 600],'w','LineWidth',3)
%     plot3([-0.375 0.875],[3 3],[0 200],'--w','LineWidth',1)
%     plot3([-0.375 0.875],[8 8],[0 200],'--w','LineWidth',1)
%     axis xy; axis tight; view(0,90);
%     xlabel('Time', 'FontSize', 18);
%     ylabel('Frequency (Hz)', 'FontSize', 18);
%     title('Contrast (Stim-Sham)')
%     colorbar
%     ylim([2 30])
%     caxis([-5 10])
%     
%     % Save figure
%     figname = [DATAOUT, '/', EEG_stim.part_num, '_TFA_CONTRAST.jpg'];
%     saveas (gcf,figname);
%         
%     % Save TFA information in struct
%     TFA_info_allpp(z).Participant = pp;
%     TFA_info_allpp(z).TFAStimT = T1;
%     TFA_info_allpp(z).TFAStimF = F1;
%     TFA_info_allpp(z).TFAStimP = P1;
%     TFA_info_allpp(z).TFAShamT = T2;
%     TFA_info_allpp(z).TFAShamF = F2;
%     TFA_info_allpp(z).TFAShamP = P2;
%     
%     z = z + 1;
%     
%     clearvars -except z TFA_info_allpp
%     
% end
% 
% save(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/TOTAL AVERAGE/TFA_INFORMATION_ALLPARTICIPANTS.mat'], 'TFA_info_allpp');

%% STATISTICS: CLUSTER-BASED PERMUTATION TEST TFA (non-baseline corrected)
close all 

load('TFA_INFORMATION_ALLPARTICIPANTS.mat');

for i = 1:size(TFA_info_allpp,2)
    tfamatstat_stim(:,:,i) = TFA_info_allpp(i).TFAStimP;
    tfamatstat_sham(:,:,i) = TFA_info_allpp(i).TFAShamP;
end

x = tfamatstat_stim(1:61,:,:);
y = tfamatstat_sham(1:61,:,:);
d = true;
p = 0.05;
num_permutations = 1000;
t = true;
num_clusters = [];

[clusters, p_values, t_sums, permutation_distribution ] = permutest(x, y, d, p, num_permutations, t, num_clusters);

signfcluster = clusters{1};
zeromat = zeros(61,81);
zeromat(signfcluster) = 200;
%figure(1)
%pcolor(zeromat);
%figure(2)
%surf(zeromat);
zerosadd = zeros(452,81);
zeromatend = vertcat(zeromat, zerosadd);
%contour(zeromatend)

%% TOTAL AVERAGE TFA NON-BASELINE CORRECTED
close all

load('TFA_INFORMATION_ALLPARTICIPANTS.mat');
load('Contour_TFA.mat')

colorax1 = [-15 15];
colorax2 = [-4 15];

for i = 1:size(TFA_info_allpp,2)
    tfamat_stim(:,:,i) = TFA_info_allpp(i).TFAStimP;
    tfamat_sham(:,:,i) = TFA_info_allpp(i).TFAShamP;
end

P1 = mean(tfamat_stim,3);
T1 = TFA_info_allpp(1).TFAStimT;
F1 = TFA_info_allpp(1).TFAStimF;

P2 = mean(tfamat_sham,3);
T2 = TFA_info_allpp(1).TFAShamT;
F2 = TFA_info_allpp(1).TFAShamF;

figure(1)
colormap jet
plot3([0 0],[-10 40],[0 600],'w','LineWidth',3)
hold on
plot3([-0.375 0.875],[4 4],[0 200],'--w','LineWidth',1)
plot3([-0.375 0.875],[8 8],[0 200],'--w','LineWidth',1)
hold on
surf(T1,F1,10*log10(abs(P1)),'EdgeColor','none');
shading interp
axis xy; axis tight; view(0,90);
xlabel('Time (ms)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
%title('Spectrogram Stimulation', 'FontSize', 12, 'FontWeight', 'bold')
c = colorbar;
ylim([2 30])
caxis(colorax1)
c.Label.String = 'Power (dB)';
c.Label.FontSize = 18;
c.Label.FontWeight = 'bold';
axis square 
set(gca,'FontSize',14)
set(gca,'xticklabel',({'-200', '0', '200', '400', '600', '800'}))

figure(2)
colormap jet
surf(T2,F2,10*log10(abs(P2)),'EdgeColor','none');
shading interp
hold on
plot3([0 0],[-10 40],[0 600],'w','LineWidth',3)
plot3([-0.375 0.875],[4 4],[0 200],'--w','LineWidth',1)
plot3([-0.375 0.875],[8 8],[0 200],'--w','LineWidth',1)
axis xy; axis tight; view(0,90);
xlabel('Time (ms)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
%title('Spectogram Sham', 'FontSize', 12, 'FontWeight', 'bold')
c = colorbar;
ylim([2 30])
caxis(colorax1)
c.Label.String = 'Power (dB)';
c.Label.FontSize = 18;
c.Label.FontWeight = 'bold';
axis square
set(gca,'FontSize',14)
set(gca,'xticklabel',({'-200', '0', '200', '400', '600', '800'}))
set(gca,'Box', 'on', 'Color', 'white')

P3 = (P1./P2);
figure(3)
colormap jet
surf(T2,F2,10*log10(abs(P3)),'EdgeColor','none');
shading interp
hold on
plot3([0 0],[-10 40],[0 600],'w','LineWidth',3)
plot3([-0.375 0.875],[4 4],[0 200],'--w','LineWidth',1)
plot3([-0.375 0.875],[8 8],[0 200],'--w','LineWidth',1)
axis xy; axis tight; view(0,90);
xlabel('Time (ms)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
%title('Contrast (Stimulation/Sham)', 'FontSize', 12, 'FontWeight', 'bold')
c = colorbar;
ylim([2 30])
caxis(colorax2)
c.Label.String = 'Relative power (Stimulation/Sham)';
c.Label.FontSize = 18;
c.Label.FontWeight = 'bold';
axis square
set(gca,'FontSize',14)
set(gca,'xticklabel',({'-200', '0', '200', '400', '600', '800'}))
hold on 
contour3(T1, F1, zeromatend, 'Color', 'k', 'LineWidth', 0.1, 'LineStyle', '-');

%% STATISTICS: CLUSTER-BASED PERMUTATION TEST TFA (baseline corrected)
close all 

load('TFA_INFORMATION_ALLPARTICIPANTS_BASELINECORRECTED.mat');

for i = 1:size(TFA_info_allpp_base,2)
    tfamatstat_stim(:,:,i) = TFA_info_allpp_base(i).TFAStimP;
    tfamatstat_sham(:,:,i) = TFA_info_allpp_base(i).TFAShamP;
end

x = tfamatstat_stim(1:61,:,:);
y = tfamatstat_sham(1:61,:,:);
d = true;
p = 0.05;
num_permutations = 1000;
t = true;
num_clusters = [];

[clusters, p_values, t_sums, permutation_distribution ] = permutest(x, y, d, p, num_permutations, t, num_clusters);

signfcluster = clusters{1};
zeromat = zeros(61,81);
zeromat(signfcluster) = 200;
%figure(1)
%pcolor(zeromat);
%figure(2)
%surf(zeromat);
zerosadd = zeros(452,81);
zeromatendbase = vertcat(zeromat, zerosadd);
%contour(zeromatend)

%% TOTAL AVERAGE TFA BASELINE CORRECTED

% load('TFA_INFORMATION_ALLPARTICIPANTS.mat');
% TFA_info_allpp_base = TFA_info_allpp;
% 
% Evoked = [33:81];
% Baseline = [1:32];
% 
% for i = 1:size(TFA_info_allpp_base,2)
%     baselinepowerstim = mean(TFA_info_allpp_base(i).TFAStimP(:,Baseline),2);
%     TFA_info_allpp_base(i).TFAStimP = TFA_info_allpp_base(i).TFAStimP./baselinepowerstim;
%     baselinepowersham = mean(TFA_info_allpp_base(i).TFAShamP(:,Baseline),2);
%     TFA_info_allpp_base(i).TFAShamP = TFA_info_allpp_base(i).TFAShamP./baselinepowersham;
% end

load('TFA_INFORMATION_ALLPARTICIPANTS_BASELINECORRECTED.mat');
load('Contour_TFA_base.mat');

colorax1 = [-4 6];
colorax2 = [1 15];

for i = 1:size(TFA_info_allpp_base,2)
    tfamat_stim(:,:,i) = TFA_info_allpp_base(i).TFAStimP;
    tfamat_sham(:,:,i) = TFA_info_allpp_base(i).TFAShamP;
end

P1 = mean(tfamat_stim,3);
T1 = TFA_info_allpp_base(1).TFAStimT;
F1 = TFA_info_allpp_base(1).TFAStimF;

P2 = mean(tfamat_sham,3);
T2 = TFA_info_allpp_base(1).TFAShamT;
F2 = TFA_info_allpp_base(1).TFAShamF;

figure(1)
colormap jet
plot3([0 0],[-10 40],[0 600],'w','LineWidth',3)
hold on
plot3([-0.375 0.875],[4 4],[0 200],'--w','LineWidth',1)
plot3([-0.375 0.875],[8 8],[0 200],'--w','LineWidth',1)
hold on
surf(T1,F1,10*log10(abs(P1)),'EdgeColor','none');
shading interp
axis xy; axis tight; view(0,90);
xlabel('Time (ms)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
%title('Spectrogram Stimulation', 'FontSize', 12, 'FontWeight', 'bold')
c = colorbar;
ylim([2 30])
caxis(colorax1)
c.Label.String = 'Power (dB)';
c.Label.FontSize = 18;
c.Label.FontWeight = 'bold';
axis square 
set(gca,'FontSize',14)
set(gca,'xticklabel',({'-200', '0', '200', '400', '600', '800'}))

figure(2)
colormap jet
surf(T2,F2,10*log10(abs(P2)),'EdgeColor','none');
shading interp
hold on
plot3([0 0],[-10 40],[0 600],'w','LineWidth',3)
plot3([-0.375 0.875],[4 4],[0 200],'--w','LineWidth',1)
plot3([-0.375 0.875],[8 8],[0 200],'--w','LineWidth',1)
axis xy; axis tight; view(0,90);
xlabel('Time (ms)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
%title('Spectogram Sham', 'FontSize', 12, 'FontWeight', 'bold')
c = colorbar;
ylim([2 30])
caxis(colorax1)
c.Label.String = 'Power (dB)';
c.Label.FontSize = 18;
c.Label.FontWeight = 'bold';
axis square
set(gca,'FontSize',14)
set(gca,'xticklabel',({'-200', '0', '200', '400', '600', '800'}))
set(gca,'Box', 'on', 'Color', 'white')

P3 = (P1./P2);
figure(3)
colormap jet
surf(T2,F2,10*log10(abs(P3)),'EdgeColor','none');
shading interp
hold on
plot3([0 0],[-10 40],[0 600],'w','LineWidth',3)
plot3([-0.375 0.875],[4 4],[0 200],'--w','LineWidth',1)
plot3([-0.375 0.875],[8 8],[0 200],'--w','LineWidth',1)
axis xy; axis tight; view(0,90);
xlabel('Time (ms)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Frequency (Hz)', 'FontSize', 18, 'FontWeight', 'bold');
%title('Contrast (Stimulation/Sham)', 'FontSize', 12, 'FontWeight', 'bold')
c = colorbar;
ylim([2 30])
caxis(colorax2)
c.Label.String = 'Relative power (Stimulation/Sham)';
c.Label.FontSize = 18;
c.Label.FontWeight = 'bold';
axis square
set(gca,'FontSize',14)
set(gca,'xticklabel',({'-200', '0', '200', '400', '600', '800'}))
hold on 
contour3(T1, F1, zeromatendbase, 'Color', 'k', 'LineWidth', 0.1, 'LineStyle', '-');

%% LOOP OVER ALL PARTICIPANTS FOR INDIVIDUAL PSDs
% close all
% 
% z = 1;
% for i = 12:35
%     pp = ['P', num2str(i)];
%     load(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/',pp, '/', pp,'_epochs_stim_nonset.mat']);
%     load(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/',pp, '/', pp,'_epochs_sham_nonset.mat']);
%     EEG_data_stim = EEG_data_stim(:,257:768);
%     EEG_data_sham = EEG_data_sham(:,257:768);
%     ERP_average_stim_evoked = mean(EEG_data_stim,1);
%     ERP_average_sham_evoked = mean(EEG_data_sham,1);
%     
%     % Pwelch for every line
%     theta_range = [4 8];
%     beta_range = [12 30];
%     total_range = [0.3 30];
%     window = 125;
%     nooverlap = [];
%     fftsize = 2^10;
%     Fs = 512;
%     
%     [Pwel_stim, Fwel_stim] = pwelch(ERP_average_stim_evoked,window,nooverlap,fftsize,Fs);
%     bandpowertheta_stim = bandpower(Pwel_stim,Fwel_stim,theta_range,'psd'); %band on theta
%     bandpowerbeta_stim = bandpower(Pwel_stim, Fwel_stim,beta_range,'psd');
%     bandpowertotal_stim = bandpower(Pwel_stim, Fwel_stim, total_range,'psd');
%     
%     [Pwel_sham, Fwel_sham] = pwelch(ERP_average_sham_evoked,window,nooverlap,fftsize,Fs);
%     bandpowertheta_sham = bandpower(Pwel_sham,Fwel_sham,theta_range,'psd'); %band on theta
%     bandpowerbeta_sham = bandpower(Pwel_sham,Fwel_sham,beta_range,'psd');
%     bandpowertotal_sham = bandpower(Pwel_sham, Fwel_stim, total_range,'psd');
%         
%     ratiopower = (Pwel_stim./Pwel_sham);
%     ratiobandpowertheta = (bandpowertheta_stim/bandpowertheta_sham);
%     ratiobandpowerbeta = (bandpowerbeta_stim/bandpowerbeta_sham);
%     
%     relativepowertheta_stim = (bandpowertheta_stim/bandpowertotal_stim)*100;
%     relativepowertheta_sham = (bandpowertheta_sham/bandpowertotal_sham)*100;
%     relativepowerbeta_stim = (bandpowertheta_stim/bandpowertotal_stim)*100;
%     relativepowerbeta_sham = (bandpowerbeta_sham/bandpowertotal_sham)*100;
%     
%     PSD_info_allpp(z).Participant = pp;
%     PSD_info_allpp(z).PwelStim = Pwel_stim;
%     PSD_info_allpp(z).PwelSham = Pwel_sham;
%     PSD_info_allpp(z).FwelStim = Fwel_stim;
%     PSD_info_allpp(z).FwelSham = Fwel_sham;
%     PSD_info_allpp(z).PwelRatio = ratiopower;
%     PSD_info_allpp(z).BandpowerThetaStim = bandpowertheta_stim;
%     PSD_info_allpp(z).BandpowerThetaSham = bandpowertheta_sham;
%     PSD_info_allpp(z).BandpowerBetaStim = bandpowerbeta_stim;
%     PSD_info_allpp(z).BandpowerBetaSham = bandpowerbeta_sham;
%     PSD_info_allpp(z).BandpowerTotalStim = bandpowertotal_stim;
%     PSD_info_allpp(z).BandpowerTotalSham = bandpowertotal_sham;
%     PSD_info_allpp(z).BandpowerThetaRatio = ratiobandpowertheta;
%     PSD_info_allpp(z).BandpowerBetaRatio = ratiobandpowerbeta;
%     PSD_info_allpp(z).RelativePowerThetaStim = relativepowertheta_stim;
%     PSD_info_allpp(z).RelativePowerThetaSham = relativepowertheta_sham;
%     PSD_info_allpp(z).RelativePowerBetaStim = relativepowerbeta_stim;
%     PSD_info_allpp(z).RelativePowerBetaSham = relativepowerbeta_sham;
%     
%     z = z + 1;
% end
%
% save(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/TOTAL AVERAGE/PSD_INFORMATION_ALLPARTICIPANTS.mat'], 'PSD_info_allpp');

%% STATISTICS: CLUSTER-BASED PERMUTATION TEST PSD

load('totalpsd_stim.mat')
load('totalpsd_sham.mat')

x = psdmat_stim(1:61,:);
y = psdmat_sham(1:61,:);
d = true;
p = 0.05;
num_permutations = 1000;
t = true;
num_clusters = [];

[clusters, p_values, t_sums, permutation_distribution ] = permutest(x, y, d, p, num_permutations, t, num_clusters);
 
%% TOTAL AVERAGE PSD (non-baseline corrected)

close all

load('PSD_INFORMATION_ALLPARTICIPANTS.mat');
theta_range = [4 8];

for i = 1:size(PSD_info_allpp,2)
    psdmat_stim(:,i) = PSD_info_allpp(i).PwelStim;
    psdmat_sham(:,i) = PSD_info_allpp(i).PwelSham;
    psdmat_ratio(:,i) = PSD_info_allpp(i).PwelRatio;
end
%save(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/TOTAL AVERAGE/totalpsd_stim.mat'], 'psdmat_stim');
%save(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/TOTAL AVERAGE/totalpsd_sham.mat'], 'psdmat_sham');
 
meanpower_stim = mean(psdmat_stim, 2);
meanpower_sham = mean(psdmat_sham,2);
meanpower_ratio = mean(psdmat_ratio,2);

% Plot Power spectrum
figure(1)
plot(PSD_info_allpp(1).FwelStim,10*log10(meanpower_stim), 'LineWidth', 1.15, 'Color', 'blue')                % Plot for the baseline window
hold on
plot(PSD_info_allpp(1).FwelSham,10*log10(meanpower_sham), 'LineWidth', 1.15, 'Color', 'red')            % Plot for the presentation window
line([theta_range(1) theta_range(1)], [-40 20], 'Color', 'black', 'LineStyle', '--'); %Create a tick mark at x = t1(i) with a height of 100
line([theta_range(2) theta_range(2)], [-40 20], 'Color', 'black', 'LineStyle', '--'); %Create a tick mark at x = t1(i) with a height of 100
xlim([2 30])
ylim([-40 10])
%title('PSD', 'FontSize', 12,'FontWeight', 'bold' )
xlabel('Frequency (Hz)','FontSize', 18,'FontWeight', 'bold')
ylabel('Power (dB)', 'FontSize', 18, 'FontWeight', 'bold')
axis square 
a = line([PSD_info_allpp(1).FwelStim(1) PSD_info_allpp(1).FwelStim(38)], [-39.15 -39.15], 'Color', 'black', 'LineWidth', 15.2);
a.Color(4) = 0.3;
legend('Stimulation','Sham', 'FontSize', 18, 'FontWeight', 'bold')
set(gca,'box','off')
set(gca,'FontSize',14)

figure(2)
plot(PSD_info_allpp(1).FwelStim,meanpower_ratio, 'LineWidth', 1.15)
%line([1 40], [1 1], 'Color', 'black','LineStyle','--'); %Create a tick mark at x = t1(i) with a height of 100
yline(1, '--', 'y = 1', 'FontSize', 12)
xlim([2 30])
ylim([-1 11])
%title('Power ratio Stimulation/Sham', 'FontSize', 12,'FontWeight', 'bold' )
xlabel('Frequency (Hz)','FontSize', 18,'FontWeight', 'bold')
ylabel('Relative power (Stimulation/Sham)','FontSize', 18, 'FontWeight', 'bold')
axis square
set(gca,'box','off')
set(gca,'FontSize',14)

%% THETA/BETA bandpower
close all

load('PSD_INFORMATION_ALLPARTICIPANTS.mat');

matrixtheta_stim = [PSD_info_allpp.BandpowerThetaStim];
matrixtheta_sham = [PSD_info_allpp.BandpowerThetaSham];
matrixbeta_stim = [PSD_info_allpp.BandpowerBetaStim];
matrixbeta_sham = [PSD_info_allpp.BandpowerBetaSham];

[a,b] = vartest2(matrixtheta_stim,matrixtheta_sham);
[aa,bb] = vartest2(matrixbeta_stim,matrixbeta_sham);
[c,d] = jbtest(matrixtheta_stim);
[f,g] = jbtest(matrixtheta_sham);
[h,i] = jbtest(matrixbeta_stim);
[j,k] = jbtest(matrixbeta_sham);

median_thetastim = median(matrixtheta_stim);
median_thetasham = median(matrixtheta_sham);
median_betastim = median(matrixbeta_stim);
median_betasham = median(matrixbeta_sham);
iqr_thetastim = iqr(matrixtheta_stim);
iqr_thetasham = iqr(matrixtheta_sham);
iqr_betastim = iqr(matrixbeta_stim);
iqr_betasham = iqr(matrixbeta_sham);

[pt,ht,statst] = signrank(matrixtheta_stim,matrixtheta_sham);
[pb,hb,statsb] = signrank(matrixbeta_stim,matrixbeta_sham);

matrixtheta_stim([11,14]) = [7.5]; %remove outliers
matrixtheta_sham(23)= [0.5]; %remove outliers
matrixbeta_sham(23)= [0.15]; %remove outliers
matrixtheta = [matrixtheta_stim' matrixtheta_sham'];
matrixbeta = [matrixbeta_stim' matrixbeta_sham'];

% figure(1)
% meanbandthetastim = mean([PSD_info_allpp.BandpowerThetaStim],2);
% stdbandthetastim = std([PSD_info_allpp.BandpowerThetaStim]);
% meanbandthetasham = mean([PSD_info_allpp.BandpowerThetaSham],2);
% stdbandthetasham = std([PSD_info_allpp.BandpowerThetaSham]);
% bar([meanbandthetastim, meanbandthetasham])
% hold on
% errorbar([1:2],[meanbandthetastim, meanbandthetasham], [stdbandthetastim, stdbandthetasham],[stdbandthetastim, stdbandthetasham], 'LineStyle', 'none', 'Color', 'black');
% set(gca,'XTick',[1,2])
% set(gca,'xticklabel',({'Stimulation', 'Sham'}))
% ylim([0,10])
% title('Theta bandpower comparison')
% xlabel('Condition', 'FontSize', 12)
% ylabel('Theta bandpower post-stimulus (a.u.)', 'FontSize', 12)

% figure(1)
% boxplot(matrixtheta, 'Colors', 'br')
% set(gca,'XTick',[1,2])
% set(gca,'xticklabel',({'Stimulation', 'Sham'}))
% %title('Theta', 'FontSize', 12, 'FontWeight', 'bold')
% xlabel('Condition', 'FontSize', 12, 'FontWeight', 'bold')
% ylabel('Theta bandpower', 'FontSize', 12, 'FontWeight', 'bold')
% ylim([0,9])

figure(1);
x = 1:2;
colors = [0 0 1; 1 0 0];
ax = axes();
hold(ax);
for i=1:2
    boxchart(x(i)*ones(size(matrixtheta(:,i))), matrixtheta(:,i), 'BoxFaceColor', colors(i,:))
end
set(gca,'XTick',[1,2], 'FontSize', 14)
set(gca,'xticklabel',({'Stimulation', 'Sham'}))
%title('Theta', 'FontSize', 12, 'FontWeight', 'bold')
xlabel('Condition', 'FontSize', 18, 'FontWeight', 'bold')
ylabel('Theta bandpower (a.u.)', 'FontSize', 18, 'FontWeight', 'bold')
ylim([0,8.5])

% figure(2)
% meanbandbetastim = mean([PSD_info_allpp.BandpowerBetaStim],2);
% stdbandbetastim = std([PSD_info_allpp.BandpowerBetaStim]);
% meanbandbetasham = mean([PSD_info_allpp.BandpowerBetaSham],2);
% stdbandbetasham = std([PSD_info_allpp.BandpowerBetaSham]);
% bar([meanbandbetastim, meanbandbetasham])
% hold on
% errorbar([1:2],[meanbandbetastim, meanbandbetasham], [stdbandbetastim, stdbandbetasham],[stdbandbetastim, stdbandbetasham], 'LineStyle', 'none', 'Color', 'black');
% set(gca,'XTick',[1,2])
% set(gca,'xticklabel',({'Stimulation', 'Sham'}))
% ylim([0,0.2])
% title('Beta bandpower comparison')
% xlabel('Condition', 'FontSize', 12)
% ylabel('Beta bandpower post-stimulus (a.u.)', 'FontSize', 12)

figure(2)
x = 1:2;
colors = [0 0 1; 1 0 0];
ax = axes();
hold(ax);
for i=1:2
    boxchart(x(i)*ones(size(matrixbeta(:,i))), matrixbeta(:,i), 'BoxFaceColor', colors(i,:))
end
set(gca,'XTick',[1,2], 'FontSize', 14)
set(gca,'xticklabel',({'Stimulation', 'Sham'}))
%title('Theta', 'FontSize', 12, 'FontWeight', 'bold')
xlabel('Condition', 'FontSize', 18, 'FontWeight', 'bold')
ylabel('Beta bandpower (a.u.)', 'FontSize', 18, 'FontWeight', 'bold')
ylim([0,0.25])

%% Relative power
% 
% figure(1)
% meanbandthetastim = mean([PSD_info_allpp.RelativePowerThetaStim],2);
% stdbandthetastim = std([PSD_info_allpp.RelativePowerThetaStim]);
% meanbandthetasham = mean([PSD_info_allpp.RelativePowerThetaSham],2);
% stdbandthetasham = std([PSD_info_allpp.RelativePowerThetaSham]);
% bar([meanbandthetastim, meanbandthetasham])
% hold on
% errorbar([1:2],[meanbandthetastim, meanbandthetasham], [stdbandthetastim, stdbandthetasham],[stdbandthetastim, stdbandthetasham], 'LineStyle', 'none', 'Color', 'black');
% set(gca,'XTick',[1,2])
% set(gca,'xticklabel',({'Stimulation', 'Sham'}))
% ylim([0,50])
% title('Relative Theta bandpower comparison')
% xlabel('Condition', 'FontSize', 12)
% ylabel('Relative Theta bandpower post-stimulus (a.u.)', 'FontSize', 12)
%
% figure(1)
% boxplotmatrixtheta(:,1) = [PSD_info_allpp.BandpowerThetaStim];
% boxplotmatrixtheta(:,2)= [PSD_info_allpp.BandpowerThetaSham];
% boxplot(boxplotmatrixtheta)
% set(gca,'XTick',[1,2])
% set(gca,'xticklabel',({'Stimulation', 'Sham'}))
% title('Relative Theta bandpower comparison')
% xlabel('Condition', 'FontSize', 12)
% ylabel('Relative Theta bandpower post-stimulus (a.u.)', 'FontSize', 12)
% 
% figure(2)
% meanbandbetastim = mean([PSD_info_allpp.RelativePowerBetaStim],2);
% stdbandbetastim = std([PSD_info_allpp.RelativePowerBetaStim]);
% meanbandbetasham = mean([PSD_info_allpp.RelativePowerBetaSham],2);
% stdbandbetasham = std([PSD_info_allpp.RelativePowerBetaSham]);
% bar([meanbandbetastim, meanbandbetasham])
% hold on
% errorbar([1:2],[meanbandbetastim, meanbandbetasham], [stdbandbetastim, stdbandbetasham],[stdbandbetastim, stdbandbetasham], 'LineStyle', 'none', 'Color', 'black');
% set(gca,'XTick',[1,2])
% set(gca,'xticklabel',({'Stimulation', 'Sham'}))
% ylim([0,10])
% title('Relative Beta bandpower comparison')
% xlabel('Condition', 'FontSize', 12)
% ylabel('Relative Beta bandpower post-stimulus (a.u.)', 'FontSize', 12)
%
% figure(2)
% boxplotmatrixbeta(:,1) = [PSD_info_allpp.BandpowerBetaStim];
% boxplotmatrixbeta(:,2)= [PSD_info_allpp.BandpowerBetaSham];
% boxplot(boxplotmatrixbeta)
% set(gca,'XTick',[1,2])
% set(gca,'xticklabel',({'Stimulation', 'Sham'}))
% title('Relative Beta bandpower comparison')
% xlabel('Condition', 'FontSize', 12)
% ylabel('Relative Beta bandpower post-stimulus (a.u.)', 'FontSize', 12)

%% THETA BANDPOWER RATIO PER PARTICIPANT
close all

load('PSD_INFORMATION_ALLPARTICIPANTS.mat');

figure(1)
%title('Distribution theta ratio participants post-stimulus')
plot([1:24],[PSD_info_allpp.BandpowerThetaRatio],'.','MarkerSize',37, 'Color', 'b')
xlabel('Subject #', 'FontSize', 18,'FontWeight', 'bold')
ylabel('Relative theta power increase','FontSize', 18,'FontWeight', 'bold')
set(gca,'XTick',[1:24], 'FontSize', 14)
set(gca,'xticklabel',({'1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19','20', '21', '22', '23', '24'}))
ylim([-1,40])
%line([1 24], [1 1], 'color', 'black','LineStyle','--');%Create a tick mark at x = t1(i) with a height of 100
%line([1 24], [9.0783 9.0783], 'color', 'black','LineStyle','--');
yline(1, '--', 'y = 1', 'FontSize', 12)
set(gca,'box','off')