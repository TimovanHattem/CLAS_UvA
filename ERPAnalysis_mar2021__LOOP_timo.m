% ANALYSIS FOR ERP ANALYSIS THETA TARGETING
%
% Timo van Hattem & Joao Patriota
% Experiment Theta_Targeting
% Updated: 1-6-2021

%% SETTINGS

clear all
close all
clc

% Set pathS
if isempty(strfind(path, ['/home/Public/Projects/sleep-and-cognition/generic/matlab-scripts/toolboxes/eeglab_v2020', pathsep]))
    addpath('/home/Public/Projects/sleep-and-cognition/generic/matlab-scripts/toolboxes/eeglab_v2020');
end
addpath(genpath('/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/Raw data'));
addpath('/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/matlab-scripts/2. ZMax Analysis/Algorithm Accuracy'); % to get e.g. 'bandpass.m'
addpath(genpath('/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results'));
fprintf('Paths added!\n')

%% LOOP OVER ALL PARTICIPANTS
z = 1;

for i = 12:35
    
    % Clear workspace
    close all
    clc
    
    % Load EEG dataset in eeglab
    eeglab;
    eegfilename = ['P', num2str(i), '_epoched_3hz_clean_STIM.set'];
    eegpath = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/P', num2str(i)];
    EEG_stim = pop_loadset('filename',eegfilename,'filepath',eegpath);
    eegfilename = ['P', num2str(i), '_epoched_3hz_clean_SHAM.set'];
    eegpath = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/P', num2str(i)];
    EEG_sham = pop_loadset('filename',eegfilename,'filepath', eegpath);
    fprintf('EEG file loaded!\n')
    
    % Setting output file path
    pp = EEG_stim.part_num;
    if strcmpi(pp,'P66')
        pp = 'P16';
    end
    DATAOUT = ['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/',pp];
    
    % Averaging trials
    average_erp_stim = mean(EEG_stim.data,3);
    average_erp_sham = mean(EEG_sham.data,3);
    
    % Plot ERP image
    figure;
    x = [-256:1:511];
    plot(x,average_erp_stim, 'b', 'LineWidth', 1.5)
    hold on
    plot(x,average_erp_sham, 'r','LineWidth', 1.5)
    xlim([-256 511])
    set(gca,'XTick',[-253,0,256,511])
    set(gca,'xticklabel',({'-500', '0', '500', '1000'}))
    title(pp)
    xlabel('Time (ms)')
    ylabel(['Amplitude (', char(181),'V)'])
    line([0 0], [-15 15], 'color', 'black')
    legend('Stimulation', 'Sham')
    
    % Save figure
    figname = [DATAOUT, '/', EEG_stim.part_num, '_ERPs_3hz.jpg'];
    saveas (gcf,figname);
    close
    
    % Save ERP information in struct
    ERP_info_allpp_3hz(z).Participant = pp;
    ERP_info_allpp_3hz(z).AverageERPStim = average_erp_stim;
    ERP_info_allpp_3hz(z).AverageERPSham = average_erp_sham;
    
    z = z + 1;
    
    clearvars -except z ERP_info_allpp_3hz
end

save(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/TOTAL AVERAGE/ERP_INFORMATION_ALLPARTICIPANTS_3hz.mat'], 'ERP_info_allpp_3hz');

%% TOTAL AVERAGE

load('ERP_INFORMATION_ALLPARTICIPANTS_3hz.mat')

% Averaging trials
for i = 1:size(ERP_info_allpp_3hz,2)
    pp = ERP_info_allpp_3hz(i).Participant;
    total_erp = ERP_info_allpp_3hz(i).AverageERPStim;
    erpmat_stim(i,:) = total_erp;
    %save(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/',pp, '/', pp,'_totalerp_stim.mat'], 'total_erp');
    total_erp = ERP_info_allpp_3hz(i).AverageERPSham;
    erpmat_sham(i,:) = total_erp;
    %save(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/',pp, '/', pp,'_totalerp_sham.mat'], 'total_erp');
end
%save(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/TOTAL AVERAGE/totalerp_stim.mat'], 'erpmat_stim');
%save(['/home/Public/Projects/sleep-and-cognition/experiments/0016.Measuring-theta-during-REM-sleep/results/ThetaTargeting_2021/STA/TOTAL AVERAGE/totalerp_sham.mat'], 'erpmat_sham');
   
totalaverage_erp_stim = mean(erpmat_stim,1);
totalaverage_erp_sham = mean(erpmat_sham,1);

%% STATISTICS: CLUSTER-BASED PERMUTATION TEST

load('totalerp_stim.mat')
load('totalerp_sham.mat')

x = erpmat_stim';
y = erpmat_sham';
d = true;
p = 0.05;
num_permutations = 1000;
t = true;
num_clusters = [];

[clusters, p_values, t_sums, permutation_distribution ] = permutest(x, y, d, p, num_permutations, t, num_clusters);
   
%% PLOT FINAL ERP IMAGE

figure;
x = [-256:1:511];
plot(x,totalaverage_erp_stim, 'b', 'LineWidth', 1.15)
hold on
plot(x,totalaverage_erp_sham, 'r','LineWidth', 1.15)
stdshade(erpmat_stim,0.09,'b',x,[])
stdshade(erpmat_sham,0.09,'r',x,[])
xlim([-256 511])
ylim([-11.5 11.5])
set(gca,'XTick',[-205,-102,0,102,205,307,409,511], 'FontSize', 14)
set(gca,'xticklabel',({'-400', '-200', '0', '200', '400', '600', '800', '1000'}))
set(gca,'box','off')
%title('Grand average')
xlabel('Time (ms)', 'FontSize', 18, 'FontWeight', 'bold')
ylabel(['Amplitude (', char(181),'V)'], 'FontSize', 18, 'FontWeight', 'bold')
line([0 0], [-11.5 11.5], 'Color', 'black', 'LineStyle', '--')
a = line([(353-256) (398-256)], [-11.1 -11.1], 'Color', 'black', 'LineWidth', 15.2);
a.Color(4) = 0.3; 
% = p < 0.001
b = line([(414-256) (461-256)], [-11.1 -11.1], 'Color', 'black', 'LineWidth', 15.2);
b.Color(4) = 0.3;
% p < 0.001
c = line([(291-256) (341-256)], [-11.1 -11.1], 'Color', 'black', 'LineWidth', 15.2);
c.Color(4) = 0.3;
% p < 0.001
d = line([(480-256) (515-256)], [-11.1 -11.1], 'Color', 'black', 'LineWidth', 15.2);
d.Color(4) = 0.3;
% p = 0.031
legend('Stimulation', 'Sham', 'FontSize', 18, 'FontWeight', 'bold')
% h=gcf;
% set(h,'Position',[50 50 1200 800]);
% set(h,'PaperOrientation','landscape');
% print(gcf, '-dpdf', 'test.pdf')
