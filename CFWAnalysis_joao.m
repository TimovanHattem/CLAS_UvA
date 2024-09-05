%% Loop to get CFW information
Fs = 512
for i = 1:size(ERP_info_allpp,2)
    clear cfs f temp_cfw
    Temp_cfw = double(ERP_info_allpp(i).AverageERPSham);
    [cfs,f] = cwt(Temp_cfw,Fs);
    CFW_info_allpp(i).Participant = ERP_info_allpp(i).Participant;
    CFW_info_allpp(i).cfs_sham = cfs;
    CFW_info_allpp(i).f_sham = f;
end

for i = 1:size(ERP_info_allpp,2)
    clear cfs f temp_cfw
    Temp_cfw = double(ERP_info_allpp(i).AverageERPStim);
    [cfs,f] = cwt(Temp_cfw,Fs);
    CFW_info_allpp(i).Participant = ERP_info_allpp(i).Participant;
    CFW_info_allpp(i).cfs_stim = cfs;
    CFW_info_allpp(i).f_stim = f;
end
%% Averaging STIM 
for ii = 1:size(CFW_info_allpp,2)
    
    
    tfamat_stim(:,:,ii) = CFW_info_allpp(ii).cfs_stim;
    tfamat_sham(:,:,ii) = CFW_info_allpp(ii).cfs_sham;

end
%%
P1 = mean(tfamat_stim,3);
P2 = mean(tfamat_sham,3);
%%
colorax = [-5 5];
close all
subplot 121

% sigLen = numel(sig);
sigLen = 768; %hardcoded because I know
t = (0:sigLen-1)/Fs;

surface(t,f,10*log10(abs(P1)));
% ylim([0 30])
axis square
xlabel('Time (s)')
ylabel('Frequency (Hz)')
zlabel('Magnitude')
title('Scalogram stim')
set(gca,'yscale','log')
shading interp
colorbar
caxis(colorax)
colormap
view([0 90])

subplot 122

% sigLen = numel(sig);
sigLen = 768; %hardcoded because I know
t = (0:sigLen-1)/Fs;

surface(t,f,10*log10(abs(P2)));
% ylim([0 30])
axis square
xlabel('Time (s)')
ylabel('Frequency (Hz)')
zlabel('Magnitude')
title('Scalogram sham')
set(gca,'yscale','log')
shading interp
colorbar
caxis(colorax)
colormap
view([0 90])