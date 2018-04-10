clear all; close all; clc;
stft_setting = 1;
set_settings_phasehpss;

SDR = []; SIR = []; SAR = [];

for ind=1:Nsongs
    
    clc; fprintf('BSS - Song %d / %d \n',ind,Nsongs);
     
    % Original Files path
    s1 = audioread(strcat(audio_path,int2str(ind),'_percu_orig.wav'));
    s2 = audioread(strcat(audio_path,int2str(ind),'_harmo_orig.wav'));
    sm = [s1 s2]';
    
    % Loop over the algorithms
    sd_aux = []; si_aux = []; sa_aux = [];
   for al=1:Nalgo
       
       % Load estimated signals
       s_estim1 =  audioread(strcat(audio_path,int2str(ind),'_percu_estim_',int2str(al),'.wav'));
       s_estim2 =  audioread(strcat(audio_path,int2str(ind),'_harmo_estim_',int2str(al),'.wav'));
       se = [s_estim1 s_estim2]';
       
       % BSS - framewise (30s windows)
       [sdr,~,sir,sar] = bss_eval_images_framewise(se,sm);
       sd_aux = [sd_aux ; sdr];
       si_aux = [si_aux ; sir];
       sa_aux = [sa_aux ; sar];
       
   end
   
   SDR = [SDR sd_aux]; SIR = [SIR si_aux]; SAR = [SAR sa_aux];
   
end

% Remove NaN
sdaux = SDR; li = isnan(sdaux(1,:));
SDR(:,li) = []; SIR(:,li) = []; SAR(:,li) = [];

% Record
save(strcat(metrics_path,'bss_phase-hpss_setting',int2str(stft_setting),'.mat'),'SDR','SIR','SAR');