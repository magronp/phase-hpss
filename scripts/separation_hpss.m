clear all; close all; clc;
stft_setting = 1;
set_settings_phasehpss;

% Separation
for ind=1:Nsongs
    
    clc; fprintf('song %d / %d \n',ind,Nsongs);
    
    % Get the data
    [sm,Sm,X,V_estim] = getdata_hpss(dataset_path,magnitudes_path,ind,Nfft,hop,Nw,context_length);
    [F,T,J] = size(Sm);
    
    % Phase retrieval
    S_estim = zeros(F,T,J,Nalgo);
    S_estim(:,:,:,1) = V_estim .* exp(1i*angle(X));
    S_estim(:,:,:,2) = pu_hpss(X,V_estim,hop,Npuiter);
    S_estim(:,:,:,3) = kam_hpss_mono(X,Nfft,hop,Fs);
    
    % Time-domain synthesis
    s_estim = zeros(J,length(sm),Nalgo);
    for al=1:Nalgo
        s_estim(:,:,al) = real(iSTFT(S_estim(:,:,:,al),Nfft,hop,Nw));
    end
    
    % Remove samples for which the estimation is irelevant
    s_estim = s_estim(:,context_length*hop+1:end,:);
    sm = sm(:,context_length*hop+1:end);
    
    % Record
    audiowrite(strcat(audio_path,int2str(ind),'_percu_orig.wav'),sm(1,:),Fs);
    audiowrite(strcat(audio_path,int2str(ind),'_harmo_orig.wav'),sm(2,:),Fs);
    for al = 1:Nalgo
        audiowrite(strcat(audio_path,int2str(ind),'_percu_estim_',int2str(al),'.wav'),squeeze(s_estim(1,:,al)),Fs);
        audiowrite(strcat(audio_path,int2str(ind),'_harmo_estim_',int2str(al),'.wav'),squeeze(s_estim(2,:,al)),Fs);
    end
    
end