function [sm,Sm,X,V_estim] = getdata_hpss(dataset_path,magnitudes_path,ind,Nfft,hop,Nw,context_length)

%%% True sources

% Time-domain sources
sources_path = strcat(dataset_path,'Sources/Test/');
aux = dir(sources_path);
song_name = aux(ind+2).name;

bb = audioread(strcat(sources_path,song_name,'/bass.wav'));
dd = audioread(strcat(sources_path,song_name,'/drums.wav'));
oo = audioread(strcat(sources_path,song_name,'/other.wav'));
vv = audioread(strcat(sources_path,song_name,'/vocals.wav'));
harmo = mean(bb+oo+vv,2)'; percu = mean(dd,2)';
sm = [percu;harmo]+eps;

% STFT and remove initial context frames
Sm = STFT(sm,Nfft,hop,Nw);
Sm = Sm(:,context_length+1:end,:);


%%% Estimated voice magnitude
V_percu_estim = readNPY(strcat(magnitudes_path,'magn_spect_',int2str(ind-1),'.npy'))';

% Reshape the STFTs so orig and estim have the same size
T = min(size(V_percu_estim,2),size(Sm,2));
V_percu_estim = V_percu_estim(:,1:T);
Sm = Sm(:,1:T,:);
X = sum(Sm,3);
Vx = abs(X);
V_estim = zeros(size(Sm)); V_estim(:,:,1) = V_percu_estim; V_estim(:,:,2) = Vx-V_percu_estim;

% Time domain original sources
sm = iSTFT(Sm,Nfft,hop,Nw);

end