% Set the settings used in the experiments

Fs = 44100;
Nsongs = 50;

% STFT parameters
Nfft = 4096;
switch stft_setting
    case 1
        Nw = 2049;
        hop = 384;
        context_length = 10;
    case 2
        Nw = 4096;
        hop = 1024;   
        context_length = 10;
end


% Paths
dataset_path = 'dataset/';
magnitudes_path = strcat('magnitude_spectrograms/setting',int2str(stft_setting),'/');
audio_path = strcat('audio_files/setting',int2str(stft_setting),'/');
metrics_path = 'metrics/';

% Algorithms
algos = {'MTN-mixphase','MTN-PUiter','KAM'};
Nalgo = length(algos);
Npuiter = 50;
