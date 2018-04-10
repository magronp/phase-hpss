%*************************************************************************
% Source separation with kernel additive modelling.
%*************************************************************************
% MATLAB implementation of the harmonic/percussive separation method 
% presented in "Kernel Spectrogram Models for Source Separation", presented
% at HSCMA'2014
% ------------------------------------------------------------------------
%
% This script performs separation of the harmonic and percussive parts of 
% music.
%
% Separation is performed through spatial kernel backfitting. Harmonic and
% percussive parts are modelled with different kernels.
%
% Supports many parameters, see the parameters section. In particular, you
% can activate a 'light" mode that performs less iterations and thus works
% faster.
%
% When using this for your work, please mention the following references:
% 
% @inproceedings{liutkusKAM14b, 
% AUTHOR = {A. Liutkus and Z. Rafii and B. Pardo and D. Fitzgerald and L. Daudet}, 
% TITLE = {{Kernel Spectrogram models for source separation}}, 
% BOOKTITLE = {IEEE. Workshop on Hands-free Speech Communication and Microphone Arrays (HSCMA 2014)},
% YEAR = {2014}, 
% MONTH = May, 
% ADDRESS = {Nancy, France}, 
% URL = {http://hal.inria.fr/hal-00959384} }
%
% and:
%
%
%@Conference{FitzgeraldKAMHPdafx,
%  Title                    = {Harmonic/Percussive Separation Using {K}ernel {A}dditive {M}odelling},
%  Author                   = {D. Fitzgerald  and A. Liutkus and Z. Rafii and B. Pardo and L. Daudet},
%  Booktitle                = {Proc. of the 7th Int. Conference on Digital Audio Effects (DAFx’14),},
%  Year                     = {2014},
% note = {submitted},
%}
%
%*************************************************************************
% Copyright (C) 2014, Inria, Antoine Liutkus
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU Affero General Public License as
%    published by the Free Software Foundation, either version 3 of the
%    License, or (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Affero General Public License for more details.
%
%    You should have received a copy of the GNU Affero General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function Se = kam_hpss_mono(X,Nfft,hop,Fs)


[F,T] = size(X);
I = 1;
J = 2;


%**************************************************************************
% Kernels
%**************************************************************************

%definition of the different frequency bands. 
fstop_bands = [250,3000,10000,22050];

%1) Stable harmonic source
%size of kernels for harmonic source. One rectangular kernel for each band
%for harmonic source, we should basically have flat kernels
harmonicStability = [[0,   0.1]; ... %bass: should be stable for really a long time
                     [0,  0.1]; ...  %medium stable for a short time
                     [5,  0.1]; ...   %mid : stable for moderate time but also in frequency a bit
                     [5,  0.4]];       %high: shoudl really be stable for a long time

%2) Percussive source
%size of kernels for percussive source. One rectangular kernel for each band
%for percussive source, we should have flatall kernels
percussiveStability = [[30,   0.02]; ... %bass drum should be a bit stable over time
                     [150,  0]; ...  
                     [200,  0.01]; ...
                     [500,0.1]]; 

niter= 6;          %number of iterations

%Creating kernels
%----------------
poscuts=floor(min(fstop_bands,Fs)/Fs*Nfft);
kernels = {};

%1) Harmonic part
kernels{end+1}={};
for band = 1:length(poscuts)
    kernels{end}{band}= buildRectangularKernel(harmonicStability(band,:),Nfft,Fs,hop);
end
%2) Percussive part
kernels{end+1} = {};
for band = 1:length(poscuts)
    kernels{end}{band}= buildRectangularKernel(percussiveStability(band,:),Nfft,Fs,hop);
end



% Initializing model
%-------------------
%PSD are initialized simply as mixtures PSD /J

S = single(repmat(sum(abs(X).^2,3),[1,1,J]))/J; 
%All spatial covariance matrices as initially identity
R = zeros(1,I,I,J); 
for j = 1:J
    R(1,:,:,j) = eye(I);
end
R = repmat(R,[F,1,1,1]);

%Performing estimation through kernel backfitting
%--------------------------------------------------------
for it=1:niter
    %Separating sources with current model
    out = posterior(X,S,'R',R,'Ymu',1);
    Ymu = out.Ymu;

    %Backfitting each source
    for j=1:J
        fprintf('Backfitting, iteration %d / %d, source %d / %d \n',it,niter,j,J);
        tempS = zeros(F,T);

        %learning spatial covariance and PSD from separated image
        [Z,Rj] = learnDSP(Ymu(:,:,:,j),1);

        %nan or inf will occur if all PSD estimates are 0, typical for voice below
        %minimal frequency
        Rj(isnan(Rj)) = 0;
        Rj(isinf(Rj)) = 0;
        R(:,:,:,j)=Rj;

        %median filter the estimated PSD Z
        if iscell(kernels{j})
            %we handle the case were we have several frequency bands
            pos = 1;
            for k = 1:length(poscuts)
                order= round(length(find(kernels{j}{k}))/2);
                tempS(pos:poscuts(k),:) = ordfilt2(Z(pos:poscuts(k),:),order,kernels{j}{k},'symmetric');
                pos = poscuts(k)+1;
            end
        else
            order= round(length(find(kernels{j}))/2);
            tempS = ordfilt2(Z,order,kernels{j},'symmetric');
        end
        S(:,:,j) = tempS;
    end

end

% Final separation
out = posterior(X,S,'R',R,'Ymu',1);
Ymu = squeeze(out.Ymu);

% Switch the order (percu first, harmo second)
Se = zeros(F,T,J);
Se(:,:,1) = Ymu(:,:,2);
Se(:,:,2) = Ymu(:,:,1);

end