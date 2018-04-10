function [P,R] = learnDSP(S,Niter)
%--------------------------------------------------------------------------
% Learn Power Spectral Density and spatial covariance matrix of a
% multichannel image
% 
% Input :
%   * S : FxTxJ : STFT or MDCT of the multichannel image
%
% Output : 
%   * P : FxT   : estimated DSP of signal
%   * R : FxJxJ : estimated spatial covariance matrix
%--------------------------------------------------------------------------
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



if nargin==1
    Niter = 15;
end
[F,T,J] = size(S);

P = max(eps,mean(abs(S).^2,3));
R = zeros(F,J,J);

for f = 1:F
    R(f,:,:) = eye(J);
end

EYE = R*1E-5;
Rss = zeros(F,T,J,J);


for j1 = 1:J
    for j2 = 1:J
        Rss(:,:,j1,j2) = S(:,:,j1).*conj(S(:,:,j2));
    end
end


Rinv = zeros(size(R));

if Niter >1 
    disp('Learning multichannel covariance and DSP...')
    str = [];
end
for it = 1:Niter
    if Niter>1
        fprintf(repmat('\b',1,length(str)));
        str=sprintf('       Iteration %2d of total %2d', it, Niter);
        fprintf('%s', str);
    end
    
    R(:,:,:) = squeeze(mean(bsxfun(@times,1./max(eps,P),Rss),2))+EYE;
    for f = 1:F
        R(isnan(R))=0;
        R(isinf(R))=0;
        Rinv(f,:,:) = pinv(squeeze(R(f,:,:)));
    end
    
    P = zeros(F,T);
    for j=1:J
        for j2 = 1:J
            P = P + real(bsxfun(@times, Rinv(:,j,j2), Rss(:,:,j2,j)));
        end
    end
    P = P/J;
end
if Niter>1
    fprintf('\n');
end
