function out = posterior(X,P,varargin)
%--------------------------------------------------------------------------
% Posterior distribution of sources/images given mixtures.
% Special cases of interest are source separation and aposteriori coding.
% 
% centered and decorrelated residuals xi along with their predicted
% variances.
%
%   * Mandatory parameters : 
%       - X : FxTxI : STFT or MDCT of the mictures. 
%       - P : FxTxJ : estimated Power Spectral Densities of the sources
%
%
%   * Optional parameters : 
%
%       MIXING PARAMETERS : 
%       - A :   FxIxJ : frequency response of mixing filters. 
%            or IxJ   : instantaneous mixing matrix
%       - R : FxIxIxJ 
%          or 1xIxIxJ : spatial covariances of the source images
%
%           note :  1) Only provide one of these. 
%                   2) Default is A = ones(F,I,J)
%
%       - U : FxIxJ   : frequency response of beamforming filters
%
%
%       ORIGINAL SOURCE SIGNALS : 
%       - S : FxTxJ   : original source signals, default : empty
%       - Y : FxTxIxJ : original image signals, default : empty
%
%       PARALLEL OPTIONS:
%       - parallel :  1 or 0 (default). Beware, using parallel may use a 
%                     LOT of memory, especially for large files
%       - nworkers : -1 if all available workers (defaults), or the number
%                    of workers to use
%
%       OUTPUTS : defines content of varargout
%       - 'Smu'    : whether to output Smu (default = 0)
%       - 'Ymu'    : whether to output Ymu (default = 0)
%       - 'SKLT'   : whether to output SKLT and SLambda (default=0)
%       - 'RPost'  : whether to output Rpost 
%       - 'PPost'  : whether to output Ppost 
%
%
% Output : out is a structure containing the following fields, if selected:
%       - 'Smu'     : FxTxJ : estimated sources 
%       - 'Ymu'     : FxTxIJ : estimated images
%       - 'SKLT'    : FxTxJ  : centered KLT of S|X (for source coding)
%       - 'SL'      : FxTxJ  : variances of SKLT   (for source coding)
%       - 'KLT'     : FxTxJxJ: KLT matrix          (for source coding)
%       - 'Rpost'   : FxIxIxJ : reestimation of R at M-step
%       - 'Ppost'   : FxTxJ   : reestimation of P at M-step
%
%--------------------------------------------------------------------------
% Copyright (C) 2014, Inria, Antoine Liutkus
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Inria nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL INRIA BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


%Initialization
%--------------------------------------------------------------
regularization = 1E-5;
minL = 1E-4;


%Size of signals
[F,T,I] = size(X);
J       = size(P,3);

%Parsing input
p = inputParser;
p.addParamValue('A', [], @(x)isnumeric(x));
p.addParamValue('R', [], @(x)isnumeric(x));
p.addParamValue('U', [], @(x)isnumeric(x));
p.addParamValue('S', [], @(x)isnumeric(x));
p.addParamValue('Y', [], @(x)isnumeric(x));
p.addParamValue('Smu', 0, @(x)isnumeric(x)||islogical(x));
p.addParamValue('Ymu', 0, @(x)isnumeric(x)||islogical(x));
p.addParamValue('SKLT', 0, @(x)isnumeric(x)||islogical(x));
p.addParamValue('Rpost', 0, @(x)isnumeric(x)||islogical(x));
p.addParamValue('Ppost', 0, @(x)isnumeric(x)||islogical(x));
p.addParamValue('parallel', 0, @(x)isnumeric(x)||islogical(x));
p.addParamValue('nworkers', -1, @(x)isnumeric(x));
p.KeepUnmatched = true;
p.parse(varargin{:})

P = P;%+eps;%max(eps,P);

% Checking consistencies
if ~isempty(p.Results.R)&&(p.Results.Smu)&&isempty(p.Results.U)
    disp('   cannot recover sources if mixing is full-rank and U==[]. Aborting.');
    return;
end
if isempty(p.Results.S)&&(p.Results.SKLT)
    disp('   cannot compute mean-removed aposteriori KLT of sources without S. aborting.');
    return;
end
if ~isempty(p.Results.A) && ~isempty(p.Results.R)
    disp('both mixing filters and spatial covariances are set. Aborting.')
    return;
end
if (~p.Results.Smu)&&(~p.Results.Ymu)&&(~p.Results.SKLT)&&(~p.Results.Rpost)
    disp('Nothing to do. Aborting.')
    return;
end

% Building images spatial covariances and mixing parameters
R = zeros(F,I,I,J,'single');
A = p.Results.A;

if ~isempty(A)
    %A is available
    if (size(A,1) == 1) || (numel(size(A))==2)
        % as instantaneous mixing
        A = shiftdim(reshape(repmat(A,1,F),I,J,F),2);
    end
elseif isempty(p.Results.R)
    % A is not available and neither is R
    A = ones(F,I,J);
end


if ~isempty(p.Results.R)
    % R is available
    if size(p.Results.R,1) == 1
        for j = 1:J
            for i1 = 1:I
                for i2 = 1:I
                    R(:,i1,i2,j) = p.Results.R(1,i1,i2,j);
                end
            end
        end
    else
        R = p.Results.R;
    end
elseif ~isempty(A)
    %building R from A
    R = zeros(F,I,I,J,'single');
    for j = 1:J
        for i1=1:I
            for i2 = 1:I
                R(:,i1,i2,j) = A(:,i1,j).*conj(A(:,i2,j));
            end
        end
    end
end


if ~isempty(p.Results.U)
    for j = 1:J
        A(:,:,j) = squeeze(sum( bsxfun(@times, permute(R(:,:,:,j), [1 3 2]), p.Results.U(:,:,j)),2));
    end
end


% Initializing matlab parallel toolbox if available
if p.Results.parallel && (exist('matlabpool'))
    % Start labs if necessary.
    sz = matlabpool('size');
    if (sz ==0)
        if p.Results.nworkers>=0
            matlabpool(p.Results.nworkers)
        else
            matlabpool
        end
    end
    % Check we got some now.
    sz = matlabpool('size');
    if (sz ==0)
        error('Failed to open parallel workers');
    end
end


%Computing mixture covariance matrices using spatial covariances
%as Cxx(f,t,:,:) = sum_j [P(f,t,j)*R(f,:,:)]
Cxx = sum( bsxfun(@times, permute(R,[1 5 2 3 4]),permute(P,[1 2 4 5 3])),5);
for i=1:I
    Cxx(:,:,i,i) = Cxx(:,:,i,i)+regularization;
end

%Computing inverses of mixture covariance matrices
invCxx = zeros(F,T,I,I,'single');
if I==1
    invCxx = 1./Cxx;
elseif I==2
    invDet = 1./(Cxx(:,:,1,1).*Cxx(:,:,2,2) - Cxx(:,:,1,2).*Cxx(:,:,2,1));
    invCxx(:,:,1,1) = invDet.*Cxx(:,:,2,2);
    invCxx(:,:,2,1) = -invDet.*Cxx(:,:,2,1);
    invCxx(:,:,1,2) = -invDet.*Cxx(:,:,1,2);
    invCxx(:,:,2,2) = invDet.*Cxx(:,:,1,1);
else
    %general case : no use of analytical expression (slow!)
    disp('   computing inverses of mix covariance matrices..');
    parfor_progress(F);
    for f=1:F
        parfor_progress;
        for t = 1:T
            invCxx(f,t,:,:)=pinv(squeeze(Cxx(f,t,:,:)));
        end
    end
    parfor_progress(0);
end

if p.Results.SKLT || p.Results.Smu
    %first comput diag[ P(f,t,:)] * A(f,:,:)'
    CssAh = zeros(F,T,J,I,'single');
    for j =1:J
        for i = 1:I
            CssAh(:,:,j,i) = bsxfun(@times, conj(A(:,i,j)),P(:,:,j));
        end
    end
    %then multiply by Cxx(f,t,:,:)^-1  to get Wiener gains
    G = zeros(F,T,J,I,'single');
    for j=1:J
        for i = 1:I
            G(:,:,j,i) = sum( bsxfun(@times, squeeze(CssAh(:,:,j,:)), squeeze(invCxx(:,:,:,i))),3);
        end
    end
    
    % Computing a posteriori mean and variance
    Smu = zeros(F,T,J,'single');
    for j = 1:J
        Smu(:,:,j) = sum( X .* squeeze(G(:,:,j,:)),3);
    end
    
    if p.Results.SKLT
        disp('   computing KLT and posterior variances..');
        SL = zeros(F,T,J,'single');
        SKLT = zeros(F,T,J,'single');
        KLT  = zeros(F,T,J,J,'single');
        parfor_progress(F); % Initialize

        for f=1:F
            parfor_progress; % Count
            Af = reshape(squeeze(A(f,:,:)),I,J);
            for t=1:T
                Cpost = diag(squeeze(P(f,t,:))) - squeeze(G(f,t,:,:))*Af*diag(squeeze(P(f,t,:)));
                [KLT(f,t,:,:), Lc, ~] = svd(Cpost);
                SL(f,t,:) = max(minL, diag(Lc));
                SKLT(f,t,:) = squeeze(KLT(f,t,:,:))'*(squeeze(p.Results.S(f,t,:))-squeeze(Smu(f,t,:)));
            end
        end
        parfor_progress(0); % Clean up

        fprintf('\n')
    end
end

if p.Results.Ymu || p.Results.Rpost || p.Results.Ppost
    Ymu = zeros(F,T,I,J,'single');
    if p.Results.Rpost
        Rpost = zeros(F,I,I,J,'single');
    end
    if p.Results.Ppost
        Ppost = zeros(F,T,J,'single');
    end
    
    
    parfor j=1:J
        Ymuj = zeros(F,T,I,'single');
        Rj = R(:,:,:,j);
        
        disp(sprintf('   separating image %d / %d and computing appropriate information..',j,J));
        G = zeros(F,T,I,I,'single');
        Cyx = bsxfun(@times,permute(Rj(:,:,:),[1 5 2 3 4]),P(:,:,j));

        for i1=1:I
            for i2=1:I
                G(:,:,i1,i2) = sum( bsxfun(@times, squeeze(Cyx(:,:,i1,:)), squeeze(invCxx(:,:,:,i2))),3);
            end
        end       
        
        for i1 = 1:I
            Ymuj(:,:,i1) = sum( bsxfun(@times,squeeze( G(:,:,i1,:)), X),3);
        end
        
        if p.Results.Rpost || p.Results.Ppost
            %First compute a posteriori covariance of images (E-step)
            KY = zeros(F,T,I,I,'single');
            for i1 = 1:I
                for i2=1:I
                    KY(:,:,i1,i2) = Ymuj(:,:,i1).*conj(Ymuj(:,:,i2));
                end
            end
            
            
            IminG = bsxfun(@minus,permute(eye(I),[3 4 1 2]), G);
            IminGCyx = zeros(F,T,I,I,'single');
            for i1 = 1 :I
                for i2 = 1:I
                    IminGCyx(:,:,i1,i2) = ...
                        sum( bsxfun(@times, squeeze(IminG(:,:,i1,:)), squeeze(Cyx(:,:,:,i2))),3);
                end
            end
            KY = KY + IminGCyx;

            %We need update of spatial covariance matrix
            % first compute inverse of a Rj
            disp('   computing inverses of spatial covariance matrices..');
            if I==1
                invRj = 1./Rj;
            elseif I==2
                invDet = 1./(Rj(:,1,1).*Rj(:,2,2) - Rj(:,1,2).*Rj(:,2,1));
                invRj(:,1,1) =  invDet.*Rj(:,2,2);
                invRj(:,2,1) = -invDet.*Rj(:,2,1);
                invRj(:,1,2) = -invDet.*Rj(:,1,2);
                invRj(:,2,2) =  invDet.*Rj(:,1,1);
            else
                %general case : no use of analytical expression (slow!)
                invRj = zeros(F,I,I,'single');
                if size(p.Results.R,1) == 1
                    invRjInst = pinv(squeeze(Rj(1,:,:)));
                    for i1 = 1:I
                        for i2 = 1:I
                            invRj(:,i1,i2) = invRjInst(i1,i2);
                        end
                    end
                else
                    for f=1:F
                        invRj(f,:,:)=pinv(squeeze(Rj(f,:,:)));
                    end
                end
            end

            %Computing 1/I*tr(invRj*KY)
            Ppost = zeros(size(P));
            for i=1:I
                Ppost(:,:,j) = Ppost(:,:,j) + real(sum(bsxfun(@times,invRj(:,i,:),KY(:,:,:,i)),3));
            end
            Ppost(:,:,j) = Ppost(:,:,j)/I;


            if p.Results.Rpost
                %We need Rpost : compute reestimate of R (M-step)
                Rpost(:,:,:,j) = 1/T*squeeze(sum( bsxfun(@times, KY, 1./Ppost(:,:,j)),2));
            end
        end        
        Ymu(:,:,:,j) = Ymuj;        
    end
end
   
if p.Results.Smu
    out.Smu = Smu;
end
if p.Results.Ymu
    out.Ymu = Ymu;
end
if p.Results.SKLT
    out.SL   = SL;
    out.SKLT = SKLT;
    out.KLT  = KLT;
end
if p.Results.Rpost
    out.Rpost = Rpost;
end
if p.Results.Ppost
    out.Ppost = Ppost;
end