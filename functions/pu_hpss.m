% Source separation procedure using phase unwrapping for initializing the
% harmonic part
%
% Ref:
% "Model-based STFT phase recovery for audio source separation",
% Paul Magron, Roland Badeau and Bertrand David,
% IEEE Transactions on Speech, Audio, and Language Processing, May 2018
%
% Inputs:
%     X : F*T mixture
%     Ve : F*T*J magnitudes of the sources
%     hop : STFT overlap (in samples)
%     Nit : number of iterations
%
% Outputs:
%     Se : estimated sources

function [Se] = pu_hpss(X,Ve,hop,Nit)

% Parameters
[F,T,J] = size(Ve);
Nfft = (F-1)*2;
Ct = 2*pi*hop/Nfft;

% Weights
G = Ve.^2./(repmat(sum(Ve.^2,3),[1 1 J])+eps);

% Init
Se = Ve .* exp(1i * repmat(angle(X),[1 1 J]));

% Loop over time frames
for t=2:T
    
    % Current frame
    St = squeeze(Se(:,t,:));
    Gt = squeeze(G(:,t,:));
    Vt = abs(St);
    Xt = X(:,t);

    % PU initial estimate for harmonic source
    f_inf = get_frequencies_qifft_frame(Vt(:,2)+eps);
    phiaux = angle(Se(:,t-1,2))+Ct*f_inf;
    St(:,2) = Vt(:,2) .* exp(1i * phiaux);

    % Iter procedure
    for iter =1:Nit
        E = Xt - sum(St,2);
        Y = St + repmat(E,[1 J]).*Gt;
        St = Y ./ (abs(Y)+eps) .* Vt;
    end

    % Update the sources in frame t
    Se(:,t,:) = St;
        
end

end


% PU: frequencies and regions of influence
function [f_inf,f_centr,f_harm] = get_frequencies_qifft_frame(v)

v = v(:)';

%Central peaks
%[~,f_centr] = findpeaks(v,'MINPEAKHEIGHT',0.01*max(v));
[~,f_centr] = findpeaks(max(v,10^-6),'MINPEAKHEIGHT',max(v)*10^-2);

Nfreq = length(f_centr);
f_harm = zeros(1,Nfreq);

if (Nfreq >0)
    % QIFFT
    for ind = 1:Nfreq
        f = f_centr(ind);
        f_harm(ind) = qint(log10(v(f-1)),log10(v(f)),log10(v(f+1)))+f;
    end

    % Regions of influence
    f_inf = zeros(length(v),1);
    deb = 1;
    index_lim = zeros(1,Nfreq-1);
    
    for ind = 1:(Nfreq-1)
        f = f_centr(ind);
        fp = f_centr(ind+1);
        fin = floor((v(fp)*f+v(f)*fp)/(v(fp)+v(f)));
        f_inf(deb:fin) = f_harm(ind);
        deb = fin+1;
        index_lim(ind) = fin;
    end

    f_inf(deb:end) = f_harm(end)-1;
    
else
    f_inf = (1:length(v))'-1;
end

end

% Quadratic Interpolated FFT
function [p,b,a] = qint(ym1,y0,yp1)

p = (yp1 - ym1)/(2*(2*y0 - yp1 - ym1));
b = y0 - 0.25*(ym1-yp1)*p;
a = 0.5*(ym1 - 2*y0 + yp1);

end