classdef Signal < handle
%-------------------------------------------------------------------------
% Class name : Signal
%-------------------------------------------------------------------------
% Description : handles basic Signal operations and transforms. 
%               * wave loading
%               * STFT, with any overlap ratio, any frames length,
%                 any weighting function
%               * MDCT and its inverse
%               * Constant Q Transform (by Jacques Prado)
%               * splitting into frames
%               * pitch detection (for harmonic sounds)
%               * onset detection
%
% Main properties that are read/write 
%   * s : signal
%   * windowLength (ms)
%   * nfft (samples)
%   * overlapRatio (>=0 and <1)
%   * S : stft data
%   * Smdct : MDCT data
%   * CQTBinsPerOctave : number of bins per octave for cQt
%   * CQTMinFreq : min freq for the cQt
%   * CQTMaxFreq : max freq for the cQt
%   * CQTAlign : where to center cQt lobe weights, either 'l' for left or
%               'c' for center
%   * weightingFunction: handle to a function taking the number of points
%     N as an argument and returning the Nx1 weighting window
%
% Main properties that are read only : 
%   * sLength : signal length
%   * nChans : number of channels
%   * nfftUtil : number of bins in the positive frequency domain
%   * framesPositions, nFrames : positions and number of frames
%   * sWin, sWeights : windowed data
%   * Sq : cQt data
%
% examples : 
%   sig = Signal('myfile.wav'); %creates signal for a wave file
%   sig.windowLength = 70;      %sets windows of 70ms
%   sig.overlapRatio = 0.8;     %sets overlap ratio 0 < overlapRatio < 1
%   sig.STFT;                   %performs STFT
%   sig.CQT;                    %same for CQT
%   sig.MDCT;                   %same for MDCT (note that overlap is set to
%                               %               50 percents)
%   mySTFT = sig.S;
%   myMDCT = sig.S;
%   myCQT = sig.Sq
%
%   %modifying s.S for example
%   ...
%
%   waveform = sig.iSTFT();     % gets inverse transform
%
% Note that :
%   * all properties are automatically set relevantly in case of
%     modifications. for example, when nfft is set, windowLength is changed
%     accordingly
%   * STFT, MDCT and CQT produce exactly aligned data
%-------------------------------------------------------------------------
%
% Copyright (c) 2014, Antoine Liutkus, Inria
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions 
% are met:
%
%    *  Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%    *  Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the 
%       distribution.
%    *  Neither the name of Inria nor the names of its 
%       contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%-------------------------------------------------------------------------
    methods 
        %Constructor
        function self = Signal(varargin)
            %signal properties
            self.singlePrecision = 0;
            self.s = [];
            self.fs = 44100;
            self.sLength = 0;
            self.nChans = 0;
            self.weightingFunction = @hamming;
            
            %STFT properties
            self.S = [];
            self.windowLength = 60;
            self.nfft = 0;
            self.nfftUtil = 0;
            self.overlapRatio = 0.5;
            self.framesPositions = [];
            self.nFrames = 0;
            self.weightingWindow = [];
            self.overlap = 0;
            
            %MDCT properties
            self.Smdct = [];
            
            %CQT properties
            self.CQTKernelComputationNeeded = 1;
            self.CQTSparseKernel = [];
            self.CQTBinsPerOctave = 12;
            self.CQTMinFreq=97.9999;
            self.CQTMaxFreq = self.fs/2;
            self.CQTAlign = 'c';

            %Windowing properties
            self.sWin = [];
            self.sWeights = [];
            
            %Properties listeners
            self.propertiesListeners = {};
            self.propertiesListeners{end+1} = addlistener(self,'s','PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));
            self.propertiesListeners{end+1} = addlistener(self,'singlePrecision','PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));
            self.propertiesListeners{end+1} = addlistener(self,'fs','PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));
            self.propertiesListeners{end+1} = addlistener(self,'weightingFunction','PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));
            self.propertiesListeners{end+1} = addlistener(self,'windowLength' ,'PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));
            self.propertiesListeners{end+1} = addlistener(self,'nfft' ,'PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));
            self.propertiesListeners{end+1} = addlistener(self,'overlapRatio','PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));
            self.propertiesListeners{end+1} = addlistener(self,'S','PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));    
            self.propertiesListeners{end+1} = addlistener(self,'CQTBinsPerOctave','PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));    
            self.propertiesListeners{end+1} = addlistener(self,'CQTMinFreq','PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));    
            self.propertiesListeners{end+1} = addlistener(self,'CQTMaxFreq','PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));    
            self.propertiesListeners{end+1} = addlistener(self,'CQTAlign','PostSet',@(src,evnt)Signal.handlePropertyEvents(self,src,evnt));    
            
            switch nargin
                case 0
                    return;
                case 1
                    switch class(varargin{1})
                        case 'char'
                            %varargin{1} is a file
                            self.LoadFromFile(varargin{1});
                        case 'Signal'
                            self = varargin{1}.copy;
                    end
                case 2
                    %varargin{1} is a signal and varargin{2} is the
                    %sampling frequency
                    self.LoadWF(varargin{1},varargin{2});

            end
        end
        
        function SigOut = getOnsets(self,fMin, fMax)
            %This function returns a vector of length nFrames, containing
            % an onset detection in the signal between frequencies fMin and
            % fMax. It uses the complex spectral difference method.
            if isempty(self.S)
                self.STFT
            end
            
            IFmin = max(1,round(fMin/self.fs*self.nfft));
            IFmax = round(min(fMax/2,fMax/self.fs*self.nfft));
            
            SigFrame = self.S(IFmin:IFmax,:,:);
            
            LenWindow = size(SigFrame,1);
            nCanaux = size(SigFrame,3);
            Nbre = size(SigFrame,2);
            SigOut = zeros(size(SigFrame,2),nCanaux);
            
            for canal = 1:nCanaux
                SigOutInst = zeros(LenWindow,1);
                CurrPhi = angle(SigFrame(:,1,canal));
                SigOut(1) = 0;
                Phi = (unwrap(angle(SigFrame(:,:,canal))));
                DevPhi = [zeros(LenWindow,2) (Phi(:,3:end) - 2*Phi(:,2:end-1) +...
                        Phi(:,1:end-2))];
                    for n = 2:Nbre
                            SigOutInst = sqrt(abs(SigFrame(:,n-1)).^2 + abs(SigFrame(:,n)).^2 -2*abs(SigFrame(:,n)).*abs(SigFrame(:,n-1)).*cos(DevPhi(:, n)));
                            SigOut(n,canal) = sum(SigOutInst);
                    end
                SigOut(:,canal) = 1/max(abs(SigOut(:,canal))).*SigOut(:,canal);
            end
        end
        
        function [pitchs, pitchCandidates, harmonicTemplates, pitchScores]= mainPitch(self, fMin, fMax,pitchStep,tolerance)
            %  [pitchs, pitchCandidates, harmonicTemplates, pitchScores]= mainPitch(self, fMin, fMax,pitchStep,tolerance)
            %  Computes dominant pitch trajectory between fMin and fMax.
            %  tolerance indicates in Hz the allowed variation of the main
            %  pitch from one frame to the next.
            % Returns:
            % * pitchs: a vector giving the detected pitch in Hz for each
            % frame.
            % * pitchCandidates: pitch frequencies used as candidates (Cx1
            %   vector)
            % * harmonicTemplates: template spectrum for each candidate
            %   (FxC) matrix
            % * pitchScores: matrix indicating the score for each candidate
            %   and each frame. (CxT) matrix
            if nargin<5
                tolerance=10;
            end
            if nargin<4
                pitchStep=1;
            end
            if isempty(self.S)
                self.STFT();
            end
            pitchCandidates=fMin:pitchStep:fMax;
            harmonicTemplates=self.KLGLOTTModel(pitchCandidates).^2;
            nPitchs = length(pitchCandidates);
            filterWidth = round(tolerance/pitchStep);
            pitchScores = zeros(nPitchs, self.nFrames);
            prevScore= misc.normalize(ones(nPitchs,1));
            for n = 1:self.nFrames
                if rem(n,round(self.nFrames/10)) == 1
                    disp(sprintf('pitch detection: %d%%',round(n/self.nFrames*100)));
                    drawnow;
                end
                likelihoods = misc.normalize(exp(harmonicTemplates'*(abs(self.S(:,n,1)).^2)));
                pitchScores(:,n) = misc.normalize(misc.normalize(misc.smooth(prevScore,filterWidth)).*likelihoods);
                prevScore = pitchScores(:,n);
            end
            [~, argPitch] = max(pitchScores,[],1);
            pitchs = pitchCandidates(argPitch);
        end
        
        %Data loaders
        function LoadFromFile(self, file)
            [self.s, self.fs] = audioread(file);
            [self.sLength, self.nChans] = size(self.s);
        end 
        function LoadWF(self, waveform, fs)
            self.s = waveform;
            self.fs = fs;
            [self.sLength, self.nChans] = size(self.s);
        end
        
        %Transforms
        function S = STFT(self)
            %Builds Signal STFT, puts it in self.S and returns it.
            weighted_frames = self.split(0);
            self.S = fft(weighted_frames);  
            if self.singlePrecision
                self.S = single(self.S);
            end
            if nargout
                S = self.S;
            end
        end
        
        function V = specgram(self,nSteps)
            %Builds Signal spectrogram, 
            % through averaging nSteps delayed spectrograms
            if nargin == 1
                nSteps = 20;
            end
            weighted_frames = self.split(0);
            V = abs(fft(weighted_frames)).^2;  
            for index = 1:nSteps
                disp(sprintf('specgram, pass %d / %d',index,nSteps));drawnow
                weighted_frames = self.split(index);
                V = V*(1- 1/(index+1)) +1/(index+1)* abs(fft(weighted_frames)).^2;  
            end
            if self.singlePrecision
                V = single(self.S);
            end
        end
        
        function s = iSTFT(self)
            %Computes inverse STFT, puts it in self.s and returns it.
            s = [];
            if isempty(self.S)
                return;
            end
            tempSignal = real(ifft(self.S));
    
            tempLength = self.framesPositions(end)+self.nfft - 1;     
            oldLength = self.sLength;       
            
            s_temp=zeros(tempLength,self.nChans);
            W = zeros(tempLength,self.nChans);

            for index_frame=1:self.nFrames
                    for index_chan=1:self.nChans
                        s_temp(self.framesPositions(index_frame):(self.framesPositions(index_frame)+self.nfft-1),index_chan)= ...
                            s_temp(self.framesPositions(index_frame):(self.framesPositions(index_frame)+self.nfft-1),index_chan) + tempSignal(:,index_frame,index_chan).*self.weightingWindow;
                        W(self.framesPositions(index_frame):(self.framesPositions(index_frame)+self.nfft-1),index_chan) = ...
                            W(self.framesPositions(index_frame):(self.framesPositions(index_frame)+self.nfft-1),index_chan) + self.weightingWindow.^2;
                    end
            end
            s_temp = (s_temp./W);      
            self.s = s_temp(1:oldLength,:);

            if nargout
                s = self.s;
            end
        end
        
        function Smdct = MDCT(self)
            self.overlapRatio = 0.5;
            self.nfft = round(self.nfft/2)*2;
            N = self.nfft;
            J = N/2;
            self.suspendListeners();
            self.s = [zeros(N,self.nChans); self.s; zeros(2*N,self.nChans)];
            self.activateListeners();
            K = floor(size(self.s,1)/J)-1;
            persistent G
            if (isempty(G)) || (size(G, 2) ~= N)
                G = zeros(J,N);
                win = sqrt(Signal.hann_nonzero(N));
                for i = 0:J-1,
                     G(i+1,:) = win.*(sqrt(2/J)*cos((2*i+1)*(2*(0:2*J-1)-J+1)*pi/(4*J)));   
                end
            end
            k = [1:K]';
            indices = ((k-1)*J*ones(1,2*J) + ones(K,1)*[1:2*J])';
            for chan = 1:self.nChans
                temp = self.s(:,chan); 
                self.Smdct(:,:,chan) = G*temp(indices);
            end
            self.nFrames = K;
            if nargout
                Smdct = self.Smdct;
            end
            self.suspendListeners();
            self.s = self.s((self.nfft+1):(self.nfft+self.sLength),:);
            self.activateListeners();
        end
        
        
   function s = iMDCT(self)
            %Computes inverse MDCT, puts it in self.s and returns it.
            s = [];
            if isempty(self.Smdct)
                return;
            end
            persistent invG
            J = self.nfft/2;
            if (isempty(invG)==1) || (size(invG, 1) ~= self.nfft)
                invG = zeros(J,self.nfft);
                win = sqrt(Signal.hann_nonzero(self.nfft));
                for i = 0:J-1,
                    invG(i+1,:) = win.*(sqrt(2/J)*cos((2*i+1)*(2*(0:2*J-1)-J+1)*pi/(4*J)));
                end
                invG = invG';
            end
            signal = zeros(J*(self.nFrames+1), self.nChans);
            for chan = 1:self.nChans
                for t = 1:self.nFrames
                    T = reshape(invG*self.Smdct(:,t,chan),self.nfft/2,2);
                    T2 = [T(:,1:2:end),zeros(J,1)]+[zeros(J,1),T(:,2:2:end)];
                    signal_frame = T2(:);                    
                    indx_frame = (1:self.nfft) + (self.nfft/2)*(t-1);
                    signal(indx_frame,chan) = signal(indx_frame,chan) + signal_frame;                    
                end;
            end
            if ~isempty(self.sLength)
                signal = signal((self.nfft+1):(self.nfft+self.sLength),:);
            end;

            self.s = signal;
            if nargout
                s = self.s;
            end
        end        
        
        
        function s = unsplit(self,tempSignal)
            tempLength = self.framesPositions(end)+self.nfft - 1;     
            oldLength = self.sLength;       
            
            s_temp=zeros(tempLength,self.nChans);
            W = zeros(tempLength,self.nChans);

            for index_frame=1:self.nFrames
                    for index_chan=1:self.nChans
                        s_temp(self.framesPositions(index_frame):(self.framesPositions(index_frame)+self.nfft-1),index_chan)= ...
                            s_temp(self.framesPositions(index_frame):(self.framesPositions(index_frame)+self.nfft-1),index_chan) + tempSignal(:,index_frame,index_chan).*self.weightingWindow;
                        W(self.framesPositions(index_frame):(self.framesPositions(index_frame)+self.nfft-1),index_chan) = ...
                            W(self.framesPositions(index_frame):(self.framesPositions(index_frame)+self.nfft-1),index_chan) + self.weightingWindow.^2;
                    end
            end
            s_temp = (s_temp./W);            
            s = s_temp(1:oldLength,:);
        end
        
        function Sq = CQT(self)
                if self.CQTKernelComputationNeeded
                    self.ComputeSparseKernel;
                end
                nfftKernel = size(self.CQTSparseKernel,1);

                %Compute positions and number of STFT frames
                stepSTFT = self.nfft-self.overlap;
                framesPositionsSTFT = 1:stepSTFT:self.sLength;
                framesPositionsSTFT((framesPositionsSTFT + self.nfft -1) > self.sLength ) = [];
                framesPositionsSTFT(end+1) = framesPositionsSTFT(end) + self.nfft;
                nFramesSTFT = length(framesPositionsSTFT);                
                framesCentersSTFT = framesPositionsSTFT + self.nfft/2;
                

                %Compute CQT frames start so that frames centers are the
                %same than for the STFT
                framesPositionsCQT = round(framesCentersSTFT - nfftKernel/2);
                tempSignal = self.s;
                if framesPositionsCQT(1) < 1
                    tempSignal = [zeros(-framesPositionsCQT(1) + 1, self.nChans) ; tempSignal];
                    framesPositionsCQT = framesPositionsCQT - framesPositionsCQT(1) + 1;
                end
                lEnd = (framesPositionsCQT(end) + nfftKernel - 1) - length(tempSignal);
                finalZeroPadding = zeros(lEnd, self.nChans);
                tempSignal = [tempSignal ; finalZeroPadding];
                win = hamming(nfftKernel);
                win = win(:);
                
                %Compute CQT
                self.Sq = zeros(size(self.CQTSparseKernel,2),nFramesSTFT,self.nChans);
                disp('Computing CQT...')
                for index=1:nFramesSTFT
                    if mod(index, 100) == 1
                        disp(sprintf('     CQT : %0.2f percents',index/nFramesSTFT*100));
                    end                        
                    for index_chan = 1:self.nChans
                        self.Sq(:,index,index_chan) = self.CQTSparseKernel' ...
                                                    * fft(tempSignal(framesPositionsCQT(index):framesPositionsCQT(index)+nfftKernel-1,index_chan).*win);
                    end
                end
                disp('finished');
                
                if nargout 
                    Sq = self.Sq;
                end
        end
        
        

        
        function [sWin, sWeights] = split(self,store,initialPos)
            %default value
            if (nargin <= 1)
                store = 1;
            end
            if nargin <= 2
                initialPos = 1;
            end
            
            %computing splitting parameters
            step = self.nfft-self.overlap;
            self.framesPositions = initialPos:step:self.sLength;
            self.framesPositions((self.framesPositions + self.nfft -1) > self.sLength ) = [];
            self.framesPositions(end+1) = self.framesPositions(end) + step;
            self.nFrames = length(self.framesPositions);
            
            %initializing splitted signal matrix
            sWin = zeros(self.nfft,self.nFrames, self.nChans);
            
            %Initializing weighting window
            win = self.weightingWindow(:)*ones(1,self.nChans);

            %Initializing weighting signal
            tempLength = self.framesPositions(end)+self.nfft - 1;     
            sWeights = zeros(tempLength,1);
            
            %For each frame, compute signal
            for index=1:(self.nFrames - 1)
                sWin(:,index,:) = self.s(self.framesPositions(index):self.framesPositions(index)+self.nfft-1,:).*win;
                sWeights(self.framesPositions(index):(self.framesPositions(index)+self.nfft-1)) = ...
                            sWeights(self.framesPositions(index):(self.framesPositions(index)+self.nfft-1)) + self.weightingWindow;                
            end
            
            %Handle last frame on which zeropadding may be necessary
            lEnd = (self.framesPositions(end) + self.nfft - 1) - self.sLength;
            finalZeroPadding = zeros(lEnd, self.nChans);
            sWin(:,end,:) = [self.s(self.framesPositions(end):end,:) ; finalZeroPadding].*win;
            
            %handle last frame for weighting signal
            sWeights(self.framesPositions(end):(self.framesPositions(end)+self.nfft-1)) = ...
                        sWeights(self.framesPositions(end):(self.framesPositions(end)+self.nfft-1)) + self.weightingWindow;                

            %stores result if necessary
            if store 
                self.sWin = sWin;
                self.sWeights = sWeights;
            end
            if nargout == 0
                clear sWin
                clear sWeights
            end
                
        end
        
        function Scomp = buildComplete(self, S)
            %builds complete spectrogram/STFT/mask using a version up to nfft/2
            Scomp = zeros(size(self.S));
            for c = 1:self.nChans
                Scomp(1:self.nfftUtil,:,c) = S(:,:,c);
                Scomp(end:-1:(end-self.nfftUtil+2),:,c) = conj(S(2:end,:,c));
            end
        end
        
        %copy
        function copied = copy(self)
            copied = Signal;
            copied.suspendListeners
            Meta = ?Signal;
            for index = 1:length(Meta.Properties)
                if strcmp(Meta.Properties{index}.Name, 'propertiesListeners')
                    continue
                end
                eval(['copied.' (Meta.Properties{index}.Name) '= eval([''self.'' (Meta.Properties{index}.Name)]);']);
            end
            copied.activateListeners
        end
        
        function stripCanal(self, channel)
            self.suspendListeners;
            toKeep = setdiff(1:self.nChans,channel);
            self.s = self.s(:,toKeep);
            if ~isempty(self.S)
                self.S = self.S(:,:,toKeep);
            end
            if ~isempty(self.Sq)
                self.Sq = self.Sq(:,:,toKeep);
            end
            if ~isempty(self.sWin)
                self.sWin = self.sWin(:,:,toKeep);
            end
            self.nChans = length(toKeep);
            self.activateListeners;
        end
    end
        
    properties (Access=public, SetObservable=true)
        %Waveform, sampling frequency
        s, fs, singlePrecision
        
        %weighting parameters
        weightingFunction
        
        %STFT parameters and data
        windowLength, nfft, overlapRatio
        S
        
        %CQT parameters
        CQTBinsPerOctave, CQTMinFreq, CQTMaxFreq, CQTAlign
        
        %MDCT
        Smdct
    end
    
    properties (SetAccess = private, GetAccess = public, SetObservable=false)
        %waveform properties
        sLength, nChans, nfftUtil

        %STFT parameters and data
        framesPositions, nFrames
        
        %Windowed data
        sWin, sWeights
        
        %CQT data
        Sq
    end
    properties (Access = private)
        %Properties listeners
        propertiesListeners
        
        %STFT 
        weightingWindow, overlap
        
        %CQT
        CQTKernelComputationNeeded, CQTSparseKernel
    end
    methods (Static)
        function handlePropertyEvents(self,src,evnt)
            switch src.Name 
                case 'singlePrecision'
                    self.S = single(self.S);
                    self.Sq = single(self.Sq);
                case 's' 
                    [self.sLength, self.nChans] = size(self.s);
                case 'fs' 
                    self.nfft = round(self.fs*self.windowLength/1000);
                    self.nfftUtil = round(self.nfft/2);
                    self.overlap = round(self.nfft*self.overlapRatio); 
                    self.CQTKernelComputationNeeded = 1;
                case 'windowLength' 
                    self.nfft = round(self.fs*self.windowLength/1000);
                    self.nfftUtil = round(self.nfft/2);
                    self.overlap = round(self.nfft*self.overlapRatio); 
                    self.CQTKernelComputationNeeded = 1;
                case 'nfft' 
                    self.windowLength = (self.nfft*1000/self.fs);
                    self.overlap = round(self.nfft*self.overlapRatio); 
                    self.CQTKernelComputationNeeded = 1;
                case 'overlapRatio'
                    self.overlap = round(self.nfft*self.overlapRatio); 
                case 'S'
                    self.nFrames = size(self.S,2);
                    self.nfft = size(self.S,1);
                    self.nfftUtil = round(self.nfft/2);
                    self.nChans = size(self.S,3);
                case 'CQTBinsPerOctave'
                    self.CQTKernelComputationNeeded = 1;
                case 'CQTMinFreq'
                    self.CQTKernelComputationNeeded = 1;
                case 'CQTMaxFreq'
                    self.CQTKernelComputationNeeded = 1;
                case 'CQTAlign'
                    self.CQTKernelComputationNeeded = 1;
            end
        
            if self.nfft
                nArgWeightingFunc = nargin(self.weightingFunction);
                if (nArgWeightingFunc<0)||(nArgWeightingFunc==1)
                    self.weightingWindow = feval(self.weightingFunction, self.nfft);
                elseif nArgWeightingFunc==2
                    self.weightingWindow = feval(self.weightingFunction, self.nfft,self.overlapRatio);
                else
                    error('weighting function badly defined.');
                end
            end
        end  
        function w = hann_nonzero(N)
            n = 0:(N-1);
            phi0 = pi / N;
            Delta = 2 * phi0;
            w = 0.5 * (1 - cos( n*Delta + phi0));
        end
    end
    methods (Access = private)
        function suspendListeners(self)
            for index = 1:length(self.propertiesListeners)
                self.propertiesListeners{index}.Enabled = false;
            end
        end
        function activateListeners(self)
            for index = 1:length(self.propertiesListeners)
                self.propertiesListeners{index}.Enabled = true;
            end
        end
       function harmonicTemplates = KLGLOTTModel(self,pitchs)
            nPitchs = length(pitchs);
            harmonicTemplates = zeros(self.nfft,nPitchs);
            AV = 1;
            Oq = 0.05;
            for index = 1:nPitchs
                %build glottal pulse pattern
                T = round(self.fs/pitchs(index));
                HH = zeros(T,1);
                t = 1:round(T*Oq);
                HH(t) = ((27*AV)/(4*T*Oq^2)*(t-1).^2-(27*AV)/(4*T^2*Oq^3)*(t-1).^3);
                nRepeat = ceil(self.nfft/T);
                H = repmat(HH, nRepeat,1);
                H = H(1:self.nfft);
                H = H - mean(H);
                F = abs(fft(hann(self.nfft).*H,self.nfft));
                harmonicTemplates(:,index) = F(:)/sum(abs(F));
            end      
        end        
        function ComputeSparseKernel(self)
            disp(sprintf('compute CQT kernel on freqs %0.2f to %0.2f with %i bins per octave ...',self.CQTMinFreq,self.CQTMaxFreq,self.CQTBinsPerOctave));
            thresh = 0.0075;
            Q = 1 / (2^(1/self.CQTBinsPerOctave)-1);
            K = ceil(self.CQTBinsPerOctave*log2(self.CQTMaxFreq/self.CQTMinFreq));
            fftLen = 2^nextpow2(ceil(Q*self.fs/self.CQTMinFreq));
            self.CQTSparseKernel = [];
            im = sqrt(-1);
            tempKernel = zeros(fftLen,1);
            for k = K:-1:1;
                len = ceil( Q * self.fs / (self.CQTMinFreq*2^((k-1)/self.CQTBinsPerOctave)));
                switch self.CQTAlign
                    case 'l'
                        tempKernel =[hamming(len,'periodic')/sum(hamming(len,'periodic')); zeros(fftLen-len,1)];
                        tempKernel = tempKernel.*exp(2*pi*im*Q*(0:fftLen-1)'/len);
                    case 'c'
                        if rem(len,2)==1
                            len=len-1;
                        end;
                        index = fftLen/2 - len/2;
                        tempKernel =[zeros(index,1) ; 
                                    hamming(len,'periodic')/sum(hamming(len,'periodic')) ;...
                                    zeros(index,1)];
                        tempKernel = tempKernel.*exp(2*pi*im*Q*(0:fftLen-1)'/len);
                    case 'r'
                        tempKernel =[zeros(fftLen-len,1) ;...
                                hamming(len,'periodic')/sum(hamming(len,'periodic'))];
                        tempKernel = tempKernel.*exp(2*pi*im*Q*(0:fftLen-1)'/len);
                end
                specKernel = fft(tempKernel);
                specKernel(find(abs(specKernel)<=thresh))=0;
                self.CQTSparseKernel = sparse([specKernel self.CQTSparseKernel]);
            end
            self.CQTSparseKernel = conj(self.CQTSparseKernel)/fftLen;
            self.CQTKernelComputationNeeded = 0;
        end
    end
end
