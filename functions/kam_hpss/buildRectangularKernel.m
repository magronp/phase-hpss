function kernel = buildRectangularKernel(shape,nfft,fs,hopSize)
nF = max(1,round(shape(1)*nfft/fs));
nT = max(1,round(shape(2)/hopSize));
kernel=ones(nF,nT);