function K = periodicKernel(period,P, V,H)
%Builds a periodic kernel, looking like this:
%10001000100010001

temp = strel('periodicline',P,[0, period]);
L = 2*P*period+1;
K = zeros(2*V+1,L);
K(V+1,:)=temp.getnhood;
K(:,P*period+1-H:P*period+1+H)=1;
