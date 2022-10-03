function [S,f,theta,k] = MLEDirSpec(Y,X,Nfft,dt,dx,a,d)

%MLEDirSpec - Maximum likelihood estimate for dir., freq., and wavenumber spectrum
%
%USAGE: [S f theta k] = MLEDirSpec(Y,X,Nfft,dt,dx,a)
%
%INPUT:
%  Y - is observation data, MxN array (equally spaced down columns)
%  X - is Nx2 obs. locations (x,y)
%  Nfft - is the  length of the signal to ensemble
%    average (with 50% overlap)
%    NOTE: additional smoothing with a convoltuion filter is also inplemented.
%    The filter length is calculated so number of total averages for each frequency
%    (including ensemble and conv. averages) is at least equal to N (number of instruments).
%  dt - is sampling interval
%  dx - uniformly spaced fake spatial sampling interval (units of length) (if dx is a vector
%    desired wavenumber estimation locations are assumed, in units of 1/length)
%  a - is desired wave angle estimation locations in radians (uniformly spaced, if scalar a
%    assumed sampling interval -pi:a:pi)
%OUTPUT:
%  S - is the returned MLE spectrum
%  f - is frequency
%  theta - is angle in radians (usual r.h. angle convention, as for sin(a))
%  k - output wavenumber vector

% C. Chickadel, Jan 2011, chickadel@apl.washington.edu


% take care of theta
if prod(size(a))==1
  theta = [-pi:a:pi];
  da = a;
else
  theta = a(:)';
  da = mean(abs(diff(theta)));
end

% make lag matricies
xL = repmat(X(:,1),1,length(X(:,1))) - repmat(X(:,1)',length(X(:,1)),1);
yL = repmat(X(:,2),1,length(X(:,2))) - repmat(X(:,2)',length(X(:,2)),1);
L = sqrt(xL.^2 + yL.^2);

% make k vector
if prod(size(dx))==1
  dk = 1/max(abs(L(:)));
  %k = [dk:dk:1/(2*dx)];
  %k = 1/max(L(:)):range([1/max(L(:)) 0.5/min(L(~~L(:)))])/24:0.5/min(L(~~L(:))); % fix to 25 points from
   k = [dk:dk:1/(2*min(L(L(:)>0)))]; % use fourier determined k's (based on max and min lags)
else
  k = dx;
  dk = mean(abs(diff(k)));
end

% size of data
[N M] = size(Y);

% make f, keep only positive frequencies
T = Nfft*dt;
f =  [0:-1:-floor(Nfft/2) (ceil(Nfft/2)-1):-1:1]'/T;
frqInd = find(f > 0);
f = f(frqInd);
Nf = length(f);

% make ensemble averages of cross-spectral matrix
Nb = floor((N-(Nfft-round(Nfft/2)))/round(Nfft/2));
Q = zeros(M,M,Nf);
taper = repmat(bartlett(Nfft),1,M);
for j = 0:Nb-1
  Yt = Y((1+j*round(Nfft/2)):(j*round(Nfft/2)+Nfft),:);
  Yfft = fft(taper.*(Yt - repmat(mean(Yt),Nfft,1))); % should demean and taper time series at least
  for jj = 1:Nf
    Qtemp(:,:,jj) = Yfft(frqInd(jj),:)'*Yfft(frqInd(jj),:);
  end
  Q = Q + Qtemp/Nb;
end

% now conv. smoother to to decreace variance in estimates
Nw = round(1.25*M-Nb);
Nw = Nw+1-mod(Nw,2);  % ensure odd length window
Nw = max(Nw,3);       % set minimum window length
smWin = bartlett(Nw);
smWin = smWin/sum(smWin);
Q = shiftdim(Q,2);
QMean = repmat(mean(Q,1),[Nf 1 1]);
Q = convn(Q-QMean,smWin,'same'); % beware ends will taper due to smoothing, try to rectify
Q = shiftdim(Q+QMean,1);         % by removing mean first and adding it in later

% make lag, k and theta matricies to compute S
[cosTheta XL K] = meshgrid(cos(theta),xL(:),k);
[sinTheta YL K] = meshgrid(sin(theta),yL(:),k);
P = exp(i*2*pi*(K.*cosTheta.*XL + K.*sinTheta.*YL));
Nk = length(k);
Ntheta = length(theta);

% loop along frequencies to get spec. estimates at all the k-theta choices
warning off MATLAB:divideByZero
for j = 1:Nf
  qInv = pinv(Q(:,:,j));
  QInv = repmat(qInv(:),[1 Ntheta Nk]);
  Stemp = 1./sum(P.*QInv);
  S(j,:,:) = real(Stemp)/(sum(abs(Stemp(:)))*da*dk*4*pi^2)*mean(diag(Q(:,:,j)));
end
warning on

%key TSA waves
%comment MLE directional frequency spectral estimator
