function [eta] = press2sse_timeseries(press,Fs,offset,maxfac,filtramp)
% This function will conver pressure gauge instrument records to
% auto-spectra of eta based on the elevation of the instrument.
% pressure -> depth
% p = ph + pd
% where the hydrostatic pressure is ph
% and the dynamics pressure is pd
% here:
% ph = -\rho*g*z        
% pd = \rho*g*\eta*Kp;
% Kp = cosh(k*(h+z))/cosh(k*h);
% where Kp is the pressure response factor,
%   z is the elevation of the inst. below MSSE, 
%   h is the depth of the water from MSSE, and
%   k is the wavenumber
% therefore, 
% eta = ((p/(\rho*g))+z)/Kp;
%%%%%%
% INPUT:
% press = structure of raw pressure gage measurments
% % % z     = elevation of sensor relative to surface (positive is up)
% % % h     = depth of water
% Fs    = frequency of data
% filtramp = default us 1000!!!!1
% OUT:
% sse   = timeseries of the sea-surface elevation
 
% given
g       = 9.81;
rho     = 1000;
 
% check to make sure z is correct
% if z>0
%     error('Issue: Depth of instrument, z should be negative (positive is upward from surface)')
% end
 
% Begin computing
etaKp       = (press/(rho*g));
% find mean water depth
hp          = nanmean(etaKp);
% subtract MSSE to get etaKp
etaKp       = etaKp - hp;
% Total depth
h           = hp + offset;
% compute spectra
% % [SeeKp,f,SeeKpc]   = pwelch(etaKp,WL*Hz,OL*Hz,[],Hz,'ConfidenceLevel',0.95);
% [SeeKp,f,SeeKpc]   = pwelch(etaKp,WL,OL,[],Hz,'ConfidenceLevel',0.95);
 
NFFT = length(etaKp);
SeeKp = fft(etaKp,NFFT);
% F = ((0:1/NFFT:1/2-1/NFFT)*Fs).';
dt = 1/Fs;
fmin = -1/(2*dt);
df = 1/(NFFT*dt);
f0 = -fmin;
f = (mod(linspace(0, 2*f0-df, NFFT)+f0,  2*f0)  - f0)';
 
% find wavenumbers
omega       = 2*pi*f;
omega2      = omega.^2;
const       = omega2*h/g;
kh          = dispersi(const);
k           = kh./h;
% compute pressure response factor
Kp          = (cosh(k*(h-hp))./cosh(kh)).^2;
% add in maxfactor
% Kpp         = ones(size(Kp));
% ifreq       = find(f<maxfac);
% Kpp(ifreq)  = Kp(ifreq);
% ifreq       = find(-f<maxfac);
% Kpp(ifreq)  = Kp(ifreq);
 
% compute expected auto-spectra of sea-surface elevation
Seetemp         = SeeKp./Kp;
% generate filter
 
ihalf = round(length(Kp)/2);
ffilt           =  f(1:filtramp:ihalf);%f(1:1000:ihalf);
[temp, iloc]    = nanmin(abs(ffilt-maxfac));
filtb           = ones(1,length(ffilt));
rampval       = [0.7089 0.2363 0.0225];
filtb(iloc-1:iloc+1) = rampval;
filtb(iloc+2:end) = 0;
 
filtbn = interp(filtb,1000);
filtbn = filtbn(1:ihalf);
filtbn(f(1:ihalf)>1.5) = 0;
filtbn(f(1:ihalf)<0.85) = 1;
filtbn = movmean(filtbn,200);
filtbn = movmean(filtbn,100);
 
filt = [filtbn fliplr(filtbn)]';
filt = filt(1:length(etaKp));
 
% [temp, iloc] = nanmin(abs(f(1:4800)-maxfac));
% filtb         = ones(1,240);
% ilocb          = round(iloc/10);
% rampval       = [0.7089 0.2363 0.0225];
% filtb(ilocb-1:ilocb+1) = rampval;
% filtb(ilocb+2:end) = 0;
 
% See(f>1.2)  = NaN;
% See(f<-1.2) = NaN;
% See         = NaN(size(Seetemp));
% See(ifreq)  = Seetemp(ifreq);
 
See = Seetemp.*filt;
See(isnan(See)) = 0;
 
eta = ifft(See,NFFT,'symmetric');

