function [See,f,Seec] = press2spectra(press,Hz,WL,OL,rho,offset,maxfac)
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
% Hz    = frequency of data
% WL    = window length in sec
% OL    = overlap length in sec
% rho   = density of water
% OUT:
% See   = auto-spectra of eta
% f     = frequency

% given
g       = 9.81;

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
[SeeKp,f,SeeKpc]   = pwelch(etaKp,WL,OL,[],Hz,'ConfidenceLevel',0.95);
% find wavenumbers
omega       = 2*pi*f;
omega2      = omega.^2;
const       = omega2*h/g;
kh          = dispersi(const);
k           = kh./h;
% compute pressure response factor
Kp          = (cosh(k*(h-hp))./cosh(kh)).^2;
% add in maxfactor
Kpp         = ones(size(Kp));
ifreq       = find(f<maxfac);
Kpp(ifreq)  = Kp(ifreq);
% compute expected auto-spectra of sea-surface elevation
Seetemp     = SeeKp./Kpp;
Seectemp    = SeeKpc./Kpp;
See         = NaN(size(Seetemp));
Seec        = NaN(size(Seectemp));
See(ifreq)  = Seetemp(ifreq);
Seec(ifreq,:) = Seectemp(ifreq,:);

