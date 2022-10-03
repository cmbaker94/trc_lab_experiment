function [sprd] = calc_spread(rad,Srad)
%  Compute the directional spread

fascend =sort (Srad,2, 'ascend');
floor = mean(fascend(1:length(Srad)/12));
Srad = Srad-floor;

% normalize input spectra
intS = trapz(rad,Srad);
Srnorm = Srad/intS;%
left = trapz(rad,abs(Srnorm).*cos(rad))^2;
right = trapz(rad,abs(Srnorm).*sin(rad))^2;
sprdrad = sqrt(2*(1-sqrt(left+right)));
sprd = rad2deg(sprdrad);