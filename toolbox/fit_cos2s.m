function [dirsprd,Sfit,dirang] = fit_cos2s(dir,Sd)
% Fit a cos-2s distribution: This code will fit a directional spectra with
% the closest cos-2s distribution and output the directional spread and
% angle. The code fits the mean angle based on the max 20 values of the
% spectra and the directional spread based on the RMSE between the
% cos-2s distribution and input spectra, Sd.
%
% INPUT: 
% dir:  directions for the directional spectra [deg] 
%       (code expects that spectra is centered around 180 deg)
% Sd:   spectra as a function of degrees
% OUTPUT:
% dirsprd:  directional spread estimated [deg]
% Sfit:     cos-2s distrubtion
% dirang:   mean direction [deg]
%
% Code created by C.M. Baker, Sept 2021

rrange = [-pi/2 pi/2]; 
rad = deg2rad(dir-180);
[val,idmin] = min(abs(rad-rrange(1)));
[val,idmax] = min(abs(rad-rrange(2)));
rad = rad(idmin:idmax);
Srad = Sd(idmin:idmax);
fascend =sort (Srad,2, 'ascend');
floor = mean(fascend(1:30));
Srad = Srad-floor;

% normalize input spectra
intS = trapz(rad,Srad);
Srnorm = Srad/intS;

% angle
[~,ith] = maxk(Srad,20);
dir = mean(rad(ith));

S = 1:85; % full range of s distributions in increments of 1
TH = dir;
pltflg = 0;
for i = 1:length(S)
    sdist = ((2.^(2*S(i)-1))/pi).*((gamma(S(i)+1)).^2)./gamma(2*S(i)+1).*cos((rad-TH)/2).^(2*S(i));
    y0 = sdist;
    y1 = Srnorm;
%     mean_err(i) = mean(abs(y0-y1)) ;    % mean absolute error
%     MSE(i) = mean((y0-y1).^2) ;        % Mean square error
    RMSE(i) = sqrt(mean((y0-y1).^2)) ; % Root mean square error
end
[val,iS] = nanmin(RMSE);
Sfit = S(iS);
Sdist = ((2.^(2*Sfit-1))/pi).*((gamma(Sfit+1)).^2)./gamma(2*Sfit+1).*cos((rad-TH)/2).^(2*Sfit);

% plot input spectra and closes fit spectra
figure
plot(rad,Srnorm)
hold on
plot(rad,Sdist)

% compute directional spread
left = trapz(rad,Sdist.*cos(rad))^2;
right = trapz(rad,Sdist.*sin(rad))^2;
sprd = sqrt(2*(1-sqrt(left+right)));

% convert spread and mean angle from radians to degrees
dirsprd = rad2deg(sprd);
dirang = rad2deg(dir);