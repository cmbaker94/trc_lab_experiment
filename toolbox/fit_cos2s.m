function [dirsprd,Sfit,dirang] = fit_cos2s(Tinfo,dir,Sd,varargin)
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
Srnorm = Srad/intS;%

% angle
if nargin > 3
    dir = deg2rad(varargin{1});
else
    if Tinfo.spread == 0
        [~,ith] = maxk(Srad,3);
        c = Srad(ith);
    else
        [~,ith] = maxk(Srad,round(length(Srad)/12));%30
        c = Srad(ith);
    end
    
    dir = mean(rad(ith));
%     dircentroid = trapz(Srnorm(ith).*rad(ith),rad(ith))/trapz(Srnorm(ith),rad(ith));
%     dir = dircentroid;
%     display('using direction as the centroid')
end

S = 1:0.05:171; % full range of s distributions in increments of 1
TH = dir;
pltflg = 0;
if Tinfo.spread >  45
    for i = 1:length(S)
        %     sdist = ((2.^(2*S(i)-1))/pi).*((gamma(S(i)+1)).^2)./gamma(2*S(i)+1).*cos((rad-TH)/2).^(2*S(i));
        % Could  try  this  verion
        %     test(i) = gamma(S(i));
        %     test(i) = gamma(S(i)+1)/2/sqrt(pi)./gamma(S(i)+1/2);
        sdist = gamma(S(i)+1)/2/sqrt(pi)./gamma(S(i)+1/2).*cos((rad-TH)/2).^(2*S(i));
        y0 = sdist;
        y1 = Srnorm;
        %     mean_err(i) = mean(abs(y0-y1)) ;    % mean absolute error
        %     MSE(i) = mean((y0-y1).^2) ;        % Mean square error
        RMSE(i) = sqrt(mean((y0-y1).^2)) ; % Root mean square error
    end
    [val,iS] = nanmin(RMSE);
    Sfit = S(iS);
    % Sdist = ((2.^(2*Sfit-1))/pi).*((gamma(Sfit+1)).^2)./gamma(2*Sfit+1).*cos((rad-TH)/2).^(2*Sfit);
    Sdist = gamma(Sfit+1)/2/sqrt(pi)./gamma(Sfit+1/2).*cos((rad-TH)/2).^(2*Sfit);
elseif Tinfo.spread <45
%     c = max(Srnorm);
%     fitfun = fittype( @(s,x) c*cos((x-TH)/2).^(2*s));
    fitfun = fittype( @(c,s,x) c*cos((x-TH)/2).^(2*s));
    [fitted_curve,gof] = fit(rad',Srnorm',fitfun);%,'StartPoint',x0)
%     [fitted_curve,gof] = fit(rad',Srad',fitfun)
    coeffvals = coeffvalues(fitted_curve);
    Sfit = coeffvals(2);
    coeffvals(1) = 1/trapz(rad,cos((rad-TH)/2).^(2*coeffvals(2)));
    Sdist = coeffvals(1)*cos((rad-TH)/2).^(2*coeffvals(2));
%     Sdist = Sdist/trapz(rad,Sdist);
%     Sfit = coeffvals(1);
%     Sdist = c*cos((rad-TH)/2).^(2*coeffvals(1));
end

% plot input spectra and closes fit spectra
figure('units','inches','position',[1 1 5 5],'color','w')
plot(rad2deg(rad),Srnorm)
hold on
plot(rad2deg(rad),Sdist,'k')
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',1*get(h1,'ticklength'));
set(h1,'fontsize',20);
xlabel('$\theta~(^{\circ})$','interpreter','latex','fontsize',20);
ylabel('$S(\theta)$ (m$^{2}/^{\circ}$)','interpreter','latex','fontsize',20);
% sname =['ex_sprd'];
%     print([Tinfo.figfolder,sname],'-dpng')

% compute directional spread
left = trapz(rad,Sdist.*cos(rad))^2;
right = trapz(rad,Sdist.*sin(rad))^2;
sprd = sqrt(2*(1-sqrt(left+right)));

% left = trapz(rad,Srad.*cos(rad))^2;
% right = trapz(rad,Srad.*sin(rad))^2;
% sprd = sqrt(2*(1-sqrt(left+right)));

% convert spread and mean angle from radians to degrees
dirsprd = rad2deg(sprd);
dirang = rad2deg(dir);