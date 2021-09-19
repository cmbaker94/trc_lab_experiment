function [dirsprd, Sfit,dirang] = fit_cos2s(dir,Sd);

rrange = [-pi/2 pi/2]; 
rad = deg2rad(dir-180);
[val,idmin] = min(abs(rad-rrange(1)));
[val,idmax] = min(abs(rad-rrange(2)));
rad = rad(idmin:idmax);
Srad = Sd(idmin:idmax);
fascend =sort (Srad,2, 'ascend');
floor = mean(fascend(1:30));
Srad = Srad-floor;

intS = trapz(rad,Srad);
Srnorm = Srad/intS;

% angle
[~,ith] = maxk(Srad,20);
dir = mean(rad(ith));

S = 1:85;
TH = dir;
pltflg = 0;
for i = 1:length(S)
    sdist = ((2.^(2*S(i)-1))/pi).*((gamma(S(i)+1)).^2)./gamma(2*S(i)+1).*cos((rad-TH)/2).^(2*S(i));
%     serror(i) = fit1dregerror(rad,sdist,rad,Srnorm,pltflg);
    y0 = sdist;
    y1 = Srnorm;
%     abs_dy(i,:) = abs(y0-y1);   % absolute error
%     relerr(i,:) = abs(y0-y1)./y0 ;  % relative error
%     pererr(i,:) = abs(y0-y1)./y0*100 ;   % percentage error
    mean_err(i) = mean(abs(y0-y1)) ;    % mean absolute error
    MSE(i) = mean((y0-y1).^2) ;        % Mean square error
    RMSE(i) = sqrt(mean((y0-y1).^2)) ; % Root mean square error
end
[val,iS] = nanmin(RMSE);
Sfit = S(iS);
Sdist = ((2.^(2*Sfit-1))/pi).*((gamma(Sfit+1)).^2)./gamma(2*Sfit+1).*cos((rad-TH)/2).^(2*Sfit);


figure
plot(rad,Srnorm)
hold on
plot(rad,Sdist)

pick = 4;
left = trapz(rad,Sdist.*cos(rad))^2;
right = trapz(rad,Sdist.*sin(rad))^2;
sprd = sqrt(2*(1-sqrt(left+right)));
dirsprd = rad2deg(sprd);
dirang = rad2deg(dir);