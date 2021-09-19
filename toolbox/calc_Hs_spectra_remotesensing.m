function [Hs,Tp,See,f,Seec,nanratio] = calc_Hs_spectra_remotesensing(z,x,y,WL,OL,range,Hz,nancutoff,plotspec)
% This function will compute the spectra for all locations in a time stack.
% INPUT:
% z     = 3D matrix of elevation
% x     = array of x-location in tank coordinates
% y     = array of y-location in tank coordinates
% WL    = window length for pwelch (in sec)
% OL    = overlap length for pwelch (in sec)
% range = range of frequencies to compute Hs over
% Hz    = frequency of data
% nancutoff = cut off where not compute Hs
% plotspec = boolean to indicate if want to plot spectras
% OUTPUT:
% Hs    = 2D matrix of Hs computed from frequency range of spectra
% Tp    = 2D matrix of Tp chosen within frequency range of spectra
% See   = 3D matrix of spectra at all locations
% f     = 3D matrix of frequency of spectra for all locations
% Seec  = 3D confidence interval for spectra
% nanratio = 2D ratio of nans in timeseries

eta = NaN.*z;
Hs  = squeeze(eta(:,:,1));
Tp  = Hs;
etatemp = randi([0, 1], [1, size(eta,3)])-0.5;
[Stemp,ftemp,Sctemp]   = pwelch(etatemp,WL*Hz,OL*Hz,[],Hz,'ConfidenceLevel',0.95);
See = NaN(size(z,1),size(z,2),length(Stemp));
f = NaN(size(z,1),size(z,2),length(ftemp));
Seec = NaN(size(z,1),size(z,2),length(Stemp),2);

if plotspec == 1
    figure('units','inches','position',[1 1 12 5],'Color','w');
end

for i = 1:size(z,1)
    for j = 1:size(z,2)
        ztemp           = squeeze(z(i,j,:));    % grab z for location
        etatemp         = ztemp-nanmean(ztemp);     % compute eta for loc
        % remove nans and interpolate
        nanz            = isnan(ztemp);             % find location of NaNs
        nant            = [1:numel(ztemp)]';       % make strings for interp
        nanratio(i,j)   = sum(nanz)/length(nant); % find ratio of NaNs
        eta(i,j,:)      = etatemp;              % storing eta
        if nanratio(i,j) <= nancutoff
            ztemp(nanz)    = interp1(nant(~nanz),ztemp(~nanz),nant(nanz)); % interpolate between the NaNs
            zcutnan   = ztemp(~isnan(ztemp)); % if nans still at the beginning or end of list, just remove these points
            zmean     = nanmean(zcutnan);
            etacut   = zcutnan-zmean;
            eta4spec = detrend(etacut);
%             trend = etacut - eta4spec;
            [Stemp,ftemp,Sctemp]   = pwelch(eta4spec,WL*Hz,OL*Hz,[],Hz,'ConfidenceLevel',0.95); % compute spectra
            
            if plotspec == 1
                clf
                semilogy(ftemp,Stemp,'k','LineWidth',2)
                xlabel('$f$ (Hz)','interpreter','latex','fontsize',20);
                ylabel('$S_{\eta\eta}$ (m$^2$/Hz)','interpreter','latex','fontsize',20);
                h1=gca;
                set(h1, 'YScale', 'log')
                set(h1,'tickdir','in','xminortick','on','yminortick','on');
                set(h1,'ticklength',1*get(h1,'ticklength'));
                set(h1,'fontsize',15);
                xlim([0 2.5])
                ylim([10^-5 10^-1])
                title(['x = ',num2str(x(j)),'m , y = ',num2str(y(i))],'interpreter','latex','fontsize',20);
                pause(0.2)
            end
            
            [Hs(i,j),Tp(i,j)]   = spec2HsTp(ftemp,Stemp,range); % compute Hs, Tp from spectra
            if size(Stemp,1) == size(See,3)
                See(i,j,:)      = Stemp;
                f(i,j,:)        = ftemp;
                Seec(i,j,:,:)   = Sctemp;
            else
                See(i,j,:)      = padarray(Stemp,size(See,3),NaN);
                f(i,j,:)        = padarray(ftemp,size(See,3),NaN);
                Seec(i,j,:,2)   = padarray(Sctemp,[size(See,3), 2],NaN);
                display(['Cut ',len(isnan(etatemp)),' NaNs from string and frequency band is',len(ftemp)])
            end
        elseif nanratio(i,j) > 0.05
            continue
        end
    end
end