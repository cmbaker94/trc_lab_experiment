function calc_directional_spec_wg(Tinfo,geom,frange)
% This code will compute the directional spectra from the wire resistance
% gauges as the bridge (i.e., offshore wave conditions). This code is
% designed to compute the spectra from a pre-defined geometery (see options
% listed below) with three different methods (SSE array approach) as well
% as compute the autospectra at each sensor.
% INPUT:
% Tinfo: Trial info
% geom : geometry number, typically 1 or 2
% frange: range of frequencies to compute the directional spectra over

% geom = 1;
% frange = [0.25 1.2];% freq to calc Sd over

addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
% addpath(genpath('E:\codes\insitu'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))
addpath(genpath('E:\code\wafo'))
%
% %% NOTES
% I tried calculating:
% WL = 256 and 128 -- not sure what is best
% method: MLM, IMLM, EMEP
% geom 1, 2 (similar but 2 is more centered for 10 deg?)
% subsampling timeseries at 10 hz or using full 100 hz  (seem similar so
% use 10 Hz)

%% STEP 0: Pre-defined typical input

smoothz = 1;
inHz = 100;
array =  1; %surfzone array  wg: 1, inner shelf array wg: 2

windowlength = [256/2 256];
degreenum = [181 361];

%% STEP 1:  Get trail details

Tinfo = trial_files(Tinfo);
% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);
% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);
[camera.time] = cam_time(Tinfo);
insHz = 100;

%% STEP 2: Load data

[~,wg] = load_insitu(Tinfo,array);
IDdef.depth = Tinfo.tide;
z = wg.z;
XL       = wg.xyz(1,:);
YL       = wg.xyz(2,:);
ZL       = mean(wg.z,1);%press.xyz(3,:);

%% STEP 3: STEP 4: Select Points

% pick sensors:
if geom == 1
    sensorpick = [4  12  13  10  6 3 1]; %  BEST SO FAR
elseif geom == 2
    sensorpick = [2 4  5 7 11  12  14]; %  BEST SO FAR 2?
elseif geom == 3
    sensorpick = [4 12 11 15 5 2 9]; % try this one next
elseif geom == 4
    sensorpick = [4  12  13  10  6 3 1 15]; %  BEST SO FAR
elseif geom == 5
    sensorpick = [2 4  5 7 11  12  14  9]; %  BEST SO FAR
end

for ip = 1:length(sensorpick)
    sensorname = ['wg',num2str(sensorpick(ip),'%02.f')];
    ix(ip)  = find(strcmp(sensorname,wg.loc));
end

%% STEP 4: Prepare data

xd = XL(ix);
yd = YL(ix);
zd = ZL(ix);

tpass = 1/10;

if smoothz == 1
    for ip = 1:length(ix)
        ztemp = detrend(pl64ta(z(:,ix(ip)),insHz*tpass));
        Z(:,ip) = ztemp;
        %     ID.data         = z(:,ix(ip)); % 6000 x 14
    end
    Z = squeeze(Z(1:1/tpass:end,:));
    Time = [0:1:size(Z,1)-1]/(1/tpass);
else
    Z  = z(:,ix);
    Time = [0:1:size(Z,1)-1]/insHz;
end

%% STEP 5: Compute pwelch properties
for isel = 1:length(sensorpick)
    sensorpickautospec = sensorpick(isel); % try this one next
    sensorname = ['wg',num2str(sensorpickautospec,'%02.f')];
    ixauto  = find(strcmp(sensorname,wg.loc));
    ztemp = detrend(pl64ta(z(:,ix(ip)),insHz*tpass));
    if smoothz == 0
        [Spoint.spec(:,isel),Spoint.freq]   = pwelch(ztemp,insHz*2^6,insHz*2^5,[],insHz,'ConfidenceLevel',0.95); % compute spectra
        [Spoint.Hs(isel),Spoint.Tp(isel)] = spec2HsTp(Spoint.freq,Spoint.spec(:,isel),frange);
    elseif smoothz == 1
        [Spoint.spec(:,isel),Spoint.freq]   = pwelch(ztemp,(1/tpass)*2^6,(1/tpass)*2^5,[],(1/tpass),'ConfidenceLevel',0.95); % compute spectra
        [Spoint.Hs(isel),Spoint.Tp(isel)] = spec2HsTp(Spoint.freq,Spoint.spec(:,isel),frange);
    end
end

%% STEP 6: Prepare  to compute the spectra

H = Tinfo.tide;
N = length(xd);
types = repmat(sensortypeid('n'),N,1);
bfs   = ones(N,1);
pos   = [xd(:),yd(:),zeros(N,1)];
% Z     = Z';

%% STEP 7: Loop through multile options for windows and degree-resolution

for iWL = 1:2
    WL = windowlength(iWL);
    for idegnum = 1:length(degreenum)
        degnum = degreenum(idegnum);
        
        for imethod = 1:3
            if imethod == 1
                dirmethod = 'MLM';
            elseif imethod == 2
                dirmethod = 'IMLM';
            elseif imethod == 3
                dirmethod = 'EMEM';
            end
            
            if smoothz == 0
                Spec    = dat2dspec([Time' Z],[pos types,bfs],H,insHz*WL,degnum,dirmethod);
            elseif smoothz == 1
                Spec    = dat2dspec([Time' Z],[pos types,bfs],H,(1/tpass)*WL,degnum,dirmethod);
            end
            
            spec = Spec.S'/9;
            dirs = Spec.theta*(180/pi);
            freqs = Spec.w/(2*pi);
            dtheta        = abs(dirs(2)-dirs(1));
            df            = abs(freqs(2)-freqs(1));
            
            [~,ilow] = nanmin(abs(freqs-frange(1)));
            [~,ihigh] = nanmin(abs(freqs-frange(2)));
            Sf      = nansum(spec,2)*dtheta;
            for f = 1:size(spec,2)
                Sd(f) = trapz(freqs(ilow:ihigh),spec(ilow:ihigh,f));
            end
            
            [Hs,Tp] = spec2HsTp(freqs,Sf,frange);
            
            %% find param
            [dirsprd_fit,Sfit,dirang] = fit_cos2s(Tinfo,dirs'+180,abs(Sd));
            %      [dirsprd(iT),Sfit(iT),dirang(iT)] = fit_cos2s_angsprd(Tinfo,dirs,abs(Sd(:,iT))');
            %  this one has a different set of iterations
            rrange = [-pi/2 pi/2];
            rad = deg2rad(dirs-180);
            [val,idmin] = min(abs(rad-rrange(1)));
            [val,idmax] = min(abs(rad-rrange(2)));
            rad = rad(idmin:idmax);
            Srad = Sd(idmin:idmax);
            dirsprd_obs = calc_spread(rad,abs(Srad)');
            
            %% restructure data  and  store
            
            makestruct.spec = spec;
            makestruct.freqs = freqs;
            makestruct.dirs = dirs;
            makestruct.Sd = Sd;
            makestruct.Sf  = Sf;
            makestruct.Hs  = Hs;
            makestruct.Tp =  Tp;
            makestruct.sprd = dirsprd_fit;
            makestruct.dirang  = dirang;
            makestruct.sprd_obs = dirsprd_obs;
            
            eval([dirmethod ,' = makestruct;'])
            clear spec freqs dirs Sd Sf Hs Tp dirsprd_fit dirang dirsprd_obs
            
        end
        
        %% Try DIWASP
        
        fsmall        = [1:0.01:3];
        fsmall4       = 10^-2*fsmall.^-4;
        
        figure('units','inches','position',[1 1 10 8],'color','w');
        subplot(211)
        semilogy(IMLM.freqs,IMLM.Sf,'k','linewidth',2);
        hold on
        semilogy(MLM.freqs,MLM.Sf,'r','linewidth',2);
        semilogy(EMEM.freqs,EMEM.Sf,'b','linewidth',2);
        semilogy(Spoint.freq,Spoint.spec,'k','linewidth',1,'LineStyle','-.');
        semilogy(fsmall,fsmall4,'b','linewidth',2);
        h2 = legend('MLM','IMLM','EMEM','Single pt.')
        set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
        grid
        xlabel('$\mathrm{f}~\mathrm{(Hz)}$','interpreter','latex')
        ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
        ylim([10^-4.5 10^-0.4])
        xlim([0 3]);
        h1=gca;
        set(h1,'fontsize',15);
        
        subplot(212)
        semilogy(MLM.dirs,MLM.Sd,'k','linewidth',2);
        hold on
        semilogy(IMLM.dirs,IMLM.Sd,'r','linewidth',2);
        semilogy(EMEM.dirs,EMEM.Sd,'b','linewidth',2);
        grid
        xlabel('$\theta~\mathrm{(^{\circ})}$','interpreter','latex')
        ylabel('$S_{\theta}~(\mathrm{m^{2}/degree})$','interpreter','latex')
        xlim([-180 180]);
        ylim([10^(-6) 10^(-3)])
        h1=gca;
        set(h1,'fontsize',15);
        
        Sname1 = [Tinfo.figfolder,'dirspec_wg_geom',num2str(geom),'_dir',num2str(diff(dirs(1:2))),'_WL',num2str(WL)];
        print(Sname1,'-dpng')
        
        psname = [Tinfo.savefolder,'dirspec_wg_geom',num2str(geom),'_dir',num2str(diff(dirs(1:2))),'_WL',num2str(WL),'.mat'];
        eval(['save -v7.3 ',psname,' sensorpick',' MLM',' IMLM',' MEM',' WL',' degnum',' Spoint'])
        clear MLM IMLM MEM WL degnum dirmethod
    end
end