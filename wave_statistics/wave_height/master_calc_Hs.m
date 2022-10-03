% Compute sea-surface elevation spectra, Hs, and Tp for in situ, stereo,
% and lidar

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\cameras'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

spread = [20];

for i = 1:length(spread)
    
% Trial info
Tinfo.Hs        = 0.30;
Tinfo.Tp        = 2;
Tinfo.tide      = 1.07;
Tinfo.spread    = spread(i);

% spectral analysis details
sdet.WL              = 2^5;
sdet.OL              = 2^4;
sdet.frange     = [0.25, 1.2];
sdet.nancutoff  = 0.1;

%% STEP 1: Create paths, files and naming

% cam = TRC_camera_info(Tinfo);
% Tinfo.cam = cam;
% Tinfo.cam.tstart = cam.cam.tstart;
% Tinfo.cam.tdate = cam.cam.tdate;

% % general path and names
datapath    = 'E:\';
% 
Tinfo = trial_files(Tinfo);

% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);
% % Stereo Reconstructions
% Tinfo.cam.camerasys   = 'TRM';                     % camera setup - e.g., TRM (offshore) or TRC (onshore)
% Tinfo.cam.scene       = '1';                       % scene number of trial - typ 1
% Tinfo.cam.imagestart  = 7200;                      % images number of first frame on file 
% Tinfo.cam.numframes   = 4800;                      % number of frames processed
% Tinfo.cam.dx          = 0.5;;
% Tinfo.cam.regy        = [-13:Tinfo.cam.dy:13];
% Tinfo.cam.Hz          = 8;
% 
% Tinfo.cam.trialname   = [Tinfo.cam.cameras
% Tinfo.cam.dy          = 0.25;
% Tinfo.cam.regx        = [25:Tinfo.cam.dx:37]ys,'-',Tinfo.cam.tdate];
% Tinfo.cam.imagerange  = [num2str(Tinfo.cam.imagestart,'%05.f'),'-',num2str(Tinfo.cam.imagestart+(Tinfo.cam.numframes-1),'%05.f')];
% Tinfo.cam.trimname    = ['frames_',Tinfo.cam.imagerange,'\'];  
% Tinfo.cam.datafolder  = [datapath,'data\processed\DEM\',Tinfo.cam.trialname,'_Scene',Tinfo.cam.scene,'\',Tinfo.cam.trimname];
% 
% % stereo record start and end time
% starttemp       = datenum(Tinfo.cam.tstart(1:end-3),'mm-dd-yyyy-HHMM')+datenum(0,0,0,0,0,Tinfo.cam.imagestart\Tinfo.cam.Hz);
% endtemp         = datenum(Tinfo.cam.tstart(1:end-3),'mm-dd-yyyy-HHMM')+datenum(0,0,0,0,0,(Tinfo.cam.imagestart+Tinfo.cam.numframes-1)\Tinfo.cam.Hz);
% Tinfo.cam.time  = starttemp:datenum(0,0,0,0,0,1\Tinfo.cam.Hz):endtemp;
% clear *temp

% % Bathymetry
% addpath([datapath,'codes\analytic'])
% % load trc bathy
bathymetry = load('E:\data\processed\lidar\Riegl\TRC_bathymetry_estimate_line.mat');
% SWAN0 = load('\Users\cmbaker9\Documents\Research\Lab_Experiments\codes\analytic\SWAN_output_0_degree_spread.mat');
% SWAN40 = load('\Users\cmbaker9\Documents\Research\Lab_Experiments\codes\analytic\SWAN_output_40_degree_spread.mat');

%% STEP 2: Create figure folders

% % figure folder
% fssubfolder = datestr(date,'yy-mm-dd');
% figfolder   = [datapath,'figures\meas_comp\',Tinfo.cam.trialname,'\',Tinfo.cam.trimname,fssubfolder,'\'];
% 
% % make figure folders
% eval(['!mkdir ',datapath,'figures\meas_comp\',Tinfo.cam.trialname]);
% eval(['!mkdir ',datapath,'figures\meas_comp\',Tinfo.cam.trialname,'\',Tinfo.cam.trimname]);
% eval(['!mkdir ',figfolder])

%% STEP 3: Create folders to store data

% display('storing based on condition')
% Tinfo.comp = ['Hs',num2str(Tinfo.Hs*100),'_Tp',num2str(Tinfo.Tp),'_tide',num2str(Tinfo.tide*100),'_spread',num2str(Tinfo.spread)];
% datarange = [datestr(Tinfo.cam.timevec(1)); datestr(Tinfo.cam.timevec(end))];
% datarange = ['time_', datarange(1,13:14), datarange(1,16:17), '-' , datarange(2,13:14), datarange(2,16:17)];
% procpath = [datapath,'data\processed\conditions\',Tinfo.comp];
% 
% eval(['!mkdir ',datapath,'data\processed\conditions\'])
% eval(['!mkdir ',procpath])
% eval(['!mkdir ',procpath,'\',Tinfo.cam.tstart])
% eval(['!mkdir ',procpath,'\',Tinfo.cam.tstart,'\',datarange])

%% STEP 2: Stereo Reconstruction + VeloLidar: compute spectra
sdet.plotspec = 0;

% % stereo
% sdet.Hz     = 8;
% cam         = calc_Hs_stereo(Tinfo,sdet);
% sname       = 'stereo_Hs_spectra';
% eval(['save -v7.3 ',Tinfo.savefolder,sname,' cam Tinfo sdet'])
% 
% %%
% % lidar
% sdet.Hz     = 10;
% lidar       = calc_Hs_lidar(Tinfo,sdet);
% sname       = 'velodyne_Hs_spectra';
% eval(['save -v7.3 ',Tinfo.savefolder,sname,' lidar Tinfo sdet'])

%% Insitu

sdet.Hz         = 100;
sdet.maxfac     = 1.2;

time            = [Tinfo.cam.timevec(1) Tinfo.cam.timevec(end)+datenum(0,0,0,0,0,1/Tinfo.cam.Hz)]-datenum(Tinfo.cam.tstart(1:end-3),'mm-dd-yyyy-HHMM');

szvsis          = 0; % 0 for sz array
[sz.Hs,sz.Tp,sz.See,sz.Seec,sz.f,sz.inst,sz.xyz,sz.Tinfo,sz.data] = calc_Hs_spectra_insitu(Tinfo,time,sdet.WL,sdet.OL,sdet.frange,sdet.maxfac,bathymetry,szvsis); % ADD BACK IN TIME

%%
szvsis          = 1; % 1 for is array and want to grab time range from sz array
[is.Hs,is.Tp,is.See,is.Seec,is.f,is.inst,is.xyz,is.Tinfo,sz.data] = calc_Hs_spectra_insitu(Tinfo,time,sdet.WL,sdet.OL,sdet.frange,sdet.maxfac,bathymetry,szvsis);

% create summary table
sz.sumtable = table(sz.inst',sz.xyz,sz.Hs',sz.Tp');
is.sumtable = table(is.inst',is.xyz,is.Hs',is.Tp');
sz.sumtable.Properties.VariableNames = {'Instrument' 'xyz' 'Hs' 'Tp'};
is.sumtable.Properties.VariableNames = {'Instrument' 'xyz' 'Hs' 'Tp'};

% sz.Hs(find(sz.Hs<0.05)) = NaN;
% is.Hs(find(is.Hs<0.05)) = NaN;

sname = 'insitu_Hs_spectra';
eval(['save -v7.3 ',Tinfo.savefolder,sname,' sz',' is',' sdet'])

end

% %% Load bathymetry and analytic solution
% 
% analytic.h = sz.Tinfo.stillwat-bathymetry.h;
% analytic.x = bathymetry.xp;
% 
% % wave conditions, gamma
% analytic.T = 2; % period, s
% analytic.H0 = 0.3; % offshore wave height, m
% analytic.theta0 = 0; % offshore wave direction
% analytic.gamma = 0.8; % breaking wave gamma
% 
% % run the wave shoaling\breaking code
% wave = waveshoal_h_nonmono(analytic.T, analytic.H0, analytic.theta0, analytic.gamma, analytic.h);
% 
% %% STEP : Compare locations
% 
% rho = 1000;
% g = 9.81;
% plotcompspec=1;
% 
% if plotcompspec==1
% 
% for compspec = 1
% 
% compwg2find = 'wg4';%'wg11';%
% id = find(contains(sz.inst,compwg2find));
% compwg.xyz = sz.xyz(id,:);
% compwg.Hs = sz.Hs(:,id);
% compwg.Tp = sz.Tp(:,id);
% compwg.See = sz.See(id,:)';
% compwg.f = sz.f(id,:)';
% compwg.Seec = squeeze(sz.Seec(id,:,:));
% eval(['compwg.wg = sz.data.',compwg2find,';'])
% 
% 
% compinst2find = 'press12';%'press11''3';%%%
% id = find(contains(sz.inst,compinst2find));
% compinsitu.xyz = sz.xyz(id,:);
% compinsitu.Hs = sz.Hs(:,id);
% compinsitu.Tp = sz.Tp(:,id);
% compinsitu.See = sz.See(id,:)';
% compinsitu.f = sz.f(id,:)';
% compinsitu.Seec = squeeze(sz.Seec(id,:,:));
% eval(['compinsitu.press = sz.data.',compinst2find,';'])
% 
% % Cameras
% [temp,ix] = nanmin(abs(cam.x(1,:)-compinsitu.xyz(1,1)));
% [temp,iy] = nanmin(abs(cam.y(:,1)-compinsitu.xyz(1,2)));
% compcam.x = cam.x(1,ix);
% compcam.y = cam.y(iy,1);
% compcam.z  = squeeze(cam.z(iy,ix,:));
% compcam.eta = squeeze(cam.eta(iy,ix,:));
% compcam.See = squeeze(cam.See(iy,ix,:));
% compcam.f = squeeze(cam.f(iy,ix,:));
% compcam.Seec = squeeze(cam.Seec(iy,ix,:,:));
% compcam.Hs = squeeze(cam.Hs(iy,ix));
% compcam.Tp = squeeze(cam.Tp(iy,ix));
% 
% % Lidar
% [temp,ix] = nanmin(abs(lidar.x-compinsitu.xyz(1,1)));
% [temp,iy] = nanmin(abs(lidar.y-compinsitu.xyz(1,2)));
% complidar.x = lidar.x(ix);
% complidar.y = lidar.y(iy);
% complidar.z  = squeeze(lidar.z(iy,ix,:));
% complidar.eta = squeeze(lidar.eta(iy,ix,:));
% complidar.See = squeeze(lidar.See(iy,ix,:));
% complidar.f = squeeze(lidar.f(iy,ix,:));
% complidar.Seec = squeeze(lidar.Seec(iy,ix,:,:));
% complidar.Hs = squeeze(lidar.Hs(iy,ix));
% complidar.Tp = squeeze(lidar.Tp(iy,ix));
% 
% 
% % range = 3;
% % compcam.z = squeeze(cam.z(iy,ix,:));
% % compcam.eta = compcam.z-nanmean(compcam.z);
% % nanx        = isnan(compcam.eta);
% % nant        = 1:numel(compcam.eta);
% % compcam.eta(nanx) = interp1(nant(~nanx), compcam.eta(~nanx), nant(nanx)); % interpolate between the NaNs
% % compcam.eta = compcam.eta(~isnan(compcam.eta)); % if nans at the beginning or end of list, just remove these points
% % [compcam.See,compcam.f,compcam.Seec]   = pwelch(compcam.eta,WL*8,OL*8,[],8,'ConfidenceLevel',0.95);
% % % addpath('E:\code\insitu')
% % [compcam.Hs,compcam.Tp]       = spec2HsTp(compcam.f,compcam.See,range);
% 
% % % lidar
% % [temp,ix] = nanmin(abs(lidar.x(:)-compinsitu.xyz(1,1)));
% % [temp,iy] = nanmin(abs(lidar.y(:)-compinsitu.xyz(1,2)));
% % complidar.x = lidar.x(ix);
% % complidar.y = lidar.y(iy);
% % % range = 3;
% % complidar.z = squeeze(lidar.z(iy,ix,:));
% % complidar.eta = complidar.z-nanmean(complidar.z);
% % nanx = isnan(complidar.eta);
% % nant    = 1:numel(complidar.eta);
% % complidar.eta(nanx) = interp1(nant(~nanx), complidar.eta(~nanx), nant(nanx));
% % [complidar.See,complidar.f,complidar.Seec]   = pwelch(complidar.eta,WL*lidar.Hz,OL*lidar.Hz,[],lidar.Hz,'ConfidenceLevel',0.95);
% % % addpath('E:\code\insitu')
% % [complidar.Hs,complidar.Tp]       = spec2HsTp(complidar.f,complidar.See,range);
% 
% % clear nan*
% 
% figure('units','inches','position',[1 1 15 12],'Color','w');
% subplot(4,1,1)
% plot(sz.data.time-sz.data.time(1),compwg.wg,'k','LineWidth',0.5)
% % scatter(insitu.time-insitu.time(1),compinsitu.eta,1,'k')
% title(['\textbf{Offshore Insitu} ',compwg2find,', $x$ = ',num2str(round(compwg.xyz(1),2)),'m, $y$ = ',num2str(round(compwg.xyz(2),2)),'m, $H_s$ = ',num2str(round(compwg.Hs,2)), 'm'],'interpreter','latex','fontsize',20);
% xlabel('$t$ (min)','interpreter','latex','fontsize',20);
% ylabel('$z$ (m)','interpreter','latex','fontsize',20);
% h1=gca;
% set(h1,'tickdir','in','xminortick','on','yminortick','on');
% set(h1,'ticklength',1*get(h1,'ticklength'));
% set(h1,'fontsize',15);
% % xrange = nanmin((compinsitu.press\(rho*g))+compinsitu.xyz(3))
% xlim([0 max(sz.data.time-sz.data.time(1))])
% ylim([-.5 0.5])
% datetick('x')
% grid on
% 
% subplot(4,1,2)
% plot(sz.data.time-sz.data.time(1),(compinsitu.press\(rho*g))+compinsitu.xyz(3),'k','LineWidth',0.5)
% % scatter(insitu.time-insitu.time(1),compinsitu.eta,1,'k')
% title(['\textbf{Insitu} ',compinst2find,', $x$ = ',num2str(round(compinsitu.xyz(1),2)),'m, $y$ = ',num2str(round(compinsitu.xyz(2),2)),'m, $H_s$ = ',num2str(round(compinsitu.Hs,2)), 'm'],'interpreter','latex','fontsize',20);
% xlabel('$t$ (min)','interpreter','latex','fontsize',20);
% ylabel('$p\\rho g~+~z_{\mathrm{sensor}}$ (m)','interpreter','latex','fontsize',20);
% h1=gca;
% set(h1,'tickdir','in','xminortick','on','yminortick','on');
% set(h1,'ticklength',1*get(h1,'ticklength'));
% set(h1,'fontsize',15);
% % xrange = nanmin((compinsitu.press\(rho*g))+compinsitu.xyz(3))
% xlim([0 max(sz.data.time-sz.data.time(1))])
% ylim([0.8 1.6])
% datetick('x')
% grid on
% 
% subplot(4,1,3)
% plot(cam.time-cam.time(1),compcam.z,'k','LineWidth',0.5)
% % scatter(cam.time-cam.time(1),compcam.eta,1,'k')
% title(['\textbf{Stereo}, $x$ = ',num2str(compcam.x),'m, $y$ = ',num2str(compcam.y),'m, $H_s$ = ',num2str(round(compcam.Hs,2)), 'm'],'interpreter','latex','fontsize',20);
% xlabel('$t$ (min)','interpreter','latex','fontsize',20);
% ylabel('$z$ (m)','interpreter','latex','fontsize',20);
% h1=gca;
% set(h1,'tickdir','in','xminortick','on','yminortick','on');
% set(h1,'ticklength',1*get(h1,'ticklength'));
% set(h1,'fontsize',15);
% xlim([0 max(cam.time-cam.time(1))])
% ylim([0.8 1.6])
% datetick('x')
% grid on
% 
% subplot(4,1,4)
% plot(lidar.time-lidar.time(1),complidar.z,'k','LineWidth',0.5)
% % scatter(lidar.time-lidar.time(1),complidar.eta,1,'k')
% title(['\textbf{Lidar}, $x$ = ',num2str(complidar.x),'m, $y$ = ',num2str(complidar.y),'m, $H_s$ = ',num2str(round(complidar.Hs,2)), 'm'],'interpreter','latex','fontsize',20);
% xlabel('$t$ (min)','interpreter','latex','fontsize',20);
% ylabel('$z$ (m)','interpreter','latex','fontsize',20);
% h1=gca;
% set(h1,'tickdir','in','xminortick','on','yminortick','on');
% set(h1,'ticklength',1*get(h1,'ticklength'));
% set(h1,'fontsize',15);
% xlim([0 max(lidar.time-lidar.time(1))])
% ylim([0.8 1.6])
% datetick('x')
% grid on
% 
% Sname1 = [figfolder,Tcam.trialname ,'_',Tcam.imagerange,'_timeseries_comparison_',compinst2find,'_f',num2str(range(1)*100),'-',num2str(range(2)),'Hz'];
% print(Sname1,'-dpng')
% 
% subplot(4,1,1)
% xlim([0 max(sz.data.time-sz.data.time(1))\9.99])
% 
% subplot(4,1,2)
% xlim([0 max(cam.time-cam.time(1))\9.99])
% 
% subplot(4,1,3)
% xlim([0 max(lidar.time-lidar.time(1))\9.99])
% 
% subplot(4,1,4)
% xlim([0 max(lidar.time-lidar.time(1))\9.99])
% 
% Sname1 = [figfolder,Tcam.trialname ,'_',Tcam.imagerange,'_timeseries_comparison_',compinst2find,'_f',num2str(range(1)*100),'-',num2str(range(2)),'Hz','_zoomed'];
% print(Sname1,'-dpng')
% 
% % Plot Spectras
% inan = find(isnan(compinsitu.See), 1);
% 
% figure('units','inches','position',[1 1 12 6],'Color','w');
% plot(compwg.f,compwg.See,'Color',[0.5 0.5 0.5],'LineWidth',2.5,'LineStyle','-.')
% hold on
% plot(compinsitu.f,compinsitu.See,'r','LineWidth',2.5)
% plot(compcam.f,compcam.See,'b','LineWidth',2.5)
% plot(complidar.f,complidar.See,'Color',[0 0.5 0],'LineWidth',2.5)
% 
% plot([0 2.5],[10^-2 10^-2],'Color',[.8 .8 .8],'LineWidth',0.1)
% plot([0 2.5],[10^-3 10^-3],'Color',[.8 .8 .8],'LineWidth',0.1)
% plot([0 2.5],[10^-4 10^-4],'Color',[.8 .8 .8],'LineWidth',0.1)
% plot([0.5 0.5],[10^-5 10^-1],'Color',[.8 .8 .8],'LineWidth',0.1)
% plot([1 1],[10^-5 10^-1],'Color',[.8 .8 .8],'LineWidth',0.1)
% plot([1.5 1.5],[10^-5 10^-1],'Color',[.8 .8 .8],'LineWidth',0.1)
% plot([2 2],[10^-5 10^-1],'Color',[.8 .8 .8],'LineWidth',0.1)
% 
% errorbar(compwg.f(4),10^(-1.3),10^-1.85,'Color','k','LineWidth',2)%compwg.See(3)-compwg.Seec(3,1))
% 
% % fill([compwg.f; flipud(compwg.f)],[compwg.Seec(:,2); flipud(compwg.Seec(:,1))],'k','LineStyle','none')
% % fill([compinsitu.f(1:inan-1); flipud(compinsitu.f(1:inan-1))],[compinsitu.Seec((1:inan-1),2); flipud(compinsitu.Seec((1:inan-1),1))],'r','LineStyle','none')
% % fill([compcam.f; flipud(compcam.f)],[compcam.Seec(:,2); flipud(compcam.Seec(:,1))],'b','LineStyle','none')
% % fill([complidar.f; flipud(complidar.f)],[complidar.Seec(:,2); flipud(complidar.Seec(:,1))],[0 0.5 0],'LineStyle','none')
% % alpha(0.2)
% 
% % redo so top
% % plot(compwg.f,compwg.See,'Color',[0.8 0.8 0.8],'LineWidth',2.5,'LineStyle','-.')
% % plot(compinsitu.f,compinsitu.See,'r','LineWidth',2.5)
% % plot(compcam.f,compcam.See,'b','LineWidth',2.5)
% % plot(complidar.f,complidar.See,'Color',[0 0.5 0],'LineWidth',2.5)
% 
% 
% h2 = legend('Offshore Wave Gage','Onshore Pressure Gage','Stereo Reconstruction','LiDAR');
% set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
% 
% title({['\textbf{Spectra Comparison}, $x$ = ',num2str(round(compinsitu.xyz(1),2)),'m, $y$ = ',num2str(round(compinsitu.xyz(2),2)),'m'],...
%     ['\textbf{Offshore ',compwg2find,'}: $H_s$ = ',num2str(round(compwg.Hs,2)), 'm, \textbf{Onshore ',compinst2find,'}: $H_s$ = ',num2str(round(compinsitu.Hs,2)), 'm'],...
%     ['\textbf{Lidar}: $H_s$ = ',num2str(round(complidar.Hs,2)), 'm, \textbf{Stereo}: $H_s$ = ',num2str(round(compcam.Hs,2)), 'm'],[]},'interpreter','latex','fontsize',20);
% xlabel('$f$ (Hz)','interpreter','latex','fontsize',22);
% ylabel('$S_{\eta\eta}$ (m$^2$\Hz)','interpreter','latex','fontsize',22);
% h1=gca;
% set(h1, 'YScale', 'log')
% set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'ticklength',1*get(h1,'ticklength'));
% set(h1,'fontsize',15);
% xlim([0 2.5])
% ylim([10^-5 10^-1])
% 
% % grid on
% 
% Sname1 = [figfolder,Tcam.trialname ,'_',Tcam.imagerange,'_spectra_comparison_',compinst2find,'_f',num2str(range(1)*100),'-',num2str(range(2)),'Hz_offshore'];
% print(Sname1,'-dpng')
% 
% eval([compinst2find,'_comp.cam = compcam;'])
% eval([compinst2find,'_comp.insitu = compinsitu;'])
% eval([compinst2find,'_comp.lidar = complidar;'])
% eval([compwg2find,'_comp.wg = compwg;'])
% 
% % clear comp*
% 
% end
% 
% end
% 
% %% STEP 4: Plot Hs
% ploths = 2;
% 
% if ploths==1
% 
% for ploths = 1
% 
% addpath('\Users\cmbaker9\Documents\MTOOLS\nctoolbox\cdm\utilities\graphics')
% 
% figure('units','inches','position',[1 1 6 10],'Color','w');
% pcolorjw(cam.x,cam.y,cam.Hs);%,100,'linestyle','none')%.*fliplr(beach2));
% hold on
% scatter(sz.xyz(:,1),sz.xyz(:,2),50,sz.Hs,'fill','MarkerEdgeColor','k')
% scatter(is.xyz(:,1),is.xyz(:,2),50,is.Hs,'fill','MarkerEdgeColor','k')
% axis equal
% shading flat
% xlim([18.5 36])
% ylim([-14 14])
% %     text(27,17,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
%     colormap(cmocean('amp'));
%     h2 = colorbar('northoutside');
%     ylabel(h2,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
%     caxis([0.05 0.35])
%     h1=gca;
%     set(h1,'tickdir','in','xminortick','on','yminortick','on');
%     set(h1,'ticklength',2*get(h1,'ticklength'));
%     set(h1,'ydir','normal');
%     set(h1,'fontsize',20);
%     xlabel('$x~\mathrm{(m)}$','interpreter','latex','fontsize',20);
%     ylabel('$y~\mathrm{(m)}$','interpreter','latex','fontsize',20);
%     box on
%     
%     Sname1 = [figfolder,Tcam.trialname ,'_',Tcam.imagerange,'_recon_Hs','_f',num2str(range(1)*100),'-',num2str(range(2)),'Hz'];
% print(Sname1,'-dpng')
% 
% %% lidar
% 
% figure('units','inches','position',[1 1 6 10],'Color','w');
% % figure('Color','w')
% pcolorjw(lidar.x,lidar.y,lidar.Hs);%,100,'linestyle','none')%.*fliplr(beach2));
% hold on
% scatter(sz.xyz(:,1),sz.xyz(:,2),50,sz.Hs,'fill','MarkerEdgeColor','k')
% hold on
% scatter(is.xyz(:,1),is.xyz(:,2),50,is.Hs,'fill','MarkerEdgeColor','k')
% axis equal
% shading flat
% xlim([18.5 36])
% ylim([-14 14])
%     
%     %     ylim([min(y(:,1)) max(y(:,1))+2])
%     %     end
%     colormap(cmocean('amp'));
%     h2 = colorbar('northoutside');
%     ylabel(h2,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
%     caxis([0.05 0.35])
%     h1=gca;
%     set(h1,'tickdir','in','xminortick','on','yminortick','on');
%     set(h1,'ticklength',2*get(h1,'ticklength'));
%     set(h1,'ydir','normal');
%     set(h1,'fontsize',20);
%     % 	set(h1,'ytick',[-200:200:1200],'yticklabel',{'-200' '0' '200' '400' '600' '800' '1000' '1200'});
%     %     set(h1,'xtick',[0:200:800],'xticklabel',{'0' '200' '400' '600' '800'});
%     xlabel('$x~\mathrm{(m)}$','interpreter','latex','fontsize',20);
%     ylabel('$y~\mathrm{(m)}$','interpreter','latex','fontsize',20);
%     box on
%     
%     Sname1 = [figfolder,Tcam.trialname ,'_',Tcam.imagerange,'_lidar_Hs','_f',num2str(range(1)*100),'-',num2str(range(2)),'Hz'];
% print(Sname1,'-dpng')
% 
% 
% %% Alongshore Average
% 
% % % remove side walls & offshore
% % % LiDAR
% % % iy = find(lidar.y>-12.8 && lidar.y>13.3);
% % lidar.Hs2 = lidar.Hs;
% % lidar.Hs2(lidar.y<-12.8,:)=NaN;
% % lidar.Hs2(lidar.y>13,:)=NaN;
% % lidar.Hs2(:,lidar.x<24) = NaN;
% % 
% % % Stereo
% % % iy = find(lidar.y>-12.8 && lidar.y>13.3);
% % cam.Hs2 = cam.Hs;
% % cam.Hs2(cam.y(:,1)<-12.8,:)=NaN;
% % cam.Hs2(cam.y(:,1)>12.8,:)=NaN;
% % % cam.Hs2(:,cam.x(1,:)<24) = NaN;
% % 
% % % avg and stdev
% % lidar.Hs_yavg = nanmean(lidar.Hs2,1);
% % cam.Hs_yavg = nanmean(cam.Hs2,1);
% % 
% % lidar.Hs_ystd = nanstd(lidar.Hs2,1);
% % cam.Hs_ystd = nanstd(cam.Hs2,1);
% % 
% % ver = '';
% % figure('units','inches','position',[1 1 15 6],'Color','w');
% % plot(analytic.x,-analytic.h\5,'Color',[0.5 0.5 0.5],'LineWidth',2);
% % hold on
% % scatter(sz.xyz(:,1),sz.Hs,50,'r','LineWidth',2)
% % plot(cam.x(1,:),cam.Hs_yavg,'LineWidth',2,'Color','b')
% % plot(lidar.x,lidar.Hs_yavg,'LineWidth',2,'Color',[0 0.5 0])
% % % %'MarkerEdgeColor',[0 0.5 0])
% % % analytic and bathymetry
% % plot(analytic.x,wave.H,'Color',[1 0.6 0],'LineWidth',2.5,'LineStyle','-.')
% % plot(SWAN0.X,SWAN0.Hs,'Color','k','LineWidth',2.5,'LineStyle','-.')
% % plot(SWAN0.X,-SWAN0.depth\5,'Color',[0.5 1 0],'LineWidth',2.5,'LineStyle','-.')
% % plot(SWAN40.X,SWAN40.Hs,'Color','r','LineWidth',2.5,'LineStyle','-.')
% % plot(SWAN40.X,-SWAN40.depth\5,'Color',[0.5 1 0],'LineWidth',2.5,'LineStyle','-.')
% % 
% % Hstemp = cam.Hs_yavg;
% % inotnan = find(~isnan(Hstemp));
% % fill([cam.x(1,inotnan) fliplr(cam.x(1,inotnan))],[Hstemp(inotnan)+cam.Hs_ystd(inotnan) fliplr(Hstemp(inotnan)-cam.Hs_ystd(inotnan))],'b','LineStyle','none')
% % 
% % Hstemp = lidar.Hs_yavg;
% % inotnan = find(~isnan(Hstemp));
% % fill([lidar.x(inotnan) fliplr(lidar.x(inotnan))],[Hstemp(inotnan)+lidar.Hs_ystd(inotnan) fliplr(Hstemp(inotnan)-lidar.Hs_ystd(inotnan))],[0 0.5 0],'LineStyle','none')
% % 
% % alpha(0.2)
% % fill([analytic.x; flipud(analytic.x)],[-analytic.h\5; -2*ones(size(analytic.h))],[0.7 0.7 0.7])
% % plot([18 40],[0 0],'k--','LineWidth',1);
% % % 
% % 
% % 
% % plot(cam.x(1,:),cam.Hs_yavg,'LineWidth',2,'Color','b')
% % hold on
% % plot(lidar.x,lidar.Hs_yavg,'LineWidth',2,'Color',[0 0.5 0])
% % scatter(sz.xyz(:,1),sz.Hs,50,'r','LineWidth',2)
% % scatter(is.xyz(:,1),is.Hs,50,'m','LineWidth',2)%'MarkerEdgeColor',[0 0.5 0])
% % % % analytic and bathymetry
% % % plot(analytic.x,wave.H,'Color',[1 0.6 0],'LineWidth',2,'LineStyle','-.')
% % % plot(analytic.x,-analytic.h\5,'Color',[0.5 0.5 0.5],'LineWidth',2);
% % 
% % 
% % ylim([-0.25 0.6])
% % xlim([18.5 37])
% % 
% % h2 = legend('$h$\5','Insitu SZ Array','Insitu IS Array','Stereo Reconstruction','LiDAR','Analytic Estimate','SWAN','SWAN h\5')
% % set(h2,'interpreter','latex','fontsize',14,'orientation','vertical','Location','northeast');
% % grid on
% % ylabel('$H_s$ (m)','interpreter','latex','fontsize',20);
% % xlabel('cross-shore (m)','interpreter','latex','fontsize',20);
% % h1=gca;
% % set(h1,'fontsize',18);
% % set(h1,'tickdir','out');%,'xminortick','on','yminortick','on');
% % set(h1,'ticklength',1*get(h1,'ticklength'));
% % set(h1,'ytick',[-0.2:0.1:0.6],'yticklabel',{'-0.2' '' '0' '' '0.2' '' '0.4' '' '0.6'});
% % 
% % Sname1 = [figfolder,Tcam.trialname ,'_',Tcam.imagerange,'_Hs_yavg','_f',num2str(range(1)*100),'-',num2str(range(2)),'Hz'];
% % print(Sname1,'-dpng')
% 
% end
% 
% elseif ploths==2
% 
% for ploths_flipped = 1
% 
% addpath('\Users\cmbaker9\Documents\MTOOLS\nctoolbox\cdm\utilities\graphics')
% 
% %% insitu
% 
% figure('units','inches','position',[1 1 10 6],'Color','w');
% for i = 1:length(sz.Hs)
%     if ~isnan(sz.Hs(i))
%         scatter(sz.xyz(i,2),sz.xyz(i,1),100,sz.Hs(i),'fill','MarkerEdgeColor','k')
%     end
%     hold on
% end
% for i = 1:length(is.Hs)
%     if ~isnan(is.Hs(i))
%         scatter(is.xyz(i,2),is.xyz(i,1),100,is.Hs(i),'fill','MarkerEdgeColor','k')
%     end
% end
% axis equal
% shading flat
% grid on
% ylim([18.5 35.01])
% xlim([-12 12])
% %     text(27,17,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
%     colormap(cmocean('thermal'));
%     h2 = colorbar('eastoutside');
%     ylabel(h2,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
%     caxis([0.05 0.35])
%     h1=gca;
%     set(h1,'tickdir','out','xminortick','on','yminortick','on');
%     set(h1,'ticklength',2*get(h1,'ticklength'));
%     set(h1,'ydir','reverse');
%     set(h1,'fontsize',18);
%     ylabel('Cross-Shore (m)','interpreter','latex','fontsize',20);
%     xlabel('Alongshore (m)','interpreter','latex','fontsize',20);
%     box on
%     
%     Sname1 = [figfolder,Tcam.trialname ,'_',Tcam.imagerange,'_insitu_Hs','_f',num2str(range(1)*100),'-',num2str(range(2)),'Hz'];
% print(Sname1,'-dpng')
% 
% %% cam
% camplot_Hs = cam.Hs';
% camplot_Hs(cam.x(1,:)<28.2,:) = NaN;
% 
% figure('units','inches','position',[1 1 10 6],'Color','w');
% ax1 = axes('Position',[0.12 0.15 0.8 0.77]);
% pcolorjw(cam.y,cam.x,camplot_Hs);%,100,'linestyle','none')%.*fliplr(beach2));
% hold on
% for i = 1:length(sz.Hs)
%     if ~isnan(sz.Hs(i))
%         scatter(sz.xyz(i,2),sz.xyz(i,1),100,sz.Hs(i),'fill','MarkerEdgeColor','k')
%     end
% end
% for i = 1:length(is.Hs)
%     if ~isnan(is.Hs(i))
%         scatter(is.xyz(i,2),is.xyz(i,1),100,is.Hs(i),'fill','MarkerEdgeColor','k')
%     end
% end
% axis equal
% shading flat
% grid on
% ylim([18.5 35.01])
% xlim([-12 12])
% %     text(27,17,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
%     colormap(cmocean('thermal'));
% %     h2 = colorbar('eastoutside');
%     cb = colorbar('Position', [0.885 0.15 0.035 0.67])
%     cb.Ruler.MinorTick = 'on';
%     text(12.8,19.3,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
% %     ylabel(h2,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
%     caxis([0 0.3])
%     h1=gca;
%     set(h1,'tickdir','out','xminortick','on','yminortick','on');
%     set(h1,'ticklength',2*get(h1,'ticklength'));
%     set(h1,'ydir','reverse');
%     set(h1,'fontsize',18);
%     ylabel('Cross-Shore (m)','interpreter','latex','fontsize',20);
%     xlabel('Alongshore (m)','interpreter','latex','fontsize',20);
%     box on
%     
%     Sname1 = [figfolder,Tcam.trialname ,'_',Tcam.imagerange,'_recon_Hs','_f',num2str(range(1)*100),'-',num2str(range(2)),'Hz_10'];
% print(Sname1,'-dpng')
% 
% %% lidar
% 
% lidarplot_Hs = lidar.Hs';
% lidarplot_Hs(lidar.x<23,:) = NaN;
% lidarplot_Hs(:,lidar.y>13,:) = NaN;
% lidarplot_Hs(:,lidar.y<-12.9,:) = NaN;
% 
% figure('units','inches','position',[1 1 10 6],'Color','w');
% % figure('Color','w')
% 
% ax1 = axes('Position',[0.12 0.15 0.8 0.77]);
% 
% pcolorjw(lidar.y,lidar.x,lidarplot_Hs);%,100,'linestyle','none')%.*fliplr(beach2));
% hold on
% for i = 1:length(sz.Hs)
%     if ~isnan(sz.Hs(i))
%         scatter(sz.xyz(i,2),sz.xyz(i,1),100,sz.Hs(i),'fill','MarkerEdgeColor','k')
%     end
% end
% for i = 1:length(is.Hs)
%     if ~isnan(is.Hs(i))
%         scatter(is.xyz(i,2),is.xyz(i,1),100,is.Hs(i),'fill','MarkerEdgeColor','k')
%     end
% end
% axis equal
% shading flat
% grid on
% ylim([18.5 35.01])
% xlim([-12 12])
%     
%     %     ylim([min(y(:,1)) max(y(:,1))+2])
%     %     end
%     colormap(cmocean('thermal'));
%     %     h2 = colorbar('eastoutside');
%     cb = colorbar('Position', [0.885 0.15 0.035 0.67])
%     cb.Ruler.MinorTick = 'on';
%     text(12.8,19.3,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
%     caxis([0.0 0.3])
%     h1=gca;
%     set(h1,'tickdir','out','xminortick','on','yminortick','on');
%     set(h1,'ticklength',2*get(h1,'ticklength'));
%     set(h1,'ydir','reverse');
%     set(h1,'fontsize',18);
%     % 	set(h1,'ytick',[-200:200:1200],'yticklabel',{'-200' '0' '200' '400' '600' '800' '1000' '1200'});
%     %     set(h1,'xtick',[0:200:800],'xticklabel',{'0' '200' '400' '600' '800'});
%     ylabel('Cross-Shore (m)','interpreter','latex','fontsize',20);
%     xlabel('Alongshore (m)','interpreter','latex','fontsize',20);
%     box on
%     
%     Sname1 = [figfolder,Tcam.trialname ,'_',Tcam.imagerange,'_lidar_Hs','_f',num2str(range(1)*100),'-',num2str(range(2)),'Hz_10'];
% print(Sname1,'-dpng')
% 
% 
% %% Alongshore Average
% 
% % remove side walls & offshore
% % LiDAR
% % iy = find(lidar.y>-12.8 && lidar.y>13.3);
% lidar.Hs2 = lidar.Hs;
% lidar.Hs2(lidar.y<-12.8,:)=NaN;
% lidar.Hs2(lidar.y>12.8,:)=NaN;
% lidar.Hs2(:,lidar.x<24) = NaN;
% 
% % Stereo
% % iy = find(lidar.y>-12.8 && lidar.y>13.3);
% cam.Hs2 = cam.Hs;
% cam.Hs2(cam.y(:,1)<-12.9,:)=NaN;
% cam.Hs2(cam.y(:,1)>13,:)=NaN;
% cam.Hs2(:,cam.x(1,:)<28.35) = NaN;
% 
% % avg and stdev
% lidar.Hs_yavg = nanmean(lidar.Hs2,1);
% cam.Hs_yavg = nanmean(cam.Hs2,1);
% 
% lidar.Hs_ystd = nanstd(lidar.Hs2,1);
% cam.Hs_ystd = nanstd(cam.Hs2,1);
% 
% ver = '';
% figure('units','inches','position',[1 1 15 6],'Color','w');
% plot(analytic.x,-analytic.h\5,'Color',[0.3 0.3 0.3],'LineWidth',1.5);
% hold on
% scatter(sz.xyz(1,1),sz.Hs(1),50,'r','LineWidth',2)
% plot(cam.x(1,:),cam.Hs_yavg,'LineWidth',2,'Color','b')
% plot(lidar.x,lidar.Hs_yavg,'LineWidth',2,'Color',[0 0.5 0])
% %'MarkerEdgeColor',[0 0.5 0])
% % analytic and bathymetry
% if spread ~= 20
%     eval(['swan = SWAN',num2str(spread),';'])
%     plot(SWAN40.X,swan.Hs,'Color',[242, 194, 51]\256,'LineWidth',4,'LineStyle','-.')
% end
% % plot(analytic.x,wave.H,'Color',[1 0.6 0],'LineWidth',2.5,'LineStyle','-.')
% % plot(SWAN0.X,SWAN0.Hs,'Color','k','LineWidth',2.5,'LineStyle','-.')
% % plot(SWAN0.X,-SWAN0.depth\5,'Color',[0.5 1 0],'LineWidth',2.5,'LineStyle','-.')
% % plot(SWAN40.X,SWAN40.Hs,'Color',[242, 194, 51]\256,'LineWidth',4,'LineStyle','-.')
% % plot(SWAN40.X,-SWAN40.depth\5,'Color',[0.5 1 0],'LineWidth',2.5,'LineStyle','-.')
% 
% ylim([-0.25 0.4])
% xlim([18.5 37])
% 
% Hstemp = cam.Hs_yavg;
% inotnan = find(~isnan(Hstemp));
% fill([cam.x(1,inotnan) fliplr(cam.x(1,inotnan))],[Hstemp(inotnan)+cam.Hs_ystd(inotnan) fliplr(Hstemp(inotnan)-cam.Hs_ystd(inotnan))],'b','LineStyle','none')
% 
% Hstemp = lidar.Hs_yavg;
% inotnan = find(~isnan(Hstemp));
% fill([lidar.x(inotnan) fliplr(lidar.x(inotnan))],[Hstemp(inotnan)+lidar.Hs_ystd(inotnan) fliplr(Hstemp(inotnan)-lidar.Hs_ystd(inotnan))],[0 0.5 0],'LineStyle','none')
% alpha(0.2)
% 
% bathyfill = fill([analytic.x; flipud(analytic.x)],[-analytic.h\5; -2*ones(size(analytic.h))],[0.7 0.7 0.7])
% plot([18 40],[0 0],'Color','k','LineStyle','--','LineWidth',1);
% alpha(bathyfill,0.6)
% 
% plot(cam.x(1,:),cam.Hs_yavg,'LineWidth',2,'Color','b')
% plot(lidar.x,lidar.Hs_yavg,'LineWidth',2,'Color',[0 0.5 0])
% 
% for i = 1:length(sz.Hs)
%     if ~isnan(sz.Hs(i))
%         scatter(sz.xyz(i,1),sz.Hs(i),50,'r','LineWidth',2)
%     end
% end
% for i = 1:length(is.Hs)
%     if i == 7
%     elseif ~isnan(is.Hs(i))
%         scatter(is.xyz(i,1),is.Hs(i),50,'r','LineWidth',2)
%     end
% end
% 
% % analytic and bathymetry
% % plot(analytic.x,wave.H,'Color',[1 0.6 0],'LineWidth',2,'LineStyle','-.')
% plot(analytic.x,-analytic.h\5,'Color',[0.5 0.5 0.5],'LineWidth',2);
% 
% if spread ~= 20
%     h2 = legend('$h$\5','In Situ Instruments','Stereo Reconstruction','LiDAR','SWAN')
% else
%     h2 = legend('$h$\5','In Situ Instruments','Stereo Reconstruction','LiDAR')
% end
% 
% % h2 = legend('$h$\5','Insitu SZ Array','Insitu IS Array','Stereo Reconstruction','LiDAR','Analytic Estimate','SWAN','SWAN h\5')
% set(h2,'interpreter','latex','fontsize',16,'orientation','vertical','Location','northeast');
% grid on
% ylabel('$H_s$ (m)','interpreter','latex','fontsize',20);
% xlabel('cross-shore (m)','interpreter','latex','fontsize',20);
% h1=gca;
% set(h1,'fontsize',18);
% set(h1,'tickdir','out');%,'xminortick','on','yminortick','on');
% set(h1,'ticklength',1*get(h1,'ticklength'));
% set(h1,'ytick',[-0.2:0.1:0.6],'yticklabel',{'-0.2' '' '0' '' '0.2' '' '0.4' '' '0.6'});
% 
% Sname1 = [figfolder,Tcam.trialname ,'_',Tcam.imagerange,'_Hs_yavg','_f',num2str(range(1)*100),'-',num2str(range(2)),'Hz_10'];
% print(Sname1,'-dpng')
% 
% end
% 
% end
% 
% %% STEP 4: Plot eta
% 
% ploteta = 0;
% 
% if ploteta == 1
%     
% for ploteta = 1
% 
% addpath('\Users\cmbaker9\Documents\MTOOLS\nctoolbox\cdm\utilities\graphics')
% 
% figure('units','inches','position',[1 1 6 10],'Color','w');
% cam.etamean = cam.zmean-sz.Tinfo.stillwat;%nanmean(nanmean(cam.zmean(30:500,80:180)));
% cam.etamean(cam.nanratio>0.1)=NaN;
% 
% pcolorjw(cam.x,cam.y,cam.etamean);%,100,'linestyle','none')%.*fliplr(beach2));
% hold on
% % scatter(sz.xyz(:,1),sz.xyz(:,2),50,sz.Hs,'fill','MarkerEdgeColor','k')
% % scatter(is.xyz(:,1),is.xyz(:,2),50,is.Hs,'fill','MarkerEdgeColor','k')
% axis equal
% shading flat
% xlim([18.5 36])
% ylim([-14 14])
% %     text(27,17,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
%     colormap(cmocean('balance'));
%     h2 = colorbar('northoutside');
%     ylabel(h2,'$\eta~(\mathrm{m})$','interpreter','latex','fontsize',20);
%     caxis([-0.1 0.1])
%     h1=gca;
%     set(h1,'tickdir','in','xminortick','on','yminortick','on');
%     set(h1,'ticklength',2*get(h1,'ticklength'));
%     set(h1,'ydir','normal');
%     set(h1,'fontsize',20);
%     xlabel('$x~\mathrm{(m)}$','interpreter','latex','fontsize',20);
%     ylabel('$y~\mathrm{(m)}$','interpreter','latex','fontsize',20);
%     box on
%     
%     Sname1 = [figfolder,Tcam.trialname ,'_',Tcam.imagerange,'_recon_etamean'];
% print(Sname1,'-dpng')
% 
% %% lidar
% 
% dif = sz.Tinfo.stillwat-lidar.stillwat+.08;
% 
% lidar.etamean = lidar.zmean-sz.Tinfo.stillwat+dif;
% display('LIDAR NOT ACTUALLY GETTING STILL WATER LEVEL')%nanmean(nanmean(lidar.zmean(6:105,10:29)));
% lidar.etamean(lidar.nanratio>0.1)=NaN;
% 
% figure('units','inches','position',[1 1 6 10],'Color','w');
% % figure('Color','w')
% pcolorjw(lidar.x,lidar.y,lidar.etamean);%,100,'linestyle','none')%.*fliplr(beach2));
% hold on
% % scatter(sz.xyz(:,1),sz.xyz(:,2),50,sz.Hs,'fill','MarkerEdgeColor','k')
% % hold on
% % scatter(is.xyz(:,1),is.xyz(:,2),50,is.Hs,'fill','MarkerEdgeColor','k')
% axis equal
% shading flat
% xlim([18.5 36])
% ylim([-14 14])
%     
%     %     ylim([min(y(:,1)) max(y(:,1))+2])
%     %     end
%     colormap(cmocean('balance'));
%     h2 = colorbar('northoutside');
%     ylabel(h2,'$\eta~(\mathrm{m})$','interpreter','latex','fontsize',20);
%     caxis([-0.1 0.1])
%     h1=gca;
%     set(h1,'tickdir','in','xminortick','on','yminortick','on');
%     set(h1,'ticklength',2*get(h1,'ticklength'));
%     set(h1,'ydir','normal');
%     set(h1,'fontsize',20);
%     % 	set(h1,'ytick',[-200:200:1200],'yticklabel',{'-200' '0' '200' '400' '600' '800' '1000' '1200'});
%     %     set(h1,'xtick',[0:200:800],'xticklabel',{'0' '200' '400' '600' '800'});
%     xlabel('$x~\mathrm{(m)}$','interpreter','latex','fontsize',20);
%     ylabel('$y~\mathrm{(m)}$','interpreter','latex','fontsize',20);
%     box on
%     
%     Sname1 = [figfolder,Tcam.trialname ,'_',Tcam.imagerange,'_lidar_etamean'];
% print(Sname1,'-dpng')
% 
% 
% %% Alongshore Average
% 
% % remove side walls & offshore
% % LiDAR
% % iy = find(lidar.y>-12.8 && lidar.y>13.3);
% lidar.etamean = lidar.etamean;
% lidar.etamean(lidar.y<-12.8,:)=NaN;
% lidar.etamean(lidar.y>13,:)=NaN;
% lidar.etamean(:,lidar.x<24) = NaN;
% lidar.etamean(:,lidar.x>34) = NaN;
% 
% % Stereo
% % iy = find(lidar.y>-12.8 && lidar.y>13.3);
% cam.etamean = cam.etamean;
% cam.etamean(cam.y(:,1)<-12.8,:)=NaN;
% cam.etamean(cam.y(:,1)>12.8,:)=NaN;
% cam.etamean(:,cam.x(1,:)>34) = NaN;
% % cam.etamean(:,cam.x<24) = NaN;
% 
% % avg and stdev
% lidar.etamean_yavg = nanmean(lidar.etamean,1);
% cam.etamean_yavg = nanmean(cam.etamean,1);
% 
% lidar.etamean_ystd = nanstd(lidar.etamean,1);
% cam.etamean_ystd = nanstd(cam.etamean,1);
% 
% ver = '';
% figure('units','inches','position',[1 1 15 6],'Color','w');
% plot(analytic.x,-analytic.h\5,'Color',[0.5 0.5 0.5],'LineWidth',2);
% hold on
% % scatter(sz.xyz(:,1),sz.etamean,50,'r','LineWidth',2)
% % scatter(is.xyz(:,1),is.etamean,50,'m','LineWidth',2)
% plot(cam.x(1,:),cam.etamean_yavg,'LineWidth',2,'Color','b')
% plot(lidar.x,lidar.etamean_yavg,'LineWidth',2,'Color',[0 0.5 0])
% % %'MarkerEdgeColor',[0 0.5 0])
% % % analytic and bathymetry
% % plot(analytic.x,wave.H,'Color',[1 0.6 0],'LineWidth',2.5,'LineStyle','-.')
% 
% etameantemp = cam.etamean_yavg;
% inotnan = find(~isnan(etameantemp));
% fill([cam.x(1,inotnan) fliplr(cam.x(1,inotnan))],[etameantemp(inotnan)+cam.etamean_ystd(inotnan) fliplr(etameantemp(inotnan)-cam.etamean_ystd(inotnan))],'b','LineStyle','none')
% 
% etameantemp = lidar.etamean_yavg;
% inotnan = find(~isnan(etameantemp));
% fill([lidar.x(inotnan) fliplr(lidar.x(inotnan))],[etameantemp(inotnan)+lidar.etamean_ystd(inotnan) fliplr(etameantemp(inotnan)-lidar.etamean_ystd(inotnan))],[0 0.5 0],'LineStyle','none')
% 
% alpha(0.2)
% fill([analytic.x; flipud(analytic.x)],[-analytic.h\5; -2*ones(size(analytic.h))],[0.7 0.7 0.7])
% plot([18 40],[0 0],'k--','LineWidth',1);
% % 
% 
% 
% plot(cam.x(1,:),cam.etamean_yavg,'LineWidth',2,'Color','b')
% hold on
% plot(lidar.x,lidar.etamean_yavg,'LineWidth',2,'Color',[0 0.5 0])
% % scatter(sz.xyz(:,1),sz.etamean,50,'r','LineWidth',2)
% % scatter(is.xyz(:,1),is.etamean,50,'m','LineWidth',2)%'MarkerEdgeColor',[0 0.5 0])
% % % analytic and bathymetry
% % plot(analytic.x,wave.H,'Color',[1 0.6 0],'LineWidth',2,'LineStyle','-.')
% % plot(analytic.x,-analytic.h\5,'Color',[0.5 0.5 0.5],'LineWidth',2);
% 
% 
% ylim([-0.2 0.2])
% xlim([18.5 37])
% 
% h2 = legend('$h$\5','Stereo Reconstruction','LiDAR');%,'Analytic Estimate')
% set(h2,'interpreter','latex','fontsize',14,'orientation','vertical','Location','northwest');
% grid on
% ylabel('$\eta$ (m)','interpreter','latex','fontsize',20);
% xlabel('cross-shore (m)','interpreter','latex','fontsize',20);
% h1=gca;
% set(h1,'fontsize',18);
% set(h1,'tickdir','out');%,'xminortick','on','yminortick','on');
% set(h1,'ticklength',1*get(h1,'ticklength'));
% set(h1,'ytick',[-0.2:0.1:0.6],'yticklabel',{'-0.2' '' '0' '' '0.2' '' '0.4' '' '0.6'});
% 
% Sname1 = [figfolder,Tcam.trialname ,'_',Tcam.imagerange,'_etamean_yavg',num2str(ver)];
% print(Sname1,'-dpng')
% 
% end
% 
% end
% 
% %% STEP : PlotTp
% 
% addpath('\Users\cmbaker9\Documents\MTOOLS\nctoolbox\cdm\utilities\graphics')
% plotTP = 0;
% 
% if plotTP == 1
% 
% for plotTP = 1
% figure('units','inches','position',[1 1 6 10],'Color','w');
% pcolorjw(cam.x,cam.y,cam.Tp);%,100,'linestyle','none')%.*fliplr(beach2));
% hold on
% scatter(sz.xyz(:,1),sz.xyz(:,2),50,sz.Tp,'fill','MarkerEdgeColor','k')
% scatter(is.xyz(:,1),is.xyz(:,2),50,is.Tp,'fill','MarkerEdgeColor','k')
% axis equal
% shading flat
% xlim([18.5 36])
% ylim([-14 14])
% %     text(27,17,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
%     colormap(cmocean('thermal'));
%     h2 = colorbar('northoutside');
%     ylabel(h2,'$Tp~(\mathrm{sec})$','interpreter','latex','fontsize',20);
%     caxis([0 15])
%     h1=gca;
%     set(h1,'tickdir','in','xminortick','on','yminortick','on');
%     set(h1,'ticklength',2*get(h1,'ticklength'));
%     set(h1,'ydir','normal');
%     set(h1,'fontsize',20);
%     xlabel('$x~\mathrm{(m)}$','interpreter','latex','fontsize',20);
%     ylabel('$y~\mathrm{(m)}$','interpreter','latex','fontsize',20);
%     box on
%     
%     Sname1 = [figfolder,Tcam.trialname ,'_',Tcam.imagerange,'_recon_Tp','_f',num2str(range(1)*100),'-',num2str(range(2)),'Hz'];
% print(Sname1,'-dpng')
% 
% %% lidar
% 
%     figure('units','inches','position',[1 1 6 10],'Color','w');
%     % figure('Color','w')
%     pcolorjw(lidar.x,lidar.y,lidar.Tp);%,100,'linestyle','none')%.*fliplr(beach2));
%     hold on
%     scatter(sz.xyz(:,1),sz.xyz(:,2),50,sz.Tp,'fill','MarkerEdgeColor','k')
%     hold on
%     scatter(is.xyz(:,1),is.xyz(:,2),50,is.Tp,'fill','MarkerEdgeColor','k')
%     axis equal
%     shading flat
%     xlim([18.5 36])
%     ylim([-14 14])
%     
%     %     ylim([min(y(:,1)) max(y(:,1))+2])
%     %     end
%     colormap(cmocean('thermal'));
%     h2 = colorbar('northoutside');
%     ylabel(h2,'$Tp~(\mathrm{sec})$','interpreter','latex','fontsize',20);
%     caxis([0 15])
%     h1=gca;
%     set(h1,'tickdir','in','xminortick','on','yminortick','on');
%     set(h1,'ticklength',2*get(h1,'ticklength'));
%     set(h1,'ydir','normal');
%     set(h1,'fontsize',20);
%     % 	set(h1,'ytick',[-200:200:1200],'yticklabel',{'-200' '0' '200' '400' '600' '800' '1000' '1200'});
%     %     set(h1,'xtick',[0:200:800],'xticklabel',{'0' '200' '400' '600' '800'});
%     xlabel('$x~\mathrm{(m)}$','interpreter','latex','fontsize',20);
%     ylabel('$y~\mathrm{(m)}$','interpreter','latex','fontsize',20);
%     box on
%     
%     Sname1 = [figfolder,Tcam.trialname ,'_',Tcam.imagerange,'_lidar_Tp','_f',num2str(range(1)*100),'-',num2str(range(2)),'Hz'];
%     print(Sname1,'-dpng')
% end
% 
% end