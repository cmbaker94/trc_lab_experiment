% Plot timeresiers

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('/Users/cmbaker9/Documents/MTOOLS'))
addpath(genpath('/Users/cmbaker9/Documents/Research/Lab_Experiments/codes/insitu'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trial info
Tinfo.Hs = 0.30;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
Tinfo.spread = 40;

if Tinfo.spread == 0
    Tinfo.cam.tstart         = '09-01-2018-2214UTC';                   % time starting collection based on spreadsheet
    Tinfo.cam.tdate          = '09-01-2018-2155UTC';      % trial date and time - format ex: 09-01-2018-2155UTC
elseif Tinfo.spread == 20
    Tinfo.cam.tstart  = '08-30-2018-2226UTC';                   % time starting collection based on spreadsheet
    Tinfo.cam.tdate   = '08-30-2018-2216UTC';      % trial date and time - format ex: 09-01-2018-2155UTC 
elseif Tinfo.spread == 40
    Tinfo.cam.tstart  = '08-30-2018-2129UTC';                   % time starting collection based on spreadsheet
    Tinfo.cam.tdate   = '08-30-2018-2119UTC';      % trial date and time - format ex: 09-01-2018-2155UTC 
end


%% STEP 1: Create paths, files and naming

% general path and names
datapath    = '/Users/cmbaker9/Documents/Research/Lab_Experiments/';

% Stereo Reconstructions
% Stereo Reconstructions
Tinfo.cam.camerasys   = 'TRM';                     % camera setup - e.g., TRM (offshore) or TRC (onshore)
Tinfo.cam.scene       = '1';                       % scene number of trial - typ 1
Tinfo.cam.imagestart  = 7200;                      % images number of first frame on file 
Tinfo.cam.numframes   = 4800;                      % number of frames processed
Tinfo.cam.dx          = 0.05;
Tinfo.cam.dy          = 0.05;
Tinfo.cam.regx        = [25:Tinfo.cam.dx:37];
Tinfo.cam.regy        = [-13:Tinfo.cam.dy:13];
Tinfo.cam.Hz          = 8;

Tinfo.cam.trialname   = [Tinfo.cam.camerasys,'-',Tinfo.cam.tdate];
Tinfo.cam.imagerange  = [num2str(Tinfo.cam.imagestart,'%05.f'),'-',num2str(Tinfo.cam.imagestart+(Tinfo.cam.numframes-1),'%05.f')];
Tinfo.cam.trimname    = ['frames_',Tinfo.cam.imagerange,'/'];  
Tinfo.cam.datafolder  = [datapath,'data/processed/DEM/',Tinfo.cam.trialname,'_Scene',Tinfo.cam.scene,'/',Tinfo.cam.trimname];

%% STEP 2: Create figure folders

% figure folder
fssubfolder = datestr(date,'yy-mm-dd');
figfolder   = [datapath,'figures/meas_comp/',Tinfo.cam.trialname,'/',Tinfo.cam.trimname,fssubfolder,'/'];

% make figure folders
eval(['!mkdir ',datapath,'figures/meas_comp/',Tinfo.cam.trialname]);
eval(['!mkdir ',datapath,'figures/meas_comp/',Tinfo.cam.trialname,'/',Tinfo.cam.trimname]);
eval(['!mkdir ',figfolder])

%% Establish timeseries range:

starttemp       = datenum(Tinfo.cam.tstart(1:end-3),'mm-dd-yyyy-HHMM')+datenum(0,0,0,0,0,Tinfo.cam.imagestart/Tinfo.cam.Hz);
endtemp         = datenum(Tinfo.cam.tstart(1:end-3),'mm-dd-yyyy-HHMM')+datenum(0,0,0,0,0,(Tinfo.cam.imagestart+Tinfo.cam.numframes-1)/Tinfo.cam.Hz);
camera.time     = starttemp:datenum(0,0,0,0,0,1/Tinfo.cam.Hz):endtemp;
clear *temp

%% Data storage location

Tinfo.comp = ['Hs',num2str(Tinfo.Hs*100),'_Tp',num2str(Tinfo.Tp),'_tide',num2str(Tinfo.tide*100),'_spread',num2str(Tinfo.spread)];
datarange = [datestr(camera.time(1)); datestr(camera.time(end))];
datarange = ['time_', datarange(1,13:14), datarange(1,16:17), '-' , datarange(2,13:14), datarange(2,16:17)];
procpath = [datapath,'data/processed/conditions/',Tinfo.comp];

datafolder = [procpath,'/',Tinfo.cam.tstart,'/',datarange,'/'];

sname = 'szarray_timeseries_allinst';

%%

load([datafolder,sname,'.mat'])

%%
rho = 1000;
g = 9.81;
plotcompspec=1;

compwg2find = 'wg4';%'wg11';%
idwg = find(contains(wg.names,compwg2find));

compinst2find = 'press3';%'press11''3';%%%
idpress = find(contains(press.name,compinst2find));




figure('units','inches','position',[1 1 15 12],'Color','w');
subplot(4,1,1)
plot(wg.time,wg.z(:,idwg),'k','LineWidth',0.5)
% scatter(insitu.time-insitu.time(1),compinsitu.eta,1,'k')
title(['\textbf{Offshore Insitu: } ',char(wg.names(idwg)),', $x$ = ',num2str(round(wg.xyz(1,idwg))),' m, $y$ = ',num2str(round(wg.xyz(2,idwg))),'m'],'interpreter','latex','fontsize',20);
% xlabel('$t$ (min)','interpreter','latex','fontsize',20);
ylabel('$z$ (m)','interpreter','latex','fontsize',20);
h1=gca;
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1*get(h1,'ticklength'));
set(h1,'fontsize',15);
% xrange = nanmin((compinsitu.press/(rho*g))+compinsitu.xyz(3))
xlim([min(wg.time) max(wg.time)])
ylim([-.3 0.3])
datetick('x')
grid on

subplot(4,1,2)
tempress = (press.press(:,idpress)/(rho*g))+0.05;
plot(press.time,tempress-nanmean(tempress),'k','LineWidth',0.5)
% scatter(insitu.time-insitu.time(1),compinsitu.eta,1,'k')
title(['\textbf{Pressure Gage}: ',char(press.name(idpress)),', $x$ = ',num2str(round(press.xyz(1,idpress))),' m, $y$ = ',num2str(round(press.xyz(2,idpress))),' m'],'interpreter','latex','fontsize',20);
% xlabel('$t$ (min)','interpreter','latex','fontsize',20);
ylabel('$p/\rho g~+~h_{\mathrm{sensor}}$ (m)','interpreter','latex','fontsize',20);
h1=gca;
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1*get(h1,'ticklength'));
set(h1,'fontsize',15);
% xrange = nanmin((compinsitu.press/(rho*g))+compinsitu.xyz(3))
xlim([min(press.time) max(press.time)])
ylim([-.2 0.2])
datetick('x')
grid on

subplot(4,1,3)
plot(camera.time,camera.z(:,idpress)-nanmean(camera.z(:,idpress)),'k','LineWidth',0.5)
% scatter(cam.time-cam.time(1),compcam.eta,1,'k')
title(['\textbf{Stereo}: $x$ = ',num2str(round(camera.xy(1,idpress))),' m, $y$ = ',num2str(round(camera.xy(2,idpress))),' m'],'interpreter','latex','fontsize',20);
% xlabel('$t$ (min)','interpreter','latex','fontsize',20);
ylabel('$z$ (m)','interpreter','latex','fontsize',20);
h1=gca;
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1*get(h1,'ticklength'));
set(h1,'fontsize',15);
datetick('x')
xlim([min(camera.time) max(lidar.time)])
ylim([-.2 0.2])
grid on

subplot(4,1,4)
plot(lidar.time,lidar.z(:,idpress)-nanmean(lidar.z(:,idpress)),'k','LineWidth',0.5)
% scatter(lidar.time-lidar.time(1),complidar.eta,1,'k')
title(['\textbf{Lidar}: $x$ = ',num2str(round(lidar.xy(1,idpress))),' m, $y$ = ',num2str(round(lidar.xy(2,idpress))),' m'],'interpreter','latex','fontsize',20);
xlabel('$t$ (min)','interpreter','latex','fontsize',20);
ylabel('$z$ (m)','interpreter','latex','fontsize',20);
h1=gca;
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1*get(h1,'ticklength'));
set(h1,'fontsize',15);
datetick('x')
xlim([min(lidar.time) max(lidar.time)])
ylim([-.2 0.2])
grid on

Sname1 = [figfolder,Tinfo.comp,'_timeseries_comparison_',compinst2find];
print(Sname1,'-dpng')

subplot(4,1,1)
xlim([min(lidar.time) (min(lidar.time)+((lidar.time(end)-lidar.time(1))/10))])

subplot(4,1,2)
xlim([min(lidar.time) (min(lidar.time)+((lidar.time(end)-lidar.time(1))/10))])

subplot(4,1,3)
xlim([min(lidar.time) (min(lidar.time)+((lidar.time(end)-lidar.time(1))/10))])

subplot(4,1,4)
xlim([min(lidar.time) (min(lidar.time)+((lidar.time(end)-lidar.time(1))/10))])

Sname1 = [figfolder,Tinfo.comp,'_timeseries_comparison_',compinst2find,'_zoomed'];
print(Sname1,'-dpng')

