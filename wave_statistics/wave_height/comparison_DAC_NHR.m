% Plot Hs and sea-surface evolutino spectra for remote sensing and in situ
% gages.

% Plot Bulk Statics
% Read in the lower resolution region and compute bulk statistics

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trial info
Tinfo.Hs = 0.25;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
Tinfo.spread = 20;

Tinfo = trial_files(Tinfo);

% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

%     [camera.time] = cam_time(Tinfo);

%% STEP 1: Load data

datapath    = 'E:\'; % general path and names
% Tinfo.comp = ['Hs',num2str(Tinfo.Hs*100),'_Tp',num2str(Tinfo.Tp),'_tide',num2str(Tinfo.tide*100),'_spread',num2str(Tinfo.spread)];
procpath = [datapath,'data\processed\conditions\',Tinfo.comp,'\',Tinfo.cam.tstart,'\'];
% procpath =

% bathymetry
bathymetry = load('E:\data\processed\lidar\Riegl\TRC_bathymetry_estimate_line.mat');
h = Tinfo.tide-bathymetry.h;
x = bathymetry.xp;

% Time series 1
% trange = Sname;
Sname = Tinfo.datarange;
load([Tinfo.savefolder,'\insitu_Hs_spectra.mat'], 'sz');
load([Tinfo.savefolder,'\insitu_Hs_spectra.mat'], 'is');
load([Tinfo.savefolder,'\velodyne_Hs_spectra.mat'], 'lidar');
load([Tinfo.savefolder,'\stereo_Hs_spectra.mat'], 'cam');

%% STEP 2: Define figure folders

% figure folder
% fssubfolder = datestr(date,'yy-mm-dd');
% % figfolder   = [datapath,'figures\meas_comp\',savefolder,'\frames_07200-11999\'];
% eval(['!mkdir ',Tinfo.figfolder])
figfolder   = Tinfo.figfolder;
eval(['!mkdir ',figfolder])

%% STEP 3: remove sensors\areas with poor data

% remove side walls & offshore
% LiDAR
lidar.Hs(lidar.y<-12.8,:)=NaN;
lidar.Hs(lidar.y>12.8,:)=NaN;
lidar.Hs(:,lidar.x<24) = NaN;
lidar.See(lidar.y<-12.8,:,:)=NaN;
lidar.See(lidar.y>12.8,:)= NaN;
lidar.See(:,lidar.x<24,:) = NaN;

% Stereo
cam.Hs(cam.y(:,1)<-12.9,:)=NaN;
cam.Hs(cam.y(:,1)>13,:)=NaN;
cam.Hs(:,cam.x(1,:)<28.35) = NaN;
cam.See(cam.y(:,1)<-12.9,:,:)=NaN;
cam.See(cam.y(:,1)>13,:,:)=NaN;
cam.See(:,cam.x(1,:)<28.35,:) = NaN;

% avg and stdev
lidar.Hs_yavg = nanmean(lidar.Hs,1);
cam.Hs_yavg = nanmean(cam.Hs,1);
lidar.See_yavg = squeeze(nanmean(lidar.See,1));
cam.See_yavg = squeeze(nanmean(cam.See,1));

lidar.Hs_ystd = nanstd(lidar.Hs,1);
cam.Hs_ystd = nanstd(cam.Hs,1);
lidar.See_ystd = squeeze(nanstd(lidar.See,1));
cam.See_ystd = squeeze(nanstd(cam.See,1));

for i = 1:length(sz.Hs)
    if sz.Hs(i)<0.05
        sz.Hs(i) = NaN;
        sz.See(i,:) = NaN;
    end
end
for i = 1:length(is.Hs)
    if is.Hs(i)<0.05
        is.Hs(i) = NaN;
        is.See(i,:) = NaN;
    end
end

%%

szlocsave = 'E:\data\processed\conditions\Hs25_Tp2_tide107_spread20\08-30-2018-1905UTC\time_1920-1929\';
islocsave = 'E:\data\processed\conditions\Hs25_Tp2_tide107_spread20\09-06-2018-2049UTC\time_2104-2113\';

szloc = 'E:\data\processed\insitu\08-30-2018-1904UTC\EN_nh_recon\';
isloc = 'E:\data\processed\insitu\09-06-2018-2048UTC\EN_nh_recon\';

for i = 1:12
    fileID = fopen([szloc,'press',num2str(i),'.txt'],'r');
    formatSpec = '%f';
    eval(['sznl.data.eta',num2str(i),' = fscanf(fileID,formatSpec);'])
end

for i = 1:12
    fileID = fopen([isloc,'press',num2str(i),'.txt'],'r');
    formatSpec = '%f';
    eval(['isnl.data.eta',num2str(i),' = fscanf(fileID,formatSpec);'])
end

%%

% spectral analysis details
sdet.WL              = 2^5;
sdet.OL              = 2^4;
sdet.frange     = [0.25, 1.2];
sdet.Hz         = 100;

istart  = 15*60*sdet.Hz;
iend    = ((15+10)*60*sdet.Hz)-1;

% sz
enames  = fieldnames(sznl.data);
pnames  = sz.inst;
for i = 1:12
    pn(i) = str2double(pnames{i}(6:end));
end

for i = 1:length(enames)
    eval(['eta = sznl.data.',enames{i},'(istart:iend);']);
    sznl.inst{i} = enames{i};
    ixyz = find(str2double(enames{i}(4:end))==pn)
    sznl.xyz(i,:) = sz.xyz(ixyz,:);
    [sznl.See(i,:),sznl.f(i,:),sznl.Seec(i,:,:)]   = pwelch(eta,sdet.WL*sdet.Hz,sdet.OL*sdet.Hz,[],sdet.Hz,'ConfidenceLevel',0.95);
    sznl.Hs_4std(i) = 4*nanstd(eta);
    [sznl.Hs(i),sznl.Tp(i)]       = spec2HsTp(sznl.f(i,:),sznl.See(i,:),sdet.frange);
end

% is
enames  = fieldnames(isnl.data);
pnames  = is.inst;
for i = 1:12
    pn(i) = str2double(pnames{i}(6:end));
end

for i = 1:length(enames)
    eval(['eta = isnl.data.',enames{i},'(istart:iend);']);
    isnl.inst{i} = enames{i};
    ixyz = find(str2double(enames{i}(4:end))==pn)
    isnl.xyz(i,:) = is.xyz(ixyz,:);
    [isnl.See(i,:),isnl.f(i,:),isnl.Seec(i,:,:)]   = pwelch(eta,sdet.WL*sdet.Hz,sdet.OL*sdet.Hz,[],sdet.Hz,'ConfidenceLevel',0.95);
    isnl.Hs_4std(i) = 4*nanstd(eta);
    [isnl.Hs(i),isnl.Tp(i)]       = spec2HsTp(isnl.f(i,:),isnl.See(i,:),sdet.frange);
end


%% STEP 4: compare spectra

for iarray = 1:2
    for ipg = 1:2
        if ipg == 1
            gages = {'wg4';'press6'}
        elseif ipg == 2
            gages = {'wg4';'press11'}
        end
        
        if iarray == 1
            insitu = sz;
            insitu_nl = sznl;
            disp('sz')
        elseif iarray == 2
            insitu = is;
            insitu_nl = isnl;
            disp('is')
        end
        grab_tseries = 0;
        spec = extract_lab_spectra(insitu,cam,lidar,gages,grab_tseries);
        specnl = extract_press_spectra_ENcorrection(insitu_nl,gages);
%         eval(['save -v7.3 ',Tinfo.savefolder,'Hs_allinst',' sz',' is'])%,' sdet'])
        
        %% STEP 4: Let's plot the comparisons
        
        % inan = find(isnan(spec1.press.See), 1);
        
        figure('units','inches','position',[1 1 12 6],'Color','w');
        plot(spec.wg.f,spec.wg.See,'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','-.')
        hold on
        plot(spec.press.f,spec.press.See,'k','LineWidth',2)
        plot(specnl.press.f,specnl.press.See,'r','LineWidth',2,'LineStyle','-.')
        plot(spec.cam.f,spec.cam.See,'b','LineWidth',2)
        plot(spec.lidar.f,spec.lidar.See,'Color',[0 0.5 0],'LineWidth',2)
        
        plot([0 2.5],[10^-20 10^-20],'Color','k','LineWidth',2)
        plot([0 2.5],[10^-20 10^-20],'Color','k','LineWidth',2,'LineStyle','-.')
        
        plot([0 2.5],[10^-2 10^-2],'Color',[.8 .8 .8],'LineWidth',0.1)
        plot([0 2.5],[10^-3 10^-3],'Color',[.8 .8 .8],'LineWidth',0.1)
        plot([0 2.5],[10^-4 10^-4],'Color',[.8 .8 .8],'LineWidth',0.1)
        plot([0.5 0.5],[10^-5 10^-1],'Color',[.8 .8 .8],'LineWidth',0.1)
        plot([1 1],[10^-5 10^-1],'Color',[.8 .8 .8],'LineWidth',0.1)
        plot([1.5 1.5],[10^-5 10^-1],'Color',[.8 .8 .8],'LineWidth',0.1)
        plot([2 2],[10^-5 10^-1],'Color',[.8 .8 .8],'LineWidth',0.1)
        
        errorbar(spec.wg.f(4),10^(-1.3),10^-1.85,'Color','k','LineWidth',2)%wg.See(3)-wg.Seec(3,1))
        
        h2 = legend('Offshore Wave Gage','Press: DAC','Press: NL','Stereo Reconstruction','LiDAR');
        set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
        
        title({['\textbf{Spectra Comparison}, $x$ = ',num2str(round(spec.press.xyz(1),2)),'m, $y$ = ',num2str(round(spec.press.xyz(2),2)),'m'],...
            ['\textbf{Offshore ',gages{1},'}: $H_s$ = ',num2str(round(spec.wg.Hs,2)), ' m, \textbf{Onshore ',gages{2},'}: $H_s$ = ',num2str(round(spec.press.Hs,2)), 'm'],...
            ['\textbf{Lidar}: $H_s$ = ',num2str(round(spec.lidar.Hs,2)),' m, \textbf{Stereo}: $H_s$ = ',num2str(round(spec.cam.Hs,2)),' m'],[]},'interpreter','latex','fontsize',20);
        xlabel('$f$ (Hz)','interpreter','latex','fontsize',22);
        ylabel('$S_{\eta\eta}$ (m$^2$/Hz)','interpreter','latex','fontsize',22);
        h1=gca;
        set(h1, 'YScale', 'log')
        set(h1,'tickdir','out','xminortick','on','yminortick','on');
        set(h1,'ticklength',1*get(h1,'ticklength'));
        set(h1,'fontsize',15);
        xlim([0 2.5])
        ylim([10^-5 10^-1])
        
        % grid on
        if iarray == 1
            sname = [figfolder,'TRM-',Tinfo.cam.tstart ,'_',Sname,'_spectra_comparison_surfzonearray_',gages{2},'_nl'];
        elseif iarray == 2
            sname = [figfolder,'TRM-',Tinfo.cam.tstart ,'_',Sname,'_spectra_comparison_innershelfarray_',gages{2},'_nl'];
        end
        print(sname,'-dpng')
        
    end
end
%% Plot Hs

ver = '';
figure('units','inches','position',[1 1 15 6],'Color','w');
plot(x,-h/5,'Color',[0.3 0.3 0.3],'LineWidth',1.5);
hold on
% S1
scatter(sz.xyz(1,1),sz.Hs(1),40,'k','LineWidth',2)
scatter(sznl.xyz(1,1),sznl.Hs(1),20,'r','LineWidth',2)
scatter(sznl.xyz(1,1),sznl.Hs_4std(1),20,'MarkerEdgeColor',[235, 168, 52]/256,'LineWidth',2)
plot(cam.x(1,:),cam.Hs_yavg,'LineWidth',2,'Color','b')
plot(lidar.x,lidar.Hs_yavg,'LineWidth',2,'Color',[0 0.5 0])

ylim([-0.25 0.4])
xlim([18.5 37])

Hstemp = cam.Hs_yavg;
inotnan = find(~isnan(Hstemp));
fill([cam.x(1,inotnan) fliplr(cam.x(1,inotnan))],[Hstemp(inotnan)+cam.Hs_ystd(inotnan) fliplr(Hstemp(inotnan)-cam.Hs_ystd(inotnan))],'b','LineStyle','none')

Hstemp = lidar.Hs_yavg;
inotnan = find(~isnan(Hstemp));
fill([lidar.x(inotnan) fliplr(lidar.x(inotnan))],[Hstemp(inotnan)+lidar.Hs_ystd(inotnan) fliplr(Hstemp(inotnan)-lidar.Hs_ystd(inotnan))],[0 0.5 0],'LineStyle','none')
alpha(0.2)

bathyfill = fill([x; flipud(x)],[-h/5; -2*ones(size(h))],[0.7 0.7 0.7])
plot([18 40],[0 0],'Color','k','LineStyle','--','LineWidth',1);
alpha(bathyfill,0.6)

plot(cam.x(1,:),cam.Hs_yavg,'LineWidth',2,'Color','b')
plot(lidar.x,lidar.Hs_yavg,'LineWidth',2,'Color',[0 0.5 0])


for i = 1:length(sz.Hs)
    if sz.Hs(i)<0.05
        sz.Hs(i) = NaN;
    end
    if ~isnan(sz.Hs(i))
        scatter(sz.xyz(i,1),sz.Hs(i),40,'k','LineWidth',2)
    end
end
for i = 1:length(is.Hs)
    if is.Hs(i)<0.05
        is.Hs(i) = NaN;
    end
    if ~isnan(is.Hs(i))
        scatter(is.xyz(i,1),is.Hs(i),40,'k','LineWidth',2)
    end
end

for i = 1:length(sznl.Hs)
    if sznl.Hs(i)<0.05
        sznl.Hs(i) = NaN;
    end
    if ~isnan(sznl.Hs(i))
        scatter(sznl.xyz(i,1),sznl.Hs(i),20,'r','LineWidth',2)
        scatter(sznl.xyz(i,1),sznl.Hs_4std(i),20,'MarkerEdgeColor',[235, 168, 52]/256,'LineWidth',2)
    end
end
for i = 1:length(isnl.Hs)
    if isnl.Hs(i)<0.05
        isnl.Hs(i) = NaN;
    end
    if ~isnan(isnl.Hs(i))
        scatter(isnl.xyz(i,1),isnl.Hs(i),20,'r','LineWidth',2)
        scatter(isnl.xyz(i,1),isnl.Hs_4std(i),20,'MarkerEdgeColor',[235, 168, 52]/256,'LineWidth',2)
    end
end

% analytic and bathymetry
% plot(analytic.x,wave.H,'Color',[1 0.6 0],'LineWidth',2,'LineStyle','-.')
plot(x,-h/5,'Color',[0.5 0.5 0.5],'LineWidth',2);
h2 = legend('$h$/5','In Situ: DAC','In Situ: NL','In Situ: NL $4*\sigma_{\eta}$','Stereo Reconstruction','LiDAR')


% h2 = legend('$h$/5','Insitu SZ Array','Insitu IS Array','Stereo Reconstruction','LiDAR','Analytic Estimate','SWAN','SWAN h/5')
set(h2,'interpreter','latex','fontsize',16,'orientation','vertical','Location','northeast');
grid on
ylabel('$H_s$ (m)','interpreter','latex','fontsize',20);
xlabel('cross-shore (m)','interpreter','latex','fontsize',20);
h1=gca;
set(h1,'fontsize',18);
set(h1,'tickdir','out');%,'xminortick','on','yminortick','on');
set(h1,'ticklength',1*get(h1,'ticklength'));
set(h1,'ytick',[-0.2:0.1:0.6],'yticklabel',{'-0.2' '' '0' '' '0.2' '' '0.4' '' '0.6'});

sname1 = [figfolder,'TRM-',Tinfo.cam.tstart ,'_',Sname,'_Hsavg_comparison_withnl'];
print(sname1,'-dpng')
% close all

%%

figure('units','inches','position',[1 1 10 6],'Color','w');
td = 0.01:0.01:(length(sznl.data.eta6)/100);
plot(td,sznl.data.eta6,'LineWidth',0.5,'LineStyle','-','Color','k')
hold on
box on
% ylim([-0.2 0.45])
xlim([0 max(td)])
h1=gca;
set(h1,'tickdir','in','xminortick','off','yminortick','off');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',20);
% set(h1,'xtick',ynum,'xticklabel',{''});
ylabel('$\eta$ (m)','interpreter','latex','fontsize',24);
xlabel('time (s)','interpreter','latex','fontsize',24);
% title(ax1,'$\sigma_{\theta}=0^{\circ}$','interpreter','latex','fontsize',28);
% text(ax1,ynum(1)+0.8, 0.41,'(a)','interpreter','latex','fontsize',24);
sname1 = [figfolder,'surfzonearray_outersz_timeseries'];
print(sname1,'-dpng')

figure('units','inches','position',[1 1 10 6],'Color','w');
td = 0.01:0.01:(length(sznl.data.eta11)/100);
plot(td,sznl.data.eta11,'LineWidth',0.5,'LineStyle','-','Color','k')
hold on
box on
% ylim([-0.2 0.45])
xlim([0 max(td)])
h1=gca;
set(h1,'tickdir','in','xminortick','off','yminortick','off');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',20);
% set(h1,'xtick',ynum,'xticklabel',{''});
ylabel('$\eta$ (m)','interpreter','latex','fontsize',24);
xlabel('time (s)','interpreter','latex','fontsize',24);
% title(ax1,'$\sigma_{\theta}=0^{\circ}$','interpreter','latex','fontsize',28);
% text(ax1,ynum(1)+0.8, 0.41,'(a)','interpreter','latex','fontsize',24);
sname1 = [figfolder,'surfzonearray_innersz_timeseries'];
print(sname1,'-dpng')