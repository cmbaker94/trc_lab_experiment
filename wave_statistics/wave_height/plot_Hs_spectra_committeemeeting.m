% Plot Hs and sea-surface evolutino spectra for remote sensing and in situ
% gages.

% Plot Bulk Statics 
% Read in the lower resolution region and compute bulk statistics

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('/Users/cmbaker9/Documents/MTOOLS'))
addpath(genpath('/Users/cmbaker9/Documents/Research/Lab_Experiments/codes/insitu'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

spread = [0, 40];

for isprd = 1:length(spread)
    iarray = 1
        for ipg = 1:2
            if ipg == 1
                gages = {'wg4';'press6'}
            elseif ipg == 2
                gages = {'wg4';'press11'}
            end

% Trial info
Tinfo.Hs = 0.3;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
Tinfo.spread = spread(isprd);
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

%% STEP 3: remove sensors/areas with poor data

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

%% STEP 4: compare spectra

%     gages = {'wg4';'press2'};
            if iarray == 1
                insitu = sz;
                disp('sz')
            elseif iarray == 2
                insitu = is;
                disp('is')
            end
grab_tseries = 0;
spec = extract_lab_spectra(insitu,cam,lidar,gages,grab_tseries);

%% STEP 4: Let's plot the comparisons

% inan = find(isnan(spec1.press.See), 1);

figure('units','inches','position',[1 1 12 6],'Color','w');
% plot(spec.wg.f,spec.wg.See,'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.')

plot(spec.press.f,spec.press.See,'k','LineWidth',2)
hold on
plot(spec.cam.f,spec.cam.See,'r','LineWidth',2)
plot(spec.lidar.f,spec.lidar.See,'Color','b','LineWidth',2)

plot([0 2.5],[10^-20 10^-20],'Color','k','LineWidth',2)
plot([0 2.5],[10^-20 10^-20],'Color','k','LineWidth',2,'LineStyle','-.')

plot([0 2.5],[10^-2 10^-2],'Color',[.8 .8 .8],'LineWidth',0.1)
plot([0 2.5],[10^-3 10^-3],'Color',[.8 .8 .8],'LineWidth',0.1)
plot([0 2.5],[10^-4 10^-4],'Color',[.8 .8 .8],'LineWidth',0.1)
plot([0.5 0.5],[10^-5 10^-1],'Color',[.8 .8 .8],'LineWidth',0.1)
plot([1 1],[10^-5 10^-1],'Color',[.8 .8 .8],'LineWidth',0.1)
plot([1.5 1.5],[10^-5 10^-1],'Color',[.8 .8 .8],'LineWidth',0.1)
plot([2 2],[10^-5 10^-1],'Color',[.8 .8 .8],'LineWidth',0.1)
plot(spec.press.f,spec.press.See,'k','LineWidth',2)
% hold on
plot(spec.cam.f,spec.cam.See,'r','LineWidth',2)
plot(spec.lidar.f,spec.lidar.See,'Color','b','LineWidth',2)

errorbar(spec.wg.f(4),10^(-1.3),10^-1.85,'Color','k','LineWidth',2)%wg.See(3)-wg.Seec(3,1))

% h2 = legend('Offshore Wave Gage','Onshore Pressure Gage','Stereo Reconstruction','LiDAR');
% set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');

% title({['\textbf{Spectra Comparison}, $x$ = ',num2str(round(spec.press.xyz(1),2)),'m, $y$ = ',num2str(round(spec.press.xyz(2),2)),'m'],...
%     ['\textbf{Offshore ',gages{1},'}: $H_s$ = ',num2str(round(spec.wg.Hs,2)), ' m, \textbf{Onshore ',gages{2},'}: $H_s$ = ',num2str(round(spec.press.Hs,2)), 'm'],...
%     ['\textbf{Lidar}: $H_s$ = ',num2str(round(spec.lidar.Hs,2)),' m, \textbf{Stereo}: $H_s$ = ',num2str(round(spec.cam.Hs,2)),' m'],[]},'interpreter','latex','fontsize',20);
xlabel('$f$ (Hz)','interpreter','latex','fontsize',26);
ylabel('$S_{\eta\eta}$ (m$^2$/Hz)','interpreter','latex','fontsize',26);
h1=gca;
set(h1, 'YScale', 'log')
set(h1,'tickdir','out','xminortick','on','yminortick','on');
set(h1,'ticklength',1*get(h1,'ticklength'));
set(h1,'fontsize',26);
xlim([0 2.5])
ylim([10^-4.5 10^-1.5])

% grid on

sname = [figfolder,'TRM-',Tinfo.cam.tstart ,'_',Sname,'_spectra_comparison_',gages{2}];
print(sname,'-dpng')

%% plot yavg spectra

figure('units','inches','position',[1 1 12 6],'Color','w');
% plot(spec.wg.f,spec.wg.See_yavg,'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.')
% hold on
plot(spec.press.f,spec.press.See_yavg,'k','LineWidth',2)
hold on
plot(spec.cam.f,spec.cam.See_yavg,'r','LineWidth',2)
plot(spec.lidar.f,spec.lidar.See_yavg,'Color','b','LineWidth',2)

inotnan = find(~isnan(spec.press.See_ystd));
% fill wasn't showing up cause there were negative values
% temp = fliplr(spec.press.See_yavg(1,inotnan)-spec.press.See_ystd(1,inotnan));
% temp2 = temp<0
% temp(temp2)=0.00005;
% fill([spec.press.f(inotnan); flipud(spec.press.f(inotnan))],[spec.press.See_yavg(1,inotnan)+spec.press.See_ystd(1,inotnan) temp],'r','LineStyle','none')
fill([spec.press.f(inotnan); flipud(spec.press.f(inotnan))],[spec.press.See_yavg(1,inotnan)+spec.press.See_ystd(1,inotnan) fliplr(spec.press.See_yavg(1,inotnan)-spec.press.See_ystd(1,inotnan))],'k','LineStyle','none')
inotnan = spec.wg.f<5;
% fill([spec.wg.f(inotnan); flipud(spec.wg.f(inotnan))],[spec.wg.See_yavg(1,inotnan)+spec.wg.See_ystd(1,inotnan) fliplr(spec.wg.See_yavg(1,inotnan)-spec.wg.See_ystd(1,inotnan))],'r','LineStyle','none')
fill([spec.cam.f; flipud(spec.cam.f)],[spec.cam.See_yavg+spec.cam.See_ystd; flipud(spec.cam.See_yavg-spec.cam.See_ystd)],'r','LineStyle','none')
fill([spec.lidar.f; flipud(spec.lidar.f)],[spec.lidar.See_yavg+spec.lidar.See_ystd; flipud(spec.lidar.See_yavg-spec.lidar.See_ystd)],'b','LineStyle','none')
alpha(0.2)

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

xlabel('$f$ (Hz)','interpreter','latex','fontsize',26);
ylabel('$S_{\eta\eta}$ (m$^2$/Hz)','interpreter','latex','fontsize',26);
h1=gca;
set(h1, 'YScale', 'log')
set(h1,'tickdir','out','xminortick','on','yminortick','on');
set(h1,'ticklength',1*get(h1,'ticklength'));
set(h1,'fontsize',26);
xlim([0 2])
ylim([10^-4 10^-1])

% grid on

sname = [figfolder,'TRM-',Tinfo.cam.tstart ,'_',Sname,'_spectra_yavg_comparison_x',num2str(round(spec.press.xyz(1))),'m_nolegend'];
print(sname,'-dpng')

h2 = legend('Onshore Pressure Gage','Stereo Reconstruction','LiDAR');
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');

title({['\textbf{Spectra Comparison}, $x$ = ',num2str(round(spec.press.xyz(1),2)),'m, $y$ = ',num2str(round(spec.press.xyz(2),2)),'m'],...
    ['\textbf{Offshore wg}: $H_s$ = ',num2str(round(spec.wg.Hs_yavg,2)), ' m, \textbf{Onshore pressure gage}: $H_s$ = ',num2str(round(spec.press.Hs_yavg,2)), 'm'],...
    ['\textbf{Lidar}: $H_s$ = ',num2str(round(spec.lidar.Hs_yavg,2)),' m, \textbf{Stereo}: $H_s$ = ',num2str(round(spec.cam.Hs_yavg,2)),' m'],[]},'interpreter','latex','fontsize',20);


sname = [figfolder,'TRM-',Tinfo.cam.tstart ,'_',Sname,'_spectra_yavg_comparison_x',num2str(round(spec.press.xyz(1))),'m'];
print(sname,'-dpng')

        end

%% Plot Hs

ver = '';
figure('units','inches','position',[1 1 15 6],'Color','w');
plot(x,-h/5,'Color',[0.3 0.3 0.3],'LineWidth',1.5);
hold on
% S1
scatter(sz.xyz(1,1),sz.Hs(1),50,'k','LineWidth',2)
plot(cam.x(1,:),cam.Hs_yavg,'LineWidth',2,'Color','r')
plot(lidar.x,lidar.Hs_yavg,'LineWidth',2,'Color','b')
scatter(lidar.x,lidar.Hs_yavg,20,'b','fill')

ylim([-0.25 0.4])
xlim([18.5 37])

Hstemp = cam.Hs_yavg;
inotnan = find(~isnan(Hstemp));
camfill = fill([cam.x(1,inotnan) fliplr(cam.x(1,inotnan))],[Hstemp(inotnan)+cam.Hs_ystd(inotnan) fliplr(Hstemp(inotnan)-cam.Hs_ystd(inotnan))],'r','LineStyle','none')
alpha(camfill,0.2)

Hstemp = lidar.Hs_yavg;
inotnan = find(~isnan(Hstemp));
lidarfill = fill([lidar.x(inotnan) fliplr(lidar.x(inotnan))],[Hstemp(inotnan)+lidar.Hs_ystd(inotnan) fliplr(Hstemp(inotnan)-lidar.Hs_ystd(inotnan))],'b','LineStyle','none')
alpha(lidarfill,0.2)

bathyfill = fill([x; flipud(x)],[-h/5; -2*ones(size(h))],[0.7 0.7 0.7])
plot([18 40],[0 0],'Color','k','LineStyle','--','LineWidth',1);
alpha(bathyfill,0.6)

plot(cam.x(1,:),cam.Hs_yavg,'LineWidth',2,'Color','r')
plot(lidar.x,lidar.Hs_yavg,'LineWidth',2,'Color','b')


for i = 1:length(sz.Hs)
    if sz.Hs(i)<0.05
        sz.Hs(i) = NaN;
    end
    if ~isnan(sz.Hs(i))
        scatter(sz.xyz(i,1),sz.Hs(i),50,'k','LineWidth',2)
    end
end
for i = 1:length(is.Hs)
    if is.Hs(i)<0.05
        is.Hs(i) = NaN;
    end
    if ~isnan(is.Hs(i))
        scatter(is.xyz(i,1),is.Hs(i),50,'k','LineWidth',2)
    end
end

% analytic and bathymetry
% plot(analytic.x,wave.H,'Color',[1 0.6 0],'LineWidth',2,'LineStyle','-.')
plot(x,-h/5,'Color',[0.5 0.5 0.5],'LineWidth',2);
grid on
ylabel('$H_s$ (m)','interpreter','latex','fontsize',20);
xlabel('cross-shore (m)','interpreter','latex','fontsize',20);
h1=gca;
set(h1,'fontsize',26);
set(h1,'tickdir','out');%,'xminortick','on','yminortick','on');
set(h1,'ticklength',1*get(h1,'ticklength'));
set(h1,'ytick',[-0.2:0.1:0.6],'yticklabel',{'-0.2' '' '0' '' '0.2' '' '0.4' '' '0.6'});

sname1 = [figfolder,'TRM-',Tinfo.cam.tstart ,'_',Sname,'_Hsavg_comparison_nolegend'];
print(sname1,'-dpng')

h2 = legend('$h$/5','In Situ Instruments','Stereo Reconstruction','LiDAR')
% h2 = legend('$h$/5','Insitu SZ Array','Insitu IS Array','Stereo Reconstruction','LiDAR','Analytic Estimate','SWAN','SWAN h/5')
set(h2,'interpreter','latex','fontsize',16,'orientation','vertical','Location','northeast');

sname1 = [figfolder,'TRM-',Tinfo.cam.tstart ,'_',Sname,'_Hsavg_comparison'];
print(sname1,'-dpng')

close all

%% Plot Hs plan view

addpath('/Users/cmbaker9/Documents/MTOOLS/nctoolbox/cdm/utilities/graphics')

figure('units','inches','position',[1 1 6 10],'Color','w');
pcolorjw(cam.x,cam.y,cam.Hs);
hold on
scatter(sz.xyz(:,1),sz.xyz(:,2),50,sz.Hs,'fill','MarkerEdgeColor','k')
scatter(is.xyz(:,1),is.xyz(:,2),50,is.Hs,'fill','MarkerEdgeColor','k')
axis equal
shading flat
xlim([18.5 36])
ylim([-14 14])
    colormap(cmocean('thermal'));
    h2 = colorbar('northoutside');
    ylabel(h2,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
%     caxis([0.05 0.35])
    caxis([0.05 0.25])
    h1=gca;
    set(h1,'tickdir','in','xminortick','on','yminortick','on');
    set(h1,'ticklength',2*get(h1,'ticklength'));
    set(h1,'ydir','normal');
    set(h1,'fontsize',20);
    xlabel('$x~\mathrm{(m)}$','interpreter','latex','fontsize',20);
    ylabel('$y~\mathrm{(m)}$','interpreter','latex','fontsize',20);
    box on
    
    Sname1 = [figfolder,'TRM-',Tinfo.cam.tstart ,'_',Sname,'_recon_Hs'];
print(Sname1,'-dpng')

%% lidar

figure('units','inches','position',[1 1 6 10],'Color','w');
% figure('Color','w')
pcolorjw(lidar.x,lidar.y,lidar.Hs);%,100,'linestyle','none')%.*fliplr(beach2));
hold on
scatter(sz.xyz(:,1),sz.xyz(:,2),50,sz.Hs,'fill','MarkerEdgeColor','k')
hold on
scatter(is.xyz(:,1),is.xyz(:,2),50,is.Hs,'fill','MarkerEdgeColor','k')
axis equal
shading flat
xlim([18.5 36])
ylim([-14 14])
   
    colormap(cmocean('thermal'));
    h2 = colorbar('northoutside');
    ylabel(h2,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
    caxis([0.05 0.35])
%     caxis([0.05 0.25])
    h1=gca;
    set(h1,'tickdir','in','xminortick','on','yminortick','on');
    set(h1,'ticklength',2*get(h1,'ticklength'));
    set(h1,'ydir','normal');
    set(h1,'fontsize',20);
    xlabel('$x~\mathrm{(m)}$','interpreter','latex','fontsize',20);
    ylabel('$y~\mathrm{(m)}$','interpreter','latex','fontsize',20);
    box on
    
    Sname1 = [figfolder,'TRM-',Tinfo.cam.tstart ,'_',Sname,'_lidar_Hs'];
print(Sname1,'-dpng')

%%

%% cam
camplot_Hs = cam.Hs';
camplot_Hs(cam.x(1,:)<28.2,:) = NaN;

figure('units','inches','position',[1 1 10 6],'Color','w');
ax1 = axes('Position',[0.12 0.15 0.8 0.77]);
pcolorjw(cam.y,cam.x,camplot_Hs);%,100,'linestyle','none')%.*fliplr(beach2));
hold on
for i = 1:length(sz.Hs)
    if ~isnan(sz.Hs(i))
        scatter(sz.xyz(i,2),sz.xyz(i,1),100,sz.Hs(i),'fill','MarkerEdgeColor','k')
    end
end
for i = 1:length(is.Hs)
    if ~isnan(is.Hs(i))
        scatter(is.xyz(i,2),is.xyz(i,1),100,is.Hs(i),'fill','MarkerEdgeColor','k')
    end
end
axis equal
shading flat
grid on
ylim([18.5 35.01])
xlim([-12 12])
%     text(27,17,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
    colormap(cmocean('thermal'));
%     h2 = colorbar('eastoutside');
    cb = colorbar('Position', [0.885 0.15 0.035 0.67])
    cb.Ruler.MinorTick = 'on';
    text(12.8,19.3,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
%     ylabel(h2,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
    caxis([0 0.21])
    h1=gca;
    set(h1,'tickdir','out','xminortick','on','yminortick','on');
    set(h1,'ticklength',2*get(h1,'ticklength'));
    set(h1,'ydir','reverse');
    set(h1,'fontsize',18);
    ylabel('Cross-Shore (m)','interpreter','latex','fontsize',20);
    xlabel('Alongshore (m)','interpreter','latex','fontsize',20);
    box on
    
    Sname1 = [figfolder,'TRM-',Tinfo.cam.tstart ,'_',Sname,'_cam_Hs_flip'];
print(Sname1,'-dpng')

%% lidar

lidarplot_Hs = lidar.Hs';
lidarplot_Hs(lidar.x<23,:) = NaN;
lidarplot_Hs(:,lidar.y>13,:) = NaN;
lidarplot_Hs(:,lidar.y<-12.9,:) = NaN;

figure('units','inches','position',[1 1 10 6],'Color','w');
% figure('Color','w')

ax1 = axes('Position',[0.12 0.15 0.8 0.77]);

pcolorjw(lidar.y,lidar.x,lidarplot_Hs);%,100,'linestyle','none')%.*fliplr(beach2));
hold on
for i = 1:length(sz.Hs)
    if ~isnan(sz.Hs(i))
        scatter(sz.xyz(i,2),sz.xyz(i,1),100,sz.Hs(i),'fill','MarkerEdgeColor','k')
    end
end
for i = 1:length(is.Hs)
    if ~isnan(is.Hs(i))
        scatter(is.xyz(i,2),is.xyz(i,1),100,is.Hs(i),'fill','MarkerEdgeColor','k')
    end
end
axis equal
shading flat
grid on
ylim([18.5 35.01])
xlim([-12 12])
    
    %     ylim([min(y(:,1)) max(y(:,1))+2])
    %     end
    colormap(cmocean('thermal'));
    %     h2 = colorbar('eastoutside');
    cb = colorbar('Position', [0.885 0.15 0.035 0.67])
    cb.Ruler.MinorTick = 'on';
    text(12.8,19.3,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
    caxis([0.0 0.21])
    h1=gca;
    set(h1,'tickdir','out','xminortick','on','yminortick','on');
    set(h1,'ticklength',2*get(h1,'ticklength'));
    set(h1,'ydir','reverse');
    set(h1,'fontsize',18);
    % 	set(h1,'ytick',[-200:200:1200],'yticklabel',{'-200' '0' '200' '400' '600' '800' '1000' '1200'});
    %     set(h1,'xtick',[0:200:800],'xticklabel',{'0' '200' '400' '600' '800'});
    ylabel('Cross-Shore (m)','interpreter','latex','fontsize',20);
    xlabel('Alongshore (m)','interpreter','latex','fontsize',20);
    box on
    
    Sname1 = [figfolder,'TRM-',Tinfo.cam.tstart ,'_',Sname,'_lidar_Hs_flip'];
print(Sname1,'-dpng')

%% insitu

lidarplot_Hs = lidar.Hs';
lidarplot_Hs(lidar.x<23,:) = NaN;
lidarplot_Hs(:,lidar.y>13,:) = NaN;
lidarplot_Hs(:,lidar.y<-12.9,:) = NaN;

figure('units','inches','position',[1 1 10 6],'Color','w');
% figure('Color','w')

ax1 = axes('Position',[0.12 0.15 0.8 0.77]);

scatter(sz.xyz(1,2),sz.xyz(1,1),150,0,'fill','MarkerEdgeColor','k')
hold on
scatter(sz.xyz(end,2),sz.xyz(end,1),150,0,'sq','fill','MarkerEdgeColor','k')
for i = 1:length(sz.Hs)
    if ~isnan(sz.Hs(i))
        if sz.xyz(i,1)> 22
            scatter(sz.xyz(i,2),sz.xyz(i,1),150,sz.Hs(i),'fill','MarkerEdgeColor','k')
        elseif sz.xyz(i,1) <= 22
            scatter(sz.xyz(i,2),sz.xyz(i,1),150,sz.Hs(i),'sq','fill','MarkerEdgeColor','k')
        end
    end
end
for i = 1:length(is.Hs)
    if ~isnan(is.Hs(i))
        if is.xyz(i,1)> 22
            scatter(is.xyz(i,2),is.xyz(i,1),150,is.Hs(i),'fill','MarkerEdgeColor','k')
        elseif is.xyz(i,1) <= 22
            scatter(is.xyz(i,2),is.xyz(i,1),150,is.Hs(i),'sq','fill','MarkerEdgeColor','k')
        end
    end
end
axis equal
shading flat
grid on
ylim([18.5 35.01])
xlim([-12 12])
    
    %     ylim([min(y(:,1)) max(y(:,1))+2])
    %     end
    colormap(cmocean('thermal'));
    %     h2 = colorbar('eastoutside');
    cb = colorbar('Position', [0.885 0.15 0.035 0.67])
    cb.Ruler.MinorTick = 'on';
    text(12.8,19.3,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',24);
    caxis([0.0 0.32])
    h1=gca;
    set(h1,'tickdir','out','xminortick','on','yminortick','on');
    set(h1,'ticklength',2*get(h1,'ticklength'));
    set(h1,'ydir','reverse');
    set(h1,'fontsize',18);
    % 	set(h1,'ytick',[-200:200:1200],'yticklabel',{'-200' '0' '200' '400' '600' '800' '1000' '1200'});
    %     set(h1,'xtick',[0:200:800],'xticklabel',{'0' '200' '400' '600' '800'});
    ylabel('Cross-Shore (m)','interpreter','latex','fontsize',24);
    xlabel('Alongshore (m)','interpreter','latex','fontsize',24);
    box on
    h2 = legend('Collocated ADV \& Pressure Gage','Wire Resistance Gage');
    set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','southeast');

    
    Sname1 = [figfolder,'TRM-',Tinfo.cam.tstart ,'_',Sname,'_insitu'];
print(Sname1,'-dpng')

end