% Figure of tank geomoetry and in situ gages. 

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\codes\insitu'))
addpath(genpath('E:\codes\trc_lab_experiment\toolbox'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

sprd = [0, 40];
tstart = [125.75+2, 57.25+2-1.125];%[125.6, 152.2]; % start time for timeseries for 0 and 40 deg
tlen = 40;%20 % timeseries seconds
%57.25 or 511.375.... or 405.375 or 543.375 or 511.350
%% STEP 2: Define figure folders

% figure folder
datapath = 'E:\';
fssubfolder = datestr(date,'yy-mm-dd');
figfolder   = [datapath,'figures\methods_manuscript\testing\',fssubfolder,'\'];
eval(['!mkdir ',figfolder])


%% Loop through two trials

for i = 1:length(sprd)
    
    %% STEP 0: Prepare file structures, folders, etc.
    Tinfo.spread = sprd(i);
    Tinfo.Hs = 0.3;
    Tinfo.Tp = 2;
    Tinfo.tide = 1.07;
    
    Tinfo = trial_files(Tinfo);
    
    % general path and names
    Tinfo.cam = TRC_camera_info(Tinfo.cam);
    
    % Data and figure storage
    [Tinfo] = wc_comp_store(Tinfo);
    
    %% STEP 1: Load Stereo Reconstruction Data
    
    % cam         = load([Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_res',num2str((Tinfo.cam.dx)*100),'cm.mat']);
    % cam.time    = Tinfo.cam.timevec; % Establish timeseries range
    [camera.time] = cam_time(Tinfo);
    
    
    
    %% Load time series
    
    [press,wg,camera,lidar] = extract_fulltimeseries_insitu(Tinfo,camera);
    
    
    %% Extract time and locations of interest
    
    time = camera.time(tstart(i)*8:tstart(i)*8+tlen*8);
    tseries.name = {'11','6'};
    tseries = extract_timestack_insitu(time,camera,press,lidar,tseries);
    eval(['T',num2str(Tinfo.spread),' = tseries;'])
    
    clear tseries camera lidar press wg Tinfo
    
end

%% Create time string for plotting

It = T0.Ito-(tstart(1)+tlen/2);
Ct = T0.Cto-(tstart(1)+tlen/2);
Lt = T0.Lto-(tstart(1)+tlen/2);

ylab = strsplit(num2str(-tlen/2:2:tlen/2));
% clear ylab
ynum = -tlen/2:2:tlen/2;
% if want every other labeled
% for i = 1:length(ynum)
%     if rem(ynum(i),2) == 0
%         ylab{i} = num2str(ynum(i))
%     else
%         ylab{i} = ['']
%     end
% end

%% Movmedian
pts = 3;

T0.off.Czf = movmedian(T0.off.Cz,pts);
T40.off.Czf = movmedian(T40.off.Cz,pts);
T0.on.Czf = movmedian(T0.on.Cz,pts);
T40.on.Czf = movmedian(T40.on.Cz,pts);

T0.off.Lzf = movmedian(T0.off.Lz,pts);
T40.off.Lzf = movmedian(T40.off.Lz,pts);
T0.on.Lzf = movmedian(T0.on.Lz,pts);
T40.on.Lzf = movmedian(T40.on.Lz,pts);

%% 

figure('units','inches','position',[1 1 18 10],'color','w')

ax1 = axes('Position',[0.095 0.48 0.42 0.45]);
plot([min(Ct) max(Ct)],[0 0],'LineWidth',1.5,'LineStyle','-.','Color',[0.7 0.7 0.7])
hold on
fill([-0.06 -0.06 0.06 0.06],[-0.3 0.5 0.5 -0.3],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
scatter(It,T0.off.Iz,10,[0.7 0.7 0.7],'fill')
% scatter(It,T0.off.Ifz,4,'k','fill')
scatter(Ct,T0.off.Cz,20,'r','fill')
scatter(Lt,T0.off.Lz,20,'b','fill')
plot(It,T0.off.Iz,'Color',[0.7 0.7 0.7],'LineWidth',1)
plot(It,T0.off.Ifz,'k','LineWidth',2)
plot(Ct,T0.off.Cz,'r','LineWidth',1)
plot(Lt,T0.off.Lz,'b','LineWidth',1)
% grid on
box on
ylim([-0.2 0.45])
xlim(round([min(Ct) max(Ct)]))
h1=gca;
set(h1,'tickdir','in','xminortick','off','yminortick','off');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',20);
set(h1,'xtick',ynum,'xticklabel',{''});
ylabel({'Outter Edge of Surf Zone';'$\eta$ (m)'},'interpreter','latex','fontsize',24);
% xlabel('time (s)','interpreter','latex','fontsize',24);
title(ax1,'$\sigma_{\theta}=0^{\circ}$','interpreter','latex','fontsize',28);
text(ax1,ynum(1)+0.5, 0.41,'(a)','interpreter','latex','fontsize',24);


ax2 = axes('Position',[0.095 0.09 0.42 0.35]);
plot([min(Ct) max(Ct)],[0 0],'LineWidth',1.5,'LineStyle','-.','Color',[0.7 0.7 0.7])
hold on
fill([-0.06 -0.06 0.06 0.06],[-0.3 0.5 0.5 -0.3],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
scatter(It,T0.on.Iz,10,[0.7 0.7 0.7],'fill')
% scatter(It,T0.on.Ifz,4,'k','fill')
scatter(Ct,T0.on.Cz,20,'r','fill')
scatter(Lt,T0.on.Lz,20,'b','fill')
plot(It,T0.on.Iz,'Color',[0.7 0.7 0.7],'LineWidth',1)
plot(It,T0.on.Ifz,'k','LineWidth',2)
plot(Ct,T0.on.Cz,'r','LineWidth',1)
plot(Lt,T0.on.Lz,'b','LineWidth',1)
% grid on
box on
ylim([-.15 .2])
xlim(round([min(Ct) max(Ct)]))
h1=gca;
set(h1,'tickdir','in','xminortick','off','yminortick','off');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',20);
set(h1,'xtick',ynum,'xticklabel',ylab);
xlabel('Time (s)','interpreter','latex','fontsize',24);
ylabel({'Inner Surf Zone';'$\eta$ (m)'},'interpreter','latex','fontsize',24);
text(ax2,ynum(1)+0.5, 0.17,'(b)','interpreter','latex','fontsize',24);


ax3 = axes('Position',[0.54 0.48 0.42 0.45]);
plot(It,T40.off.Ifz,'k','LineWidth',1)
hold on
scatter(It,T40.off.Iz,10,[0.7 0.7 0.7],'fill')
scatter(Ct,T40.off.Cz,20,'r','fill')
scatter(Lt,T40.off.Lz,20,'b','fill')
plot([min(Ct) max(Ct)],[0 0],'LineWidth',1.5,'LineStyle','-.','Color',[0.7 0.7 0.7])
fill([-0.06 -0.06 0.06 0.06],[-0.3 0.5 0.5 -0.3],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
fill([-0.06 -0.06 0.06 0.06]-6.625,[-0.3 0.5 0.5 -0.3],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
scatter(It,T40.off.Iz,10,[0.7 0.7 0.7],'fill')
% scatter(It,T40.off.Ifz,4,'k','fill')
scatter(Ct,T40.off.Cz,20,'r','fill')
scatter(Lt,T40.off.Lz,20,'b','fill')
plot(It,T40.off.Iz,'Color',[0.7 0.7 0.7],'LineWidth',1)
plot(It,T40.off.Ifz,'k','LineWidth',2)
plot(Ct,T40.off.Cz,'r','LineWidth',1)
plot(Lt,T40.off.Lz,'b','LineWidth',1)
% grid on
box on
ylim([-.2 .45])
xlim(round([min(Ct) max(Ct)]))
h1=gca;
set(h1,'tickdir','in','xminortick','off','yminortick','off');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',20);
set(h1,'yticklabel',[]);
set(h1,'xtick',ynum,'xticklabel',{''});
title(ax3,'$\sigma_{\theta}=40^{\circ}$','interpreter','latex','fontsize',28);
text(ax3,ynum(1)+0.5, 0.41,'(c)','interpreter','latex','fontsize',24);
h2 = legend('press, DAC','press, NHR','Stereo','Lidar','interpreter','latex','fontsize',22)
% set(h2, 'Position', [0.68 0.86 0.10 0.05],'orientation','vertical')
% set(h2, 'Position', [0.8614 0.778 0.08 0.15],'orientation','vertical')
set(h2, 'Position', [0.862 0.778 0.08 0.15],'orientation','vertical')

ax4 = axes('Position',[0.54 0.09 0.42 0.35]);
plot([min(Ct) max(Ct)],[0 0],'LineWidth',1.5,'LineStyle','-.','Color',[0.7 0.7 0.7])
hold on
fill([-0.06 -0.06 0.06 0.06],[-0.3 0.5 0.5 -0.3],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
fill([-0.06 -0.06 0.06 0.06]-6.625,[-0.3 0.5 0.5 -0.3],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
scatter(It,T40.on.Iz,10,[0.7 0.7 0.7],'fill')
% scatter(It,T40.on.Ifz,4,'k','fill')
scatter(Ct,T40.on.Cz,20,'r','fill')
scatter(Lt,T40.on.Lz,20,'b','fill')
plot(It,T40.on.Iz,'Color',[0.7 0.7 0.7],'LineWidth',1)
plot(It,T40.on.Ifz,'k','LineWidth',2)
plot(Ct,T40.on.Cz,'r','LineWidth',1)
plot(Lt,T40.on.Lz,'b','LineWidth',1)
% grid on
box on
ylim([-.15 .2])
xlim(round([min(Ct) max(Ct)]))
h1=gca;
set(h1,'tickdir','in','xminortick','off','yminortick','off');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',20);
set(h1,'yticklabel',[]);
set(h1,'xtick',ynum,'xticklabel',ylab);
xlabel('Time (s)','interpreter','latex','fontsize',24);
text(ax4,ynum(1)+0.5, 0.17,'(d)','interpreter','latex','fontsize',24);

sname = [figfolder,'timeseries_nofilt'];
print(sname,'-dpng')
 
%% Store mat files

figure('units','inches','position',[1 1 18 10],'color','w')

ax1 = axes('Position',[0.095 0.48 0.42 0.45]);
plot([min(Ct) max(Ct)],[0 0],'LineWidth',1.5,'LineStyle','-.','Color',[0.7 0.7 0.7])
hold on
fill([-0.06 -0.06 0.06 0.06],[-0.3 0.5 0.5 -0.3],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
scatter(It,T0.off.Iz,10,[0.7 0.7 0.7],'fill')
% scatter(It,T0.off.Ifz,4,'k','fill')
scatter(Ct,T0.off.Czf,20,'r','fill')
scatter(Lt,T0.off.Lzf,20,'b','fill')
plot(It,T0.off.Iz,'Color',[0.7 0.7 0.7],'LineWidth',1)
plot(It,T0.off.Ifz,'k','LineWidth',2)
plot(Ct,T0.off.Czf,'r','LineWidth',1)
plot(Lt,T0.off.Lzf,'b','LineWidth',1)
% grid on
box on
ylim([-0.2 0.45])
xlim(round([min(Ct) max(Ct)]))
h1=gca;
set(h1,'tickdir','in','xminortick','off','yminortick','off');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',20);
set(h1,'xtick',ynum,'xticklabel',{''});
ylabel({'Outter Edge of Surf Zone';'$\eta$ (m)'},'interpreter','latex','fontsize',24);
% xlabel('time (s)','interpreter','latex','fontsize',24);
title(ax1,'$\sigma_{\theta}=0^{\circ}$','interpreter','latex','fontsize',28);
text(ax1,ynum(1)+0.5, 0.41,'(a)','interpreter','latex','fontsize',24);


ax2 = axes('Position',[0.095 0.09 0.42 0.35]);
plot([min(Ct) max(Ct)],[0 0],'LineWidth',1.5,'LineStyle','-.','Color',[0.7 0.7 0.7])
hold on
fill([-0.06 -0.06 0.06 0.06],[-0.3 0.5 0.5 -0.3],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
scatter(It,T0.on.Iz,10,[0.7 0.7 0.7],'fill')
% scatter(It,T0.on.Ifz,4,'k','fill')
scatter(Ct,T0.on.Czf,20,'r','fill')
scatter(Lt,T0.on.Lzf,20,'b','fill')
plot(It,T0.on.Iz,'Color',[0.7 0.7 0.7],'LineWidth',1)
plot(It,T0.on.Ifz,'k','LineWidth',2)
plot(Ct,T0.on.Czf,'r','LineWidth',1)
plot(Lt,T0.on.Lzf,'b','LineWidth',1)
% grid on
box on
ylim([-.15 .2])
xlim(round([min(Ct) max(Ct)]))
h1=gca;
set(h1,'tickdir','in','xminortick','off','yminortick','off');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',20);
set(h1,'xtick',ynum,'xticklabel',ylab);
xlabel('Time (s)','interpreter','latex','fontsize',24);
ylabel({'Inner Surf Zone';'$\eta$ (m)'},'interpreter','latex','fontsize',24);
text(ax2,ynum(1)+0.5, 0.17,'(b)','interpreter','latex','fontsize',24);


ax3 = axes('Position',[0.54 0.48 0.42 0.45]);
plot(It,T40.off.Ifz,'k','LineWidth',1)
hold on
scatter(It,T40.off.Iz,10,[0.7 0.7 0.7],'fill')
scatter(Ct,T40.off.Czf,20,'r','fill')
scatter(Lt,T40.off.Lzf,20,'b','fill')
plot([min(Ct) max(Ct)],[0 0],'LineWidth',1.5,'LineStyle','-.','Color',[0.7 0.7 0.7])
fill([-0.06 -0.06 0.06 0.06],[-0.3 0.5 0.5 -0.3],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
fill([-0.06 -0.06 0.06 0.06]-6.625,[-0.3 0.5 0.5 -0.3],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
scatter(It,T40.off.Iz,10,[0.7 0.7 0.7],'fill')
% scatter(It,T40.off.Ifz,4,'k','fill')
scatter(Ct,T40.off.Czf,20,'r','fill')
scatter(Lt,T40.off.Lzf,20,'b','fill')
plot(It,T40.off.Iz,'Color',[0.7 0.7 0.7],'LineWidth',1)
plot(It,T40.off.Ifz,'k','LineWidth',2)
plot(Ct,T40.off.Czf,'r','LineWidth',1)
plot(Lt,T40.off.Lzf,'b','LineWidth',1)
% grid on
box on
ylim([-.2 .45])
xlim(round([min(Ct) max(Ct)]))
h1=gca;
set(h1,'tickdir','in','xminortick','off','yminortick','off');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',20);
set(h1,'yticklabel',[]);
set(h1,'xtick',ynum,'xticklabel',{''});
title(ax3,'$\sigma_{\theta}=40^{\circ}$','interpreter','latex','fontsize',28);
text(ax3,ynum(1)+0.5, 0.41,'(c)','interpreter','latex','fontsize',24);
h2 = legend('press, DAC','press, NHR','Stereo','Lidar','interpreter','latex','fontsize',22)
% set(h2, 'Position', [0.68 0.86 0.10 0.05],'orientation','vertical')
% set(h2, 'Position', [0.8614 0.778 0.08 0.15],'orientation','vertical')
set(h2, 'Position', [0.862 0.778 0.08 0.15],'orientation','vertical')

ax4 = axes('Position',[0.54 0.09 0.42 0.35]);
plot([min(Ct) max(Ct)],[0 0],'LineWidth',1.5,'LineStyle','-.','Color',[0.7 0.7 0.7])
hold on
fill([-0.06 -0.06 0.06 0.06],[-0.3 0.5 0.5 -0.3],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
fill([-0.06 -0.06 0.06 0.06]-6.625,[-0.3 0.5 0.5 -0.3],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
scatter(It,T40.on.Iz,10,[0.7 0.7 0.7],'fill')
% scatter(It,T40.on.Ifz,4,'k','fill')
scatter(Ct,T40.on.Czf,20,'r','fill')
scatter(Lt,T40.on.Lzf,20,'b','fill')
plot(It,T40.on.Iz,'Color',[0.7 0.7 0.7],'LineWidth',1)
plot(It,T40.on.Ifz,'k','LineWidth',2)
plot(Ct,T40.on.Czf,'r','LineWidth',1)
plot(Lt,T40.on.Lzf,'b','LineWidth',1)
% grid on
box on
ylim([-.15 .2])
xlim(round([min(Ct) max(Ct)]))
h1=gca;
set(h1,'tickdir','in','xminortick','off','yminortick','off');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',20);
set(h1,'yticklabel',[]);
set(h1,'xtick',ynum,'xticklabel',ylab);
xlabel('Time (s)','interpreter','latex','fontsize',24);
text(ax4,ynum(1)+0.5, 0.17,'(d)','interpreter','latex','fontsize',24);

sname = [figfolder,'timeseries_filt'];
print(sname,'-dpng')
 