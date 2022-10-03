%% EOF analysis velocities
%  Trying to assess the dominant modes from the in situ gages.

%% Housekeeping

close all
clear all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

%% User Input

% define trial information
Tinfo.Hs        = 0.3; % wave height
Tinfo.Tp        = 2; % peak period
Tinfo.tide      = 1.07; % still water elevation
Tinfo.spread    = 40; % directional spread
time_range = [datenum(0,0,0,0,0,11364/8) datenum(0,0,0,0,0,12960/8)];

%% Standard trial/processing information

% functions to define trial structures
Tinfo       = trial_files(Tinfo); % read in trial information
Tinfo.cam   = TRC_camera_info(Tinfo.cam); % create cam matrix
camera.time	= cam_time(Tinfo);

% data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

%% Insitu

% load insitu data
F1          = matfile([Tinfo.datapath,'data/processed/insitu/',Tinfo.sz.tdate,'/',Tinfo.sz.tdate,'-insitu.mat']);
% press       = F1.press;
% wg          = F1.wg;
Tcam.sz_insitu    = F1.Tinfo;
vel       = F1.vel;
ptemp   	= F1.press;
Hz          = ptemp.Hz;
clear ptemp

t_inst = 0:datenum(0,0,0,0,0,1/Hz):((length(vel.u1)-1)*datenum(0,0,0,0,0,1/Hz));
% find overlapping range for insitu
[temp,istart] = nanmin(abs(time_range(1)-t_inst));
[temp,iend] = nanmin(abs(time_range(2)-t_inst));
inst.time = t_inst(istart:iend)';
istart = istart - Tinfo.offsets(1); 
iend = iend - Tinfo.offsets(1);
display('not sure if offset is right...')


velnames  = fieldnames(vel.xyzd);

% compute bulk statistics
for i = 1:length(velnames)
    % load pressure data
    instnum = str2double(velnames{i}(4:end));
    eval(['Xszh((instnum*2)-1,:)     = detrend(vel.u',velnames{i}(4:end),'(istart:iend));'])
    eval(['Xszh((instnum*2),:)     = detrend(vel.v',velnames{i}(4:end),'(istart:iend));'])
    eval(['Xsz((instnum*3)-2,:)     = detrend(vel.u',velnames{i}(4:end),'(istart:iend));'])
    eval(['Xsz((instnum*3)-1,:)     = detrend(vel.v',velnames{i}(4:end),'(istart:iend));'])
    eval(['Xsz((instnum*3),:)       = detrend(vel.w',velnames{i}(4:end),'(istart:iend));'])
    eval(['Xszf((instnum*3)-2,:)    = detrend(vel.u',velnames{i}(4:end),');'])
    eval(['Xszf((instnum*3)-1,:)    = detrend(vel.v',velnames{i}(4:end),');'])
    eval(['Xszf((instnum*3),:)      = detrend(vel.w',velnames{i}(4:end),');'])
    instname{instnum}             = velnames{i};
    eval(['xyz(instnum,:) = vel.xyzd.',velnames{i},'(1,1:3);']) 
end

% A = rand(size(Xsz,1),size(Xsz,2));
% A = wgn(size(Xsz,1),size(Xsz,2),-6);
% [R, t1, t2, t3, t4] = EOF(A);


[EOFs,lambdah,contribution,PCs] = EOF_analysis(Xszh);
[EOFs,lambdahs,contribution,PCs] = EOF_analysis(Xszh(1:16,:));
[EOFs,lambdahi,contribution,PCs] = EOF_analysis(Xszh(17:end,:));

[Y,sh] = calc_pcp(Xszh);
[Y,shs] = calc_pcp(Xszh(1:16,:)); 
[Y,shi] = calc_pcp(Xszh(17:end,:)); 
[Y,s] = calc_pcp(Xsz);
[Y,ss] = calc_pcp(Xsz(1:24,:)); 
[Y,si] = calc_pcp(Xsz(25:end,:)); 

%%

figure('units','inches','position',[1 1 6 5],'Color','w');
plot([1:length(sh.lambda)],sh.lambda,'k','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','k')
hold on
plot([1:length(shs.lambda)],shs.lambda,'r','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','r')
plot([1:length(shi.lambda)],shi.lambda,'b','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','b')
ylabel('Fraction of variance explained','interpreter','latex','fontsize',16);
h1=gca;
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1.1*get(h1,'ticklength'));
set(h1,'fontsize',16);
xlabel('Mode','interpreter','latex','fontsize',16);
title('Horizontal Velocities','interpreter','latex','fontsize',16);
ylim([0 0.5])
box on
h2 = legend('all gauges','outer surf zone','inner surf zone','interpreter','latex','fontsize',16)
set(h2,'orientation','vertical','Location','northeast');%'Position',[0.89 0.68 0.1 0.26])

    Sname1 = [Tinfo.figfolder,'sz_svd_uv'];
print(Sname1,'-dpng')

%%

figure('units','inches','position',[1 1 6 5],'Color','w');
plot([1:length(sh.lambda)],sh.lambda,'k','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','k')
hold on
plot([1:length(shs.lambda)],shs.lambda,'r','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','r')
plot([1:length(shi.lambda)],shi.lambda,'b','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','b')
plot([1:length(sh.lambda)],sh.lambda,'k','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','k')
plot([1:length(s.lambda)],s.lambda,'k','LineWidth',0.5,'Marker','x','Markersize',5,'MarkerFaceColor','k')
hold on
plot([1:length(ss.lambda)],ss.lambda,'r','LineWidth',0.5,'Marker','x','Markersize',5,'MarkerFaceColor','r')
plot([1:length(si.lambda)],si.lambda,'b','LineWidth',0.5,'Marker','x','Markersize',5,'MarkerFaceColor','b')
ylabel('Fraction of variance explained','interpreter','latex','fontsize',16);
h1=gca;
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1.1*get(h1,'ticklength'));
set(h1,'fontsize',16);
xlabel('Mode','interpreter','latex','fontsize',16);
title('3D Velocities','interpreter','latex','fontsize',16);
ylim([0 0.5])
box on
h2 = legend('all gauges','outer surf zone','inner surf zone','$u,v$','$u,v,w$','interpreter','latex','fontsize',16)
set(h2,'orientation','vertical','Location','northeast');%'Position',[0.89 0.68 0.1 0.26])

    Sname1 = [Tinfo.figfolder,'sz_svd_uvw'];
print(Sname1,'-dpng')

%%

figure('units','inches','position',[1 1 6 5],'Color','w');
plot([1:length(lambdah)],lambdah,'k','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','k')
hold on
plot([1:length(lambdahs)],lambdahs,'r','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','r')
plot([1:length(lambdahi)],lambdahi,'b','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','b')
plot([1:length(sh.lambda)],sh.lambda,'k','LineWidth',0.5,'Marker','x','Markersize',5,'MarkerFaceColor','k')
hold on
plot([1:length(shs.lambda)],shs.lambda,'r','LineWidth',0.5,'Marker','x','Markersize',5,'MarkerFaceColor','r')
plot([1:length(shi.lambda)],shi.lambda,'b','LineWidth',0.5,'Marker','x','Markersize',5,'MarkerFaceColor','b')
ylabel('Fraction of variance explained','interpreter','latex','fontsize',16);
h1=gca;
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1.1*get(h1,'ticklength'));
set(h1,'fontsize',16);
xlabel('Mode','interpreter','latex','fontsize',16);
title('Horizontal Velocities: Different code','interpreter','latex','fontsize',16);
ylim([0 0.5])
box on
h2 = legend('all gauges','outer surf zone','inner surf zone','CB code','other script','interpreter','latex','fontsize',16)
set(h2,'orientation','vertical','Location','northeast');%'Position',[0.89 0.68 0.1 0.26])

    Sname1 = [Tinfo.figfolder,'sz_svd_uv_difcode'];
print(Sname1,'-dpng')



%%

% load insitu data
F1          = matfile([Tinfo.datapath,'data/processed/insitu/',Tinfo.is.tdate,'/',Tinfo.is.tdate,'-insitu.mat']);

% press       = F1.press;
% wg          = F1.wg;
Tcam.is_insitu    = F1.Tinfo;
vel       = F1.vel;
ptemp   	= F1.press;
Hz          = ptemp.Hz;
clear ptemp

t_inst = 0:datenum(0,0,0,0,0,1/Hz):((length(vel.u1)-1)*datenum(0,0,0,0,0,1/Hz));
% find overlapping range for insitu
[temp,istart] = nanmin(abs(time_range(1)-t_inst));
[temp,iend] = nanmin(abs(time_range(2)-t_inst));
inst.time = t_inst(istart:iend)';
istart = istart - Tinfo.offsets(1);
iend = iend - Tinfo.offsets(1);


velnames  = fieldnames(vel.xyzd);

% compute bulk statistics
for i = 1:length(velnames)
    % load pressure data
    instnum = str2double(velnames{i}(4:end));
  	eval(['Xish((instnum*2)-1,:)     = vel.u',velnames{i}(4:end),'(istart:iend);'])
    eval(['Xish((instnum*2),:)     = vel.v',velnames{i}(4:end),'(istart:iend);'])
    eval(['Xis((instnum*3)-2,:)     = vel.u',velnames{i}(4:end),'(istart:iend);'])
    eval(['Xis((instnum*3)-1,:)     = vel.v',velnames{i}(4:end),'(istart:iend);'])
    eval(['Xis((instnum*3),:)       = vel.w',velnames{i}(4:end),'(istart:iend);'])
    eval(['Xisf((instnum*3)-2,:)    = vel.u',velnames{i}(4:end),';'])
    eval(['Xisf((instnum*3)-1,:)    = vel.v',velnames{i}(4:end),';'])
    eval(['Xisf((instnum*3),:)      = vel.w',velnames{i}(4:end),';'])
    instname{instnum}             = velnames{i};
    eval(['xyz(instnum,:) = vel.xyzd.',velnames{i},'(1,1:3);']) 
    
end

[EOFs,lambdah,contribution,PCs] = EOF_analysis(Xish);
[EOFs,lambdahs,contribution,PCs] = EOF_analysis(Xish(7:end,:));
[EOFs,lambdahi,contribution,PCs] = EOF_analysis(Xish(1:6,:));

[Y,sh] = calc_pcp(Xish);
[Y,shs] = calc_pcp(Xish(7:end,:)); 
[Y,shi] = calc_pcp(Xish(1:6,:)); 
[Y,s] = calc_pcp(Xis);
[Y,ss] = calc_pcp(Xis(10:end,:)); 
[Y,si] = calc_pcp(Xis(1:9,:)); 

%%

figure('units','inches','position',[1 1 6 5],'Color','w');
plot([1:length(sh.lambda)],sh.lambda,'k','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','k')
hold on
plot([1:length(shs.lambda)],shs.lambda,'r','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','r')
plot([1:length(shi.lambda)],shi.lambda,'b','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','b')
ylabel('Fraction of variance explained','interpreter','latex','fontsize',16);
h1=gca;
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1.1*get(h1,'ticklength'));
set(h1,'fontsize',16);
xlabel('Mode','interpreter','latex','fontsize',16);
title('Horizontal Velocities','interpreter','latex','fontsize',16);
ylim([0 0.5])
box on
h2 = legend('all gauges','surfzone edge','inner shelf','interpreter','latex','fontsize',16)
set(h2,'orientation','vertical','Location','northeast');%'Position',[0.89 0.68 0.1 0.26])

    Sname1 = [Tinfo.figfolder,'is_svd_uv'];
print(Sname1,'-dpng')

%%

figure('units','inches','position',[1 1 6 5],'Color','w');
plot([1:length(sh.lambda)],sh.lambda,'k','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','k')
hold on
plot([1:length(shs.lambda)],shs.lambda,'r','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','r')
plot([1:length(shi.lambda)],shi.lambda,'b','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','b')
plot([1:length(sh.lambda)],sh.lambda,'k','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','k')
plot([1:length(s.lambda)],s.lambda,'k','LineWidth',0.5,'Marker','x','Markersize',5,'MarkerFaceColor','k')
hold on
plot([1:length(ss.lambda)],ss.lambda,'r','LineWidth',0.5,'Marker','x','Markersize',5,'MarkerFaceColor','r')
plot([1:length(si.lambda)],si.lambda,'b','LineWidth',0.5,'Marker','x','Markersize',5,'MarkerFaceColor','b')
ylabel('Fraction of variance explained','interpreter','latex','fontsize',16);
h1=gca;
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1.1*get(h1,'ticklength'));
set(h1,'fontsize',16);
xlabel('Mode','interpreter','latex','fontsize',16);
title('3D Velocities','interpreter','latex','fontsize',16);
ylim([0 0.5])
box on
h2 = legend('all gauges','surfzone edge','inner shelf','$u,v$','$u,v,w$','interpreter','latex','fontsize',16)
set(h2,'orientation','vertical','Location','northeast');%'Position',[0.89 0.68 0.1 0.26])

    Sname1 = [Tinfo.figfolder,'is_svd_uvw'];
print(Sname1,'-dpng')

%%

figure('units','inches','position',[1 1 6 5],'Color','w');
plot([1:length(lambdah)],lambdah,'k','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','k')
hold on
plot([1:length(lambdahs)],lambdahs,'r','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','r')
plot([1:length(lambdahi)],lambdahi,'b','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','b')
plot([1:length(sh.lambda)],sh.lambda,'k','LineWidth',0.5,'Marker','x','Markersize',5,'MarkerFaceColor','k')
hold on
plot([1:length(shs.lambda)],shs.lambda,'r','LineWidth',0.5,'Marker','x','Markersize',5,'MarkerFaceColor','r')
plot([1:length(shi.lambda)],shi.lambda,'b','LineWidth',0.5,'Marker','x','Markersize',5,'MarkerFaceColor','b')
ylabel('Fraction of variance explained','interpreter','latex','fontsize',16);
h1=gca;
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1.1*get(h1,'ticklength'));
set(h1,'fontsize',16);
xlabel('Mode','interpreter','latex','fontsize',16);
title('Horizontal Velocities: Different code','interpreter','latex','fontsize',16);
ylim([0 0.5])
box on
h2 = legend('all gauges','surfzone edge','inner shelf','CB code','other script','interpreter','latex','fontsize',16)
set(h2,'orientation','vertical','Location','northeast');%'Position',[0.89 0.68 0.1 0.26])

    Sname1 = [Tinfo.figfolder,'is_svd_uv_difcode'];
print(Sname1,'-dpng')



% %%
% [Lh, EOFs, EC, error, norms] = EOF1(Xszh);
% [L, EOFs, EC, error, norms] = EOF1(Xsz);
% [Li, EOFs, EC, error, norms] = EOF1(Xsz(1:24,:));
% [Ls, EOFs, EC, error, norms] = EOF1(Xsz(25:end,:));
% [Lhi, EOFs, EC, error, norms] = EOF1(Xszh(1:16,:));
% [Lhs, EOFs, EC, error, norms] = EOF1(Xszh(17:end,:));
% 
% figure('units','inches','position',[1 1 6 5],'Color','w');
% plot([1:length(L)],L,'k','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','k')
% hold on
% plot([1:length(Ls)],Ls,'r','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','r')
% plot([1:length(Li)],Li,'b','LineWidth',0.5,'Marker','o','Markersize',5,'MarkerFaceColor','b')
% 
% ylabel('$\sigma_{\theta,\mathrm{c}} (^{\circ})$','interpreter','latex','fontsize',24);
% h1=gca;
% set(h1,'tickdir','in','xminortick','on','yminortick','on');
% set(h1,'ticklength',1.1*get(h1,'ticklength'));
% set(h1,'fontsize',22);
% xlabel('$\sigma_{\theta}$ at wavemaker $(^{\circ})$','interpreter','latex','fontsize',22);
% ylim([8 30])
% box on
% % xlim([-1 41])
% % text(-0.2,2.2,'(a)','interpreter','latex','fontsize',24);
% 
% set(h1,'xtick',[0:10:40],'xticklabel',{'0', '10', '20', '30', '40'});
% h2 = legend('Stereo','Imagery','$H_s = 0.25$ m','$H_s = 0.30$ m','interpreter','latex','fontsize',20)
% set(h2,'orientation','vertical','Location','northwest');%'Position',[0.89 0.68 0.1 0.26])
% 
%     Sname1 = [figfolder,'angle_std'];
% print(Sname1,'-dpng')