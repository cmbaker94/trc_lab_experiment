% function threshold_combined_crest_length(Tinfo,gammab,IMthresh,minarea,samprate,pltflag,ixlim,iylim,idxdy,ixsel)% This code will plot rectified images using code from
% D_gridGenExpampleRect.m in the CIRN-Quantitative-Coastal-Imaging-Toolbox

close all
clear all
clc

% Set up paths and clear workspace
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\cameras'))
addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\CIRN-Quantitative-Coastal-Imaging-Toolbox/X_CoreFunctions/'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

idxdy=0.05;
ixlim=[25 33.5];%25
iylim=[-14 14];
imlen=0;
% 
% % % threshold images
% gammab = 0.2;
% IMthresh = 0.35;%0.4   or 0.3 for hs=0.25 0.7%0.6;%0.6 % image threshold
minareaIMST = 40;%150     30;%800 % min area of object
minarea = 100;
samprate =8; %Hz
pltflag = 1;
ixsel = [27 31.5];
% 
% sprd = [0, 20, 30, 40, 0, 10, 20, 30, 40]%, 40, 40, 20, 40];%[30,  40];%
% Hs = [0.3, 0.3, 0.3, 0.3, 0.25, 0.25, 0.25, 0.25, 0.25]%, 0.2, 0.3, 0.3, 0.3]; %[0.25, 0.25];%
% Tp = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3];
% h = [1.07, 1.07, 1.07, 1.07, 1.07, 1.07, 1.07, 1.07, 1.07, 1.07, 1.00, 1.07, 1.07];


sprd = [20, 40];%[30,  40];%
Hs = [0.3, 0.3]; %[0.25, 0.25];%
Tp = [3, 3];
h = [1.07, 1.07];

for op = 1:length(sprd)
    Tinfo.Hs = Hs(op);
    Tinfo.Tp = Tp(op);
    Tinfo.tide = h(op);
    Tinfo.spread = sprd(op);
    Tinfo = trial_files(Tinfo);

% filter data
Tinfo.filt = 1; % 1 if true

%% STEP 1: Create paths, files and naming
% general path and names
datapath    = 'E:/';

Tinfo = trial_files(Tinfo);

% Stereo Reconstructions
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);


%% Load stereo

subname = '';
Tinfo.cam.regy = iylim;
% F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_idxdy',num2str((idxdy)*100),'cm_resy',num2str((idxdy)*100),'cm',subname,'.mat'];
if Tinfo.filt == 1
    F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_resx',num2str((Tinfo.cam.dx)*100),'cm_resy',num2str((Tinfo.cam.dx)*100),'cm_filtered.mat'];
elseif Tinfo.filt == 0
    F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_resx',num2str((Tinfo.cam.dx)*100),'cm_resy',num2str((Tinfo.cam.dx)*100),'cm.mat'];
end


stereo = load(F1,'x','y','z');
mssetemp = nanmean(stereo.z,3);
msse(:,:,op) = mssetemp;
display(op)
end

%%


% xtemp = stereo.x(1,:);
% ytemp = stereo.y(:,1);
% 
% figure('units','inches','position',[1 1 16 10],'color','w')
% 
% subplot(251)
% imagesc(xtemp,ytemp,msse(:,:,5))
% ylim([-14 14])
% xlim([27 33])
% % colormap(cmocean('hot'))
% caxis([1.02 1.13])
% h1=gca;
% set(h1,'fontsize',14);
% set(h1,'tickdir','in');%,'xminortick','on','yminortick','on');
% set(h1,'ticklength',1*get(h1,'ticklength'));
% title('$\sigma_{\theta}=0^{\circ}$','interpreter','latex','fontsize',20);
% ylabel('$H_s=$ 0.25 m','interpreter','latex','fontsize',24);
% 
% subplot(252)
% imagesc(xtemp,ytemp,msse(:,:,6))
% ylim([-14 14])
% xlim([27 33])
% % colormap('gray')
% caxis([1.02 1.13])
% h1=gca;
% set(h1,'fontsize',14);
% set(h1,'tickdir','in');%,'xminortick','on','yminortick','on');
% set(h1,'ticklength',1*get(h1,'ticklength'));
% title('$\sigma_{\theta}=10^{\circ}$','interpreter','latex','fontsize',20);
% 
% subplot(253)
% imagesc(xtemp,ytemp,msse(:,:,7))
% ylim([-14 14])
% xlim([27 33])
% % colormap('gray')
% caxis([1.02 1.13])
% h1=gca;
% set(h1,'fontsize',14);
% set(h1,'tickdir','in');%,'xminortick','on','yminortick','on');
% set(h1,'ticklength',1*get(h1,'ticklength'));
% title('$\sigma_{\theta}=20^{\circ}$','interpreter','latex','fontsize',20);
% 
% subplot(254)
% imagesc(xtemp,ytemp,msse(:,:,8))
% ylim([-14 14])
% xlim([27 33])
% % colormap('gray')
% caxis([1.02 1.13])
% h1=gca;
% set(h1,'fontsize',14);
% set(h1,'tickdir','in');%,'xminortick','on','yminortick','on');
% set(h1,'ticklength',1*get(h1,'ticklength'));
% title('$\sigma_{\theta}=30^{\circ}$','interpreter','latex','fontsize',20);
% 
% subplot(255)
% imagesc(xtemp,ytemp,msse(:,:,9))
% ylim([-14 14])
% xlim([27 33])
% % colormap('gray')
% caxis([1.02 1.13])
% h1=gca;
% set(h1,'fontsize',14);
% set(h1,'tickdir','in');%,'xminortick','on','yminortick','on');
% set(h1,'ticklength',1*get(h1,'ticklength'));
% title('$\sigma_{\theta}=40^{\circ}$','interpreter','latex','fontsize',20);
% 
% subplot(256)
% imagesc(xtemp,ytemp,msse(:,:,1))
% ylim([-14 14])
% xlim([27 33])
% % colormap('gray')
% caxis([1.02 1.13])
% h1=gca;
% set(h1,'fontsize',14);
% set(h1,'tickdir','in');%,'xminortick','on','yminortick','on');
% set(h1,'ticklength',1*get(h1,'ticklength'));
% ylabel('$H_s=$ 0.30 m','interpreter','latex','fontsize',24);
% 
% subplot(258)
% imagesc(xtemp,ytemp,msse(:,:,2))
% ylim([-14 14])
% xlim([27 33])
% % colormap('gray')
% caxis([1.02 1.13])
% h1=gca;
% set(h1,'fontsize',14);
% set(h1,'tickdir','in');%,'xminortick','on','yminortick','on');
% set(h1,'ticklength',1*get(h1,'ticklength'));
% 
% subplot(259)
% imagesc(xtemp,ytemp,msse(:,:,3))
% ylim([-14 14])
% xlim([27 33])
% % colormap('gray')
% caxis([1.02 1.13])
% h1=gca;
% set(h1,'fontsize',14);
% set(h1,'tickdir','in');%,'xminortick','on','yminortick','on');
% set(h1,'ticklength',1*get(h1,'ticklength'));
% 
% subplot(2,5,10)
% imagesc(xtemp,ytemp,msse(:,:,4))
% ylim([-14 14])
% xlim([27 33])
% % colormap('gray')
% caxis([1.02 1.13])
% h1=gca;
% set(h1,'fontsize',14);
% set(h1,'tickdir','in');%,'xminortick','on','yminortick','on');
% set(h1,'ticklength',1*get(h1,'ticklength'));
% hc = colorbar('Position', [0.93 0.12 0.015 0.8]);
% 
% Sname = ['E:\figures\methods_manuscript\testing\','stereo_msse'];
% print(Sname,'-dpng')

%%

xtemp = stereo.x(1,:);
ytemp = stereo.y(:,1);

figure('units','inches','position',[1 1 16 10],'color','w')

subplot(121)
imagesc(xtemp,ytemp,msse(:,:,1))
ylim([-14 14])
xlim([27 33])
% colormap(cmocean('hot'))
caxis([1.02 1.13])
h1=gca;
set(h1,'fontsize',14);
set(h1,'tickdir','in');%,'xminortick','on','yminortick','on');
set(h1,'ticklength',1*get(h1,'ticklength'));
title('$\sigma_{\theta}=20^{\circ}$','interpreter','latex','fontsize',20);
ylabel('$H_s$ = 0.3 m, $T_p=$ 3 s','interpreter','latex','fontsize',24);

subplot(122)
imagesc(xtemp,ytemp,msse(:,:,2))
ylim([-14 14])
xlim([27 33])
% colormap('gray')
caxis([1.02 1.13])
h1=gca;
set(h1,'fontsize',14);
set(h1,'tickdir','in');%,'xminortick','on','yminortick','on');
set(h1,'ticklength',1*get(h1,'ticklength'));
title('$\sigma_{\theta}=40^{\circ}$','interpreter','latex','fontsize',20);

Sname = ['E:\figures\methods_manuscript\testing\','stereo_msse_tp3'];
print(Sname,'-dpng')

