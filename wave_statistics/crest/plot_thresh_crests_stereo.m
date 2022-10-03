% plot ylen

% Plot XYZ coordinates of point cloud exported from Photoscan.
% This will make a series of several plots of sea surface elevation and
% compute the sea surface elevation as a function of time. 

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\cameras'))
addpath(genpath('E:\codes\trc_lab_experiment\toolbox'))
% addpath('/Users/cmbaker9/Documents/MTOOLS/nctoolbox/cdm/utilities/graphics')

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% wave conditions: 
Tinfo.Hs = 0.25;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
sprd = [10 20 30 40];
Tinfo.filt = 1;

for isprd = 1:4
    Tinfo.spread = sprd(isprd);
dx =0.25;% for binning
samprate = 8; %Hz
Wtank = 26.25;%26.5; % width of the tank

%%%%%%%%%%%%%%%%%% END OF USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%

%% Prepare file structures, folders, etc.

Tinfo = trial_files(Tinfo);

% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

%% STEP 2: Define figure folders

figfolder = [Tinfo.figfolder,'wave_crests\'];
eval(['!mkdir ',figfolder])

%% Load stereo reconstruction data
gamma = 0.45;% image threshold
minarea = 200; % min area of object
xreg=[26.5 31.5];
subname = ['_x',num2str(round(xreg(1))),'to',num2str(round(xreg(2))),'_waveedge'];
% stereo_thresh = load([Tinfo.savefolder,'ylen_stereo_gamma',num2str(gamma*10),'_minarea',num2str(minarea),'_samprate',num2str(samprate),subname,'.mat']);
stereo_thresh = load([Tinfo.savefolder,'ylen_stereo_gamma',num2str(gamma*10),'_minarea',num2str(minarea),'_samprate',num2str(samprate),subname,'_selectedregions.mat']);

%% Load stereo

if Tinfo.filt == 1
    F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(-14)),'to',num2str(14),'m_resx',num2str((Tinfo.cam.dx)*100),'cm_resy',num2str((Tinfo.cam.dx)*100),'cm_filtered.mat'];
elseif Tinfo.filt == 0
    F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(14)),'to',num2str(14),'m_resx',num2str((Tinfo.cam.dx)*100),'cm_resy',num2str((Tinfo.cam.dx)*100),'cm.mat'];
end
    stereo = load(F1,'x','y','z');


%% 
time = 0:0.125:(60*10);

v = VideoWriter([figfolder,'wavecrests_threshold.avi']);
v.FrameRate=8; %?????
open(v)

% start plotting
figure('units','inches','position',[1 1 5 12],'color','w'); % open figure

for i = 1:8*60
%     I = double(rgb2gray(imread(fullfile(odir,['c2_',L(i,:),'.tiff']))));
%     I = imread(fullfile(odir,['c2_',L(i,:),'.tiff']));
%     
%     ax1 = axes('Position',[0.07 0.095 0.4 0.8]);
%     imagesc(ax1,image.X(1,:),image.Y(:,1),I);
%     axis image
%     axis equal
%     ylim([-13,13]);
%     xlim(ixlim);
%     set(ax1,'YDir','normal')
%     
%     ax2 = axes('Position',[0.07 0.095 0.4 0.8]);
%     Ithresh = double(image_thresh.Ithresh(:,:,i));
% % %     Ibound  = bwboundaries(image_thresh.Ithresh(:,:,i));
% %     Ibound = activecontour(Ithresh,ones(size(Ithresh)),300);
%     Ithresh(Ithresh==0)=NaN;
% %     Ithresh(124:180,25:56)=NaN;
%     h2 = pcolor(image_thresh.X,image_thresh.Y,Ithresh)
%     set(h2,'EdgeColor','none')
%     alpha(h2,0.5)
%     hold on
% %     h2b = visboundaries(Ibound,'Color','r');
% %     h2b.XData = image_thresh.X;
% %     h2b.YData = image_thresh.Y;
% %     alpha(h2b,1)
% %     plot([27.7 27.7],[-15 15],'k','LineStyle','--','LineWidth',2)
%     axis equal
%     xlim(ixlim)
%     ylim([-13 13])
%     
%     linkaxes([ax1,ax2])
%     %Hide the top axes
%     ax2.Visible = 'off';
%     ax2.XTick = [];
%     ax2.YTick = [];
%     %Give each one its own colormap
%     colormap(ax2,'hot')
%     ylabel(ax1,'$y$~(m)','interpreter','latex','fontsize',20);
%     xlabel(ax1,'$x$~(m)','interpreter','latex','fontsize',20);
%     h1=ax1;
%     set(h1,'tickdir','out','xminortick','off','yminortick','off');
%     set(h1,'ticklength',1*get(h1,'ticklength'));
%     set(h1,'fontsize',16);
    
    
    ax3 = axes('Position',[0.06 0.095 0.8 0.8]);
    h3 = pcolor(ax3,stereo.x,stereo.y,flipud(stereo.z(:,:,i))-Tinfo.tide-.03);
    set(h3,'EdgeColor','none')
    axis equal
    xlim([25 35])
    ylim([-13 13])
%     axis equal
%     hold on
    ax4 = axes('Position',[0.06 0.095 0.8 0.8]);
    Sthresh = flipud(double(stereo_thresh.Ithresh(:,:,i)));
    Sthresh(Sthresh==0)=NaN;
    h4 = pcolor(stereo_thresh.X,stereo_thresh.Y,Sthresh)
    set(h4,'EdgeColor','none')
    alpha(h4,0.5)
    hold on
%     plot([27.7 27.7],[-15 15],'k','LineStyle','--','LineWidth',2)
    axis equal
    xlim([25 35])
    ylim([-13 13])
%     axis equal
    % scatter([press.xyz(2,9) press.xyz(2,3)],[press.xyz(1,9) press.xyz(1,3)],40,'r','fill')
    % text(ax2,press.xyz(2,9)-0.35,press.xyz(1,9)-0.6,'p06','interpreter','latex','fontsize',15)
    % text(ax2,press.xyz(2,3)-0.35,press.xyz(1,3)-0.6,'p11','interpreter','latex','fontsize',15)
    shading flat
    
    linkaxes([ax3,ax4])
    %Hide the top axes
    ax4.Visible = 'off';
    ax4.XTick = [];
    ax4.YTick = [];
    %Give each one its own colormap
    colormap(ax4,'hot')
    %Then add colorbars and get everything lined up
%     set([ax3,ax4],'Position',[.17 .11 .685 .815]);
    
    colormap(ax3,flipud(cmocean('deep')));
    hc = colorbar(ax3,'Position', [0.76 0.095 0.05 0.74]);
    caxis(ax3,[-0.15001 0.15001])
%     axis equal
    %     grid on
    box on
%     xlim(ixlim)
%     ylim([-13 13])
%     axis equal
    h1=ax3;
%     set(h1,'ydir','reverse')
    set(h1,'tickdir','out','xminortick','off','yminortick','off');
    set(h1,'ticklength',1*get(h1,'ticklength'));
    set(h1,'fontsize',16);
    set(h1,'ytick',[-10:5:10],'yticklabel',{''});
    text(ax3,35.6,12.3,'$\eta$ (m)','interpreter','latex','fontsize',20);
    ylabel('$y$ (m)','interpreter','latex','fontsize',20);
    xlabel(ax3,'$x$ (m)','interpreter','latex','fontsize',20);
%     text(ax2,25.7,14.5,'(a)','interpreter','latex','fontsize',20);
%     text(ax2,27.5,15,'Stereo','interpreter','latex','fontsize',20);
%     text(ax3,25,14,'(b) Stereo Reconstruction','interpreter','latex','fontsize',20);
    set(h1,'ytick',[-10:5:10],'yticklabel',{'-10' '-5' '0' '5' '10'});
%     text(ax3,25,14,'(a) Image','interpreter','latex','fontsize',20);
    text(ax3,25,15.5,['Time: ',num2str(time(i)),' s'],'interpreter','latex','fontsize',20);

    
%     Sname = [figfolder,'wavecrest_threshold_'];
%     print(Sname,'-dpng')

    writeVideo(v,getframe(gcf))  
   
    
%     pause(2)
    clf
%     count = count+icount;
end

close(v) 
end
% disp('saved movie')