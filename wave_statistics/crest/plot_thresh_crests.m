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
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))
% addpath('/Users/cmbaker9/Documents/MTOOLS/nctoolbox/cdm/utilities/graphics')

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% wave conditions: 
Tinfo.Hs = 0.25;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
sprd = 40;%[0 10 20 30 40];
Tinfo.filt = 1;

for isprd = 1:length(sprd)
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

thresh = 0.35; % image threshold
minarea = 200; % min area of object
xreg = [26.5 31.5];
subname = ['_x',num2str(round(xreg(1))),'to',num2str(round(xreg(2))),'_ycrest'];
% image_thresh = load([Tinfo.savefolder,'ylen_thresh',num2str(thresh*100),'_minarea',num2str(minarea),'_samprate',num2str(samprate),subname,'.mat']);
image_thresh = load([Tinfo.savefolder,'ylen_thresh',num2str(thresh*100),'_minarea',num2str(minarea),'_samprate',num2str(samprate),subname,'_selectedregions.mat']);

gamma = 0.35;% image threshold
minarea = 200; % min area of object
xreg=[26.5 31.5];
idxdy = 0.02;
subname = ['_x',num2str(round(xreg(1))),'to',num2str(round(xreg(2))),'_ycrest_res',num2str(idxdy)];
% stereo_thresh = load([Tinfo.savefolder,'ylen_stereo_gamma',num2str(gamma*10),'_minarea',num2str(minarea),'_samprate',num2str(samprate),subname,'.mat']);
stereo_thresh = load([Tinfo.savefolder,'ylen_stereo_gamma',num2str(gamma*10),'_minarea',num2str(minarea),'_samprate',num2str(samprate),subname,'_selectedregions.mat']);

%% Load stereo

% if Tinfo.filt == 1
%     F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(-14)),'to',num2str(14),'m_resx',num2str((Tinfo.cam.dx)*100),'cm_resy',num2str((Tinfo.cam.dx)*100),'cm_filtered.mat'];
% elseif Tinfo.filt == 0
%     F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(14)),'to',num2str(14),'m_resx',num2str((Tinfo.cam.dx)*100),'cm_resy',num2str((Tinfo.cam.dx)*100),'cm.mat'];
% end
%     stereo = load(F1,'x','y','z');

tiffpath = ['F:\Metashape\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1\',Tinfo.cam.trimname,'dems\'];
% imagepath = ['H:\TRC_Fall_Experiment\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1_JPEG\'];
% tiffpath = ['D:\TRC_Fall_Experiment\photoscan\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1\',Tinfo.cam.trimname,'\dems\'];

%% Load insitu data

[pg,~,~] = load_insitu(Tinfo);

%% Prep images

idxdy=0.02;
ixlim=[25 35];%[18 35];
iylim=[-14 14];
geoname = ['x',num2str(ixlim(1)),'to',num2str(round(ixlim(2))),'_y',num2str(round(iylim(2))),'_res',num2str(idxdy*100),'cm'];
odir = [Tinfo.savefolder(1:74),'orthos\',geoname];
cnames = ls([odir,'\c2_*.tiff']);
L = cnames(:,4:8);
xtemp = ixlim(1):idxdy:ixlim(2);
ytemp = iylim(1):idxdy:iylim(2);
image.Y = repmat(ytemp',1,length(xtemp));
image.X = repmat(xtemp,length(ytemp),1);

%% 
time = 0:0.125:(60*10);
resy = 0.05;

v = VideoWriter([figfolder,'wavecrests_threshold.avi']);
v.FrameRate=8; %?????
open(v)

% start plotting
figure('units','inches','position',[1 1 9 10],'color','w'); % open figure

for i = 1:size(L,1)
    tifffile    = [tiffpath,Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_frames_',Tinfo.cam.imagerange,'_DEM_',num2str(i),'.tif'];
    % read in point cloud
    [xtemp,ytemp,ztemp,resxtemp,resytemp] = read_tiff(tifffile);
    ytemp = fliplr(ytemp);
    ztemp = flipud(ztemp);
    
    % need to make a 2d for finding nans
    ytemp = repmat(ytemp',1,length(xtemp));
    xtemp = repmat(xtemp,length(ytemp),1);
    
%     I = double(rgb2gray(imread(fullfile(odir,['c2_',L(i,:),'.tiff']))));
    I = imread(fullfile(odir,['c2_',L(i,:),'.tiff']));
    
    ax1 = axes('Position',[0.07 0.095 0.4 0.8]);
    imagesc(ax1,image.X(1,:),image.Y(:,1),I);
    axis image
    axis equal
    ylim([-13,13]);
    xlim(ixlim);
    set(ax1,'YDir','normal')
    
    ax2 = axes('Position',[0.07 0.095 0.4 0.8]);
    Ithresh = double(image_thresh.Ithresh(:,:,i));
     Ithresh(Ithresh==0)=NaN;
    h2 = pcolor(image_thresh.X,image_thresh.Y,Ithresh)
    set(h2,'EdgeColor','none')
    alpha(h2,0.5)
    hold on
    
    [Ln,n] = bwlabel(image_thresh.Ithresh(:,:,i));
    for ni = 1:n
        obtemp = zeros(size(Ithresh));
        obtemp(Ln==ni)=1;
        count = 0;
        for ro = 1 : size(obtemp, 1)
            if sum(obtemp(ro,:))>0
                count = count+1;
                B(1) = (find(obtemp(ro, :), 1, 'first')*resy)+xreg(1);
                B(2) = (find(obtemp(ro, :), 1, 'last')*resy)+xreg(1);
                C(count, 1) = mean(B);
                C(count, 2) = (ro*resy)+iylim(1);
            end
            clear B
        end
        plot(C(:,1),C(:,2),'r','LineWidth',2);%,'Marker','o','MarkerSize',2)
        clear C
    end
    props = regionprops(image_thresh.Ithresh(:,:,i), 'MajorAxisLength', 'Orientation', 'Centroid', 'MinorAxisLength');
    t = linspace(0,2*pi,50);
        for ni = 1:length(props)
            a = (props(ni).MajorAxisLength*resy)/2;
            b = (props(ni).MinorAxisLength*resy)/2;
            Xc = props(ni).Centroid(1)*resy+xreg(1);
            Yc = props(ni).Centroid(2)*resy+iylim(1);
            phi = deg2rad(-props(ni).Orientation);
            xe = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
            ye = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
            hc = plot(xe,ye,'y','Linewidth',2)
            hc.Color = [hc.Color 0.75];
        end    
    axis equal
    xlim(ixlim)
    ylim([-13 13])
    
    linkaxes([ax1,ax2])
    %Hide the top axes
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    %Give each one its own colormap
    colormap(ax2,'hot')
    ylabel(ax1,'$y$~(m)','interpreter','latex','fontsize',20);
    xlabel(ax1,'$x$~(m)','interpreter','latex','fontsize',20);
    h1=ax1;
    set(h1,'tickdir','out','xminortick','off','yminortick','off');
    set(h1,'ticklength',1*get(h1,'ticklength'));
    set(h1,'fontsize',16);
    set(h1,'ytick',[-10:5:10],'yticklabel',{'-10' '-5' '0' '5' '10'})
    text(ax1,25,14,'(a) Image','interpreter','latex','fontsize',20);
    text(ax1,25,15.5,['Time: ',num2str(time(i)),' s'],'interpreter','latex','fontsize',20);

    
    ax3 = axes('Position',[0.46 0.095 0.4 0.8]);
    h3 = pcolor(ax3,xtemp,ytemp,flipud(ztemp)-Tinfo.tide-.03);
%     h3 = pcolor(ax3,stereo.x,stereo.y,flipud(stereo.z(:,:,i))-Tinfo.tide-.03)
    set(h3,'EdgeColor','none')
    axis equal
    xlim(ixlim)
    ylim([-13 13])
%     axis equal
%     hold on
    ax4 = axes('Position',[0.46 0.095 0.4 0.8]);
    Sthresh = flipud(double(stereo_thresh.Ithresh(:,:,i)));
    Sthresh(Sthresh==0)=NaN;
    h4 = pcolor(stereo_thresh.X,stereo_thresh.Y,Sthresh)
    set(h4,'EdgeColor','none')
    alpha(h4,0.5)
    hold on
    [Ln,n] = bwlabel(flipud(stereo_thresh.Ithresh(:,:,i)));
    for ni = 1:n
        obtemp = zeros(size(Sthresh));
        obtemp(Ln==ni)=1;
        count = 0;
        for ro = 1 : size(obtemp, 1)
            if sum(obtemp(ro,:))>0
                count = count+1;
                B(1) = (find(obtemp(ro, :), 1, 'first')*resy)+xreg(1);
                B(2) = (find(obtemp(ro, :), 1, 'last')*resy)+xreg(1);
                C(count, 1) = mean(B);
                C(count, 2) = (ro*resy)+min(stereo_thresh.Y,[],'all');
            end
            clear B
        end
        plot(C(:,1),C(:,2),'r','LineWidth',2);%,'Marker','o','MarkerSize',2)
        clear C
    end
    props = regionprops(flipud(stereo_thresh.Ithresh(:,:,i)), 'MajorAxisLength', 'Orientation', 'Centroid', 'MinorAxisLength');
    t = linspace(0,2*pi,50);
        for ni = 1:length(props)
            a = (props(ni).MajorAxisLength*resy)/2;
            b = (props(ni).MinorAxisLength*resy)/2;
            Xc = props(ni).Centroid(1)*resy+xreg(1);
            Yc = props(ni).Centroid(2)*resy+min(stereo_thresh.Y,[],'all');
            phi = deg2rad(-props(ni).Orientation);
            xe = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
            ye = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
            hc = plot(xe,ye,'y','Linewidth',2)
            hc.Color = [hc.Color 0.75];
        end   
    
    
%     plot([27.7 27.7],[-15 15],'k','LineStyle','--','LineWidth',2)
    axis equal
    xlim(ixlim)
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
    hc = colorbar(ax3,'Position', [0.86 0.095 0.03 0.74]);
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
%     ylabel('$y$ (m)','interpreter','latex','fontsize',20);
    xlabel(ax3,'$x$ (m)','interpreter','latex','fontsize',20);
%     text(ax2,25.7,14.5,'(a)','interpreter','latex','fontsize',20);
%     text(ax2,27.5,15,'Stereo','interpreter','latex','fontsize',20);
    text(ax3,25,14,'(b) Stereo Reconstruction','interpreter','latex','fontsize',20);
    
    
    Sname = [figfolder,'wavecrest_threshold_',L(i,:)];
    print(Sname,'-dpng')

    writeVideo(v,getframe(gcf))  
   
    
%     pause(2)
    clf
%     count = count+icount;
end

close(v) 

%%
time = 0:0.125:(60*10);

v = VideoWriter([figfolder,'wavecrests_threshold_flipped.avi']);
v.FrameRate=8; %?????
open(v)

%
% start plotting
figure('units','inches','position',[1 1 12 9],'color','w'); % open figure

for i = 1:8*60*1.5%size(L,1)
    tifffile    = [tiffpath,Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_frames_',Tinfo.cam.imagerange,'_DEM_',num2str(i),'.tif'];
    % read in point cloud
    [xtemp,ytemp,ztemp,resxtemp,resytemp] = read_tiff(tifffile);
    ytemp = fliplr(ytemp);
    ztemp = flipud(ztemp);
    
    % need to make a 2d for finding nans
    ytemp = repmat(ytemp',1,length(xtemp));
    xtemp = repmat(xtemp,length(ytemp),1);
    
    I = imread(fullfile(odir,['c2_',L(i,:),'.tiff']));
    
    ax1 = axes('Position',[0.095 0.5 0.8 0.4]);
    imagesc(ax1,image.Y(:,1),image.X(1,:),imrotate(I,270));
    hold on
    scatter(pg.xyz(2,:),pg.xyz(1,:),50,'r','fill','MarkerEdgeColor','k','LineWidth',1)
    plot([-13 13],[27.6 27.6],'r','LineStyle','--','LineWidth',2)
    axis image
    axis equal
    xlim([-13,13]);
    ylim([26 35]);
    set(ax1,'YDir','reverse')
    hold on
    ax2 = axes('Position',[0.095 0.5 0.8 0.4]);
    Ithresh = double(image_thresh.Ithresh(:,:,i));
    Ithresh(Ithresh==0)=NaN;
%     h2 = pcolor(image_thresh.Y',image_thresh.X',fliplr(Ithresh'))
%     set(h2,'EdgeColor','none')
%     alpha(h2,0.5)
    hold on
    axis equal
    ylim([26 35])
    xlim([-13 13])
    
    linkaxes([ax1,ax2])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    colormap(ax2,'hot')
    h1=ax1;
    set(h1,'fontsize',16);
    ylabel(ax1,'Cross-shore~(m)','interpreter','latex','fontsize',19);
    set(h1,'tickdir','out','xminortick','off','yminortick','off');
    set(h1,'ticklength',0.6*get(h1,'ticklength'));
    set(h1,'xtick',[-10:5:10],'xticklabel',{''});
%     letfill = fill(ax1,[-13 -7.3 -7.3 -13],[26 26 27 27],[0.98 0.98 0.98],'LineStyle','none');
%     alpha(letfill,0.8)
%     text(ax1,-12.8,26.5,'Rectified Images','interpreter','latex','fontsize',19);
    text(ax1,-13,25.3,['Time: ',num2str(round(time(i))),' s'],'interpreter','latex','fontsize',19);
    set(ax2,'YDir','reverse')
    
    ax3 = axes('Position',[0.095 0.1 0.8 0.4]);
    h3 = pcolor(ax3,ytemp,xtemp,ztemp-Tinfo.tide-.03);
    hold on
    scatter(pg.xyz(2,:),pg.xyz(1,:),50,'r','fill','MarkerEdgeColor','k','LineWidth',1)
    plot([-13 13],[27.6 27.6],'r','LineStyle','--','LineWidth',2)
    set(h3,'EdgeColor','none')
    axis equal
    ylim([26 35])
    xlim([-13 13])
    hold on
    ax4 = axes('Position',[0.095 0.1 0.8 0.4]);
    Sthresh = flipud(double(stereo_thresh.Ithresh(:,:,i)));
    Sthresh(Sthresh==0)=NaN;
%     h4 = pcolor(stereo_thresh.Y',stereo_thresh.X',fliplr(Sthresh'))
%     set(h4,'EdgeColor','none')
%     alpha(h4,0.5)
    hold on
    axis equal
    ylim([26 35])
    xlim([-13 13])
    shading flat
    
    linkaxes([ax3,ax4])
    ax4.Visible = 'off';
    ax4.XTick = [];
    ax4.YTick = [];
    colormap(ax4,'hot')
    
    colormap(ax3,flipud(cmocean('deep')));
    hc = colorbar(ax3,'Position', [0.91 0.11 0.03 0.32]);
    caxis(ax3,[-0.13001 0.13001])
    box on
    h1=ax3;
    set(h1,'ydir','reverse')
    set(h1,'tickdir','out','xminortick','off','yminortick','off');
    set(h1,'ticklength',0.6*get(h1,'ticklength'));
    set(h1,'fontsize',16);
    text(ax3,13.3,26.5,'$\eta$ (m)','interpreter','latex','fontsize',19);
    xlabel(ax3,'Alongshore (m)','interpreter','latex','fontsize',19);
    ylabel(ax3,'Cross-shore (m)','interpreter','latex','fontsize',19);
%     letfill = fill(ax3,[-13 -5.3 -5.3 -13],[26 26 27 27],[0.98 0.98 0.98],'LineStyle','none');
%     alpha(letfill,0.8)
%     text(ax3,-12.8,26.5,'Stereo Reconstructions','interpreter','latex','fontsize',19);
    set(ax4,'ydir','reverse')
    
    
    Sname = [figfolder,'wavecrest_threshold_flipped_',L(i,:)];
    print(Sname,'-dpng')

    writeVideo(v,getframe(gcf))  
   
    
%     pause(2)
    clf
%     count = count+icount;
end

close(v) 
end
% disp('saved movie')