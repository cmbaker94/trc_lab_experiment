
% function threshold_image_crest_length(Tinfo,thresh,minarea,samprate,pltflag,ixlim,iylim,idxdy,ixsel)% This code will plot rectified images using code from
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

%%
% Trial info
Tinfo.Hs = 0.3;%Hs(op);
Tinfo.Tp = 2;%Tp(op);
Tinfo.tide = 1.07;% h(op);
Tinfo.spread = 30;%sprd(op);
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

%%

bathymetry = load('E:\data\processed\lidar\Riegl\TRC_bathymetry_estimate_line.mat');
xh = bathymetry.xp;
h = Tinfo.tide-movmean(bathymetry.h,50);

subname = '';
Tinfo.cam.regy = iylim;
% F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_resx',num2str((idxdy)*100),'cm_resy',num2str((idxdy)*100),'cm',subname,'.mat'];
if Tinfo.filt == 1
    F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_resx',num2str((Tinfo.cam.dx)*100),'cm_resy',num2str((Tinfo.cam.dx)*100),'cm_filtered.mat'];
elseif Tinfo.filt == 0
    F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_resx',num2str((Tinfo.cam.dx)*100),'cm_resy',num2str((Tinfo.cam.dx)*100),'cm.mat'];
end
stereo = load(F1,'x','y','z');

[val,idxmin] = min(abs(xreg(1)-x(1,:)));
[val,idxmax] = min(abs(xreg(2)-x(1,:)));
[val,idymin] = min(abs(-13.1-y(:,1)));
[val,idymax] = min(abs(13.2-y(:,1)));

X = x(idymin:idymax,idxmin:idxmax);
Y = y(idymin:idymax,idxmin:idxmax);
Z = squeeze(z(idymin:idymax,idxmin:idxmax,:));

msse = nanmean(nanmean(Z,3),1);
msse = repmat(msse,size(Z,1),1);
resx = X(1,2)-X(1,1);
resy = Y(2,1)-Y(1,1);
clear x y z

remap = 0;
% if remap == 1
%     resx = 0.2;
%     resy = 0.25;
%     xtemp = xreg(1):resx:xreg(2);
%     ytemp = Y(1,1):resy:Y(end,1);
%     y = repmat(ytemp',1,length(xtemp));
%     x = repmat(xtemp,length(ytemp),1);
% end

[val,idhmin] = min(abs(xreg(1)-xh));
[val,idhmax] = min(abs(xreg(2)-xh));
hmat = repmat(h(idhmin:idhmax)',size(Z,1),1);
thresh = (gamma/2).*hmat;

%% Load data
eval(['!mkdir ',Tinfo.figfolder,'images'])

odir = [Tinfo.savefolder(1:74),'orthos\','x',num2str(ixlim(1)),'to',num2str(round(ixlim(2))),'_y',num2str(round(iylim(2))),'_res',num2str(idxdy*100),'cm_stereo'];

imageno =  Tinfo.cam.imagestart:Tinfo.cam.Hz/samprate:(Tinfo.cam.imagestart+Tinfo.cam.numframes-1);
yclen = [];
cang = [];
yalen = [];
yim = [];
major_ax = [];
minor_ax = [];
orientation = [];
centroid = [];
yctlen = [];

xreg = ixsel;
resx = idxdy;
resy = idxdy;

xtemp = ixlim(1):resx:ixlim(2)-resx;
ytemp = iylim(1):resy:iylim(2)-resy;

[~,idxmin] = min(abs(ixsel(1)-xtemp));
[~,idxmax] = min(abs(ixsel(2)-xtemp));
xtemp = xtemp(idxmin:idxmax);

Y = repmat(ytemp',1,length(xtemp));
X = repmat(xtemp,length(ytemp),1);

imfreq = Tinfo.cam.Hz/samprate;

%%
if pltflag == 1
	figure('units','inches','position',[1 1 5 8],'color','w')
end

for i = 1:length(imageno)  
    imagefile = fullfile(odir,['c2_',sprintf('%05d',imageno(i)),'.tiff']);
    IM = double(rgb2gray(imread(imagefile)))/255;
    
    Irbw = smoothdata(IM,1,'movmedian',[3 3],'omitnan');
    Irbw = smoothdata(Irbw,2,'movmedian',[1 1],'omitnan');
    Irbw = Irbw(:,idxmin:idxmax);
    
    % apply the threshold to
    Ib   = imbinarize(Irbw,thresh);
   
    Iob = bwareaopen(Ib,minarea,8);
    Ithresh(:,:,i) = Iob;
%     Iob = bwareafilt(Iob,5);
    [L,n] = bwlabel(Iob);
    
    shapeheight = [];
    for ni = 1:n
        obtemp = 0.*Iob;
        obtemp(L==ni)=1;
        props = regionprops(obtemp, 'BoundingBox', 'MajorAxisLength', 'Orientation', 'Centroid', 'MinorAxisLength');
        aedge(1)    = (props.BoundingBox(2)*resy)+Y(1);
        aedge(2)    = aedge(1)+(props.BoundingBox(4)*resy);
        yalen        = [yalen; aedge];
        yclen       = [yclen; props.MajorAxisLength*resy];
        cang        = [cang; props.Orientation];
        yim         = [yim imageno(i)];
        major_ax    = [major_ax props.MajorAxisLength*resy];
        minor_ax	= [minor_ax props.MinorAxisLength*resy];
        orientation = [orientation props.Orientation];
        ctemp       = props.Centroid*resy+[xreg(1) iylim(1)];
        centroid    = [centroid; ctemp];

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
        xycent{length(yim)} = C(1:floor(length(C)/5):end,:);
        d = diff(C);
        totlen = sum(sqrt(sum(d.^2,2)));
        yctlen = [yctlen; totlen];
        clear C
    end
    if pltflag == 1
%         figure('units','inches','position',[1 1 5 8],'color','w')
        % plot
        ax2 = axes('Position',[0.15 0.1 0.4 0.8]);
        pcolor(X,Y,Irbw)       
        %     imagesc(nu,nv,Irbw)
        hold on
        props = regionprops(Iob, 'BoundingBox', 'MajorAxisLength', 'Orientation', 'Centroid', 'MinorAxisLength');
        
        t = linspace(0,2*pi,50);
        for ni = 1:length(props)
            a = (props(ni).MajorAxisLength*resy)/2;
            b = (props(ni).MinorAxisLength*resy)/2;
            Xc = props(ni).Centroid(1)*resy+xreg(1);
            Yc = props(ni).Centroid(2)*resy+iylim(1);
            phi = deg2rad(-props(ni).Orientation);
            x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
            y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
            plot(x,y,'r','Linewidth',5)
        end    
        shading interp
        colormap('gray');
        %     hc = colorbar('Position', [0.88 0.1 0.03 0.73]);
        caxis([0 1])
        axis equal
        box on
        xlim(ixlim)
        ylim(iylim)
        h1=gca;
%         set(h1,'ydir','reverse')
        set(h1,'tickdir','out','xminortick','on','yminortick','on');
        set(h1,'ticklength',1*get(h1,'ticklength'));
        set(h1,'fontsize',16);
        set(h1,'ytick',[-10:5:10],'yticklabel',{'-10' '-5' '0' '5' '10'});
        ylabel('$y$ (m)','interpreter','latex','fontsize',20);
        xlabel('$x$ (m)','interpreter','latex','fontsize',20);
        
        ax3 = axes('Position',[0.46 0.1 0.45 0.8]);
        pcolor(X,Y,double(Iob)) 
        hold on
        shading faceted
        colormap('gray');
        % hc = colorbar('Position', [0.575 0.1 0.02 0.3]);
        caxis([0 1])
        axis equal
        %     grid on
        box on
        shading interp
        xlim(ixsel)
        ylim(iylim)
        h1=gca;
%         set(h1,'ydir','reverse')
        set(h1,'tickdir','out','xminortick','on','yminortick','on');
        set(h1,'ticklength',1*get(h1,'ticklength'));
        set(h1,'fontsize',16);
        set(h1,'ytick',[-10:5:10],'yticklabel',{'' '' '' '' ''});
        text(ax2,13.95,26.25,'$\eta$ (m)','interpreter','latex','fontsize',20);
        xlabel('$x$ (m)','interpreter','latex','fontsize',20);
        text(ax2,25,17.5,['Threshold~ = ',num2str(thresh)],'interpreter','latex','fontsize',15);
        text(ax2,25,16.5,['Min Area~~ = ',num2str(minarea),' pix'],'interpreter','latex','fontsize',15);
        text(ax2,25,15.5,['Frame ~~~~~ = ',num2str(imageno(i))],'interpreter','latex','fontsize',15);

        sname = ['thresh_crest_',num2str(i,'%04.f')];
        print([Tinfo.figfolder,'images/',sname],'-dpng')
        pause(0.1)
        clf
    end
    clear IM Irbw
end

figure;
histogram(abs(yalen(:,1)-yalen(:,2)))
sname = 'histogram_image_alongshore';
print([Tinfo.figfolder,sname],'-dpng')

figure;
histogram(yclen)
sname = 'histogram_image_crest';
print([Tinfo.figfolder,sname],'-dpng')

figure;
histogram(yctlen)
sname = 'histogram_image_crest_transect';
print([Tinfo.figfolder,sname],'-dpng')

figure;
histogram(cang)
sname = 'histogram_image_orient';
print([Tinfo.figfolder,sname],'-dpng')

elipgeom.major_ax = major_ax;
elipgeom.minor_ax = minor_ax;
elipgeom.orientation = orientation;
elipgeom.centroid = centroid;

subname = ['_x',num2str(round(xreg(1))),'to',num2str(round(xreg(2))),'_ycrest'];
psname = [Tinfo.savefolder,'ylen_thresh',num2str(thresh*100),'_minarea',num2str(minarea),'_samprate',num2str(samprate),subname,'.mat'];
eval(['save -v7.3 ',psname,' yalen',' yim',' resx',' resy',' xreg',' yclen',' cang',' yctlen',' elipgeom',' xycent']);
% display('Not saving ylen')

psname = [Tinfo.savefolder,'ylen_thresh',num2str(thresh*100),'_minarea',num2str(minarea),'_samprate',num2str(samprate),subname,'_selectedregions.mat'];
eval(['save -v7.3 ',psname,' X',' Y',' Ithresh']);

% end
