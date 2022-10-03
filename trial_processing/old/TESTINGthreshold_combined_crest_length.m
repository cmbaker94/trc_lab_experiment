% function threshold_combined_crest_length(Tinfo,gammab,IMthreshamp,minareaIMST,minarea,samprate,pltflag,ixlim,iylim,idxdy,ixsel)
% GOAL: Identify breaking crests and measure properties.
% This code will threshold imagery based on pixel intensity, constrain
% these estimates with stereo, and then compute crest properties.

close all
clear all
clc
idxdy=0.05;
ixlim=[25 33.5];%25
iylim=[-14 14];
imlen=0;
% 
% % % threshold images
gammab = 0.24;
minareaIMST = 40;% min area of crest for stereo and imagery
minarea = 100; % min area of crest combined
IMthreshamp =  1.5; % adjust the mean cross-shore pixel intensity for thresholding imagery
IMthresh = 0.3;
samprate = 8;  %Hz
pltflag = 1;
ixsel = [27 31.5];

Tinfo.Hs = 0.25;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
Tinfo.spread = 0;
Tinfo = trial_files(Tinfo);

% filter data
Tinfo.filt = 1; % 1 if true

% Set up paths and clear workspace
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\cameras'))
addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\CIRN-Quantitative-Coastal-Imaging-Toolbox/X_CoreFunctions/'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

%% STEP 1: Create paths, files and naming
% general path and names
datapath    = 'E:/';

Tinfo = trial_files(Tinfo);

% Stereo Reconstructions
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

%% Load bathymetry

bathymetry = load('E:\data\processed\lidar\Riegl\TRC_bathymetry_estimate_line.mat');
xh = bathymetry.xp;
h = Tinfo.tide-movmean(bathymetry.h,50);
figure; plot(xh,Tinfo.tide-bathymetry.h); hold on; plot(xh,h)

%% Load stereo

subname = '';
Tinfo.cam.regy = iylim;
if Tinfo.filt == 1
    F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_resx',num2str((Tinfo.cam.dx)*100),'cm_resy',num2str((Tinfo.cam.dx)*100),'cm_filtered.mat'];
elseif Tinfo.filt == 0
    F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_resx',num2str((Tinfo.cam.dx)*100),'cm_resy',num2str((Tinfo.cam.dx)*100),'cm.mat'];
end

if pltflag == 1
    stereo.orig = load(F1,'x','y','z');
    
    [~,idxmin] = min(abs(ixsel(1)-stereo.orig.x(1,:)));
    [~,idxmax] = min(abs(ixsel(2)-stereo.orig.x(1,:)));
    [~,idymin] = min(abs(-10-stereo.orig.y(:,1)));
    [~,idymax] = min(abs(10-stereo.orig.y(:,1)));
    
    stereo.x = stereo.orig.x(:,idxmin:idxmax);
    stereo.y = stereo.orig.y(:,idxmin:idxmax);
    stereo.z = squeeze(stereo.orig.z(:,idxmin:idxmax,:));
    stereo.orig.msse = movmean(stereo.orig.z,20*8,3,'omitnan');
else
    stereo = load(F1,'x','y','z');
    
    [~,idxmin] = min(abs(ixsel(1)-stereo.x(1,:)));
    [~,idxmax] = min(abs(ixsel(2)-stereo.x(1,:)));
    [~,idymin] = min(abs(-10-stereo.y(:,1)));
    [~,idymax] = min(abs(10-stereo.y(:,1)));
    
    stereo.x = stereo.x(:,idxmin:idxmax);
    stereo.y = stereo.y(:,idxmin:idxmax);
    stereo.z = squeeze(stereo.z(:,idxmin:idxmax,:));
end

mssemvavg = 10*8; % moving mean over 80 frames (10 sec or 5 waves)
mssetemp = movmean(stereo.z,10*8,3,'omitnan'); % moving mean over 10 waves
for i = 1:mssemvavg/2
    mssetemp(:,:,i) = mssetemp(:,:,mssemvavg/2);
    mssetemp(:,:,end-i+1) = mssetemp(:,:,end-mssemvavg/2);
end
stereo.msse = mssetemp;

% % check moving mean for water surface elevation
% figure
% for i = 1:300;%size(mssetemp,3)
%     pcolor(mssetemp(:,:,i))
%     shading interp
%     colorbar
%     caxis([1 1.2])
%     pause(0.5)
% end

[val,idhmin] = min(abs(ixsel(1)-xh));
[val,idhmax] = min(abs(ixsel(2)-xh));
hmat = repmat(h(idhmin:idhmax)',size(stereo.z,1),1);
SRthresh = (gammab/2).*hmat;

%% Load imagery
eval(['!mkdir ',Tinfo.figfolder,'images'])

odir = [Tinfo.savefolder(1:74),'orthos\','x',num2str(ixlim(1)),'to',num2str(round(ixlim(2))),'_y',num2str(round(iylim(2))),'_res',num2str(idxdy*100),'cm_stereo'];

imageno =  Tinfo.cam.imagestart:Tinfo.cam.Hz/samprate:(Tinfo.cam.imagestart+Tinfo.cam.numframes-1);

xtemp = ixlim(1):idxdy:ixlim(2);
ytemp = iylim(1):idxdy:iylim(2);

[~,idxmin] = min(abs(ixsel(1)-xtemp));
[~,idxmax] = min(abs(ixsel(2)-xtemp));
xorig = xtemp;
xtemp = xtemp(idxmin:idxmax);

Y = repmat(ytemp',1,length(xtemp));
X = repmat(xtemp,length(ytemp),1);

imfreq = Tinfo.cam.Hz/samprate;

%% Find threshold imagery
I = [];
display('change back from moving mean to moving median?')
for i = 1:length(imageno)
    imagefile = fullfile(odir,['c2_',sprintf('%05d',imageno(i)),'.tiff']);
    IM = double(rgb2gray(imread(imagefile)))/255;
    if i==1
        I = IM;
    else
        I = I+IM;
    end
end

Iavg = nanmean(I/length(imageno),3);
% IMthresh = 0.2;


% 
% imreg = I(:,70:120,1)/length(imageno);
% imreg(imreg>0.45)=NaN;
% IMthresh = squeeze(nanmean(imreg,2));
% IMthresh(IMthresh<0.22)=0.22;
% IMthresh = repmat(IMthresh,1,length(xtemp))*IMthreshamp;
%clear IM I

%% prep

Iprop.yclen = [];
Iprop.cang = [];
Iprop.yalen = [];
Iprop.yim = [];
Iprop.elipgeom.major_ax = [];
Iprop.elipgeom.minor_ax = [];
Iprop.elipgeom.orientation = [];
Iprop.elipgeom.centroid = [];
Iprop.yctlen = [];
Iprop.exitregflag = [];
Iprop.crestends = [];
Ithresh = [];
IMprop = Iprop;
SRprop = Iprop;

%%

v = VideoWriter([Tinfo.figfolder,'wavecrests_threshold.avi']);
v.FrameRate=8; %?????
open(v)

if pltflag == 1
	figure('units','inches','position',[1 1 9.5 8],'color','w')
end

[~,ioff] = nanmin(abs(ixsel(1)-X(1,:)));
display('change back from moving mean to moving median?')
for i = 150:length(imageno)  
    
    imagefile = fullfile(odir,['c2_',sprintf('%05d',imageno(i)),'.tiff']);
    IM = (double(rgb2gray(imread(imagefile)))/255)-Iavg;
    IM(IM<0)=0;
    Iob = imagethresh(IM,IMthresh,minareaIMST,idxmin,idxmax); % threshold imagery
    Zim = squeeze(stereo.z(:,:,i))-squeeze(stereo.msse(:,:,i));
    Ioz = stereothresh(Zim,SRthresh,minareaIMST); % threshold stereo
        
    Icomb = Iob;
    Icomb(flipud(Ioz)==0)=0;
    [Icomb,Iprop] = procthreshprop(Icomb,Iprop,minarea,Y,ixsel,iylim,idxdy,imageno,i);
    [~,IMprop] = procthreshprop(Iob,IMprop,minarea,Y,ixsel,iylim,idxdy,imageno,i);
    [~,SRprop] = procthreshprop(Ioz,SRprop,minarea,Y,ixsel,iylim,idxdy,imageno,i);
    Ioz = bwareaopen(flipud(Ioz),minareaIMST);
    Iob = bwareaopen(Iob,minareaIMST);
    
    Ithresh(:,:,i) = Icomb;
%     Ioz = flipud(Icomb);
%     Iob = Icomb;
    if i < 4*60*8
    if pltflag == 1
        iI = (Iprop.yim==imageno(i));
        numends = sum(Iprop.crestends(iI));
        numcrests = sum(iI);
        ilen = Iprop.yctlen(iI);
        ilen(Iprop.exitregflag(iI)==1) = NaN;
        avglen  = nanmean(ilen);
        
        ax1 = axes('Position',[0.07 0.095 0.27 0.8]);
        imagesc(ax1,xorig,Y(:,1),imread(imagefile));
        hold on
        plot([ixsel(1) ixsel(1)],[Y(1,1) Y(end,1)],'Color',[196/265, 118/265, 245/265, 0.7],'LineWidth',2,'LineStyle','-.')
        plot([ixsel(2) ixsel(2)],[Y(1,1) Y(end,1)],'Color',[196/265, 118/265, 245/265, 0.7],'LineWidth',2,'LineStyle','-.')
        axis image
        axis equal
        ylim([-13,13]);
        xlim(ixlim);
        set(ax1,'YDir','normal')
        
        ax2 = axes('Position',[0.07 0.095 0.27 0.8]);
        Itemp = double(Iob);
        Itemp(Itemp==0)=NaN;
        h2 = pcolor(X,Y,Itemp);
        set(h2,'EdgeColor','none')
        alpha(h2,0.5)
        hold on
%         pltcrestparam(Iob,Itemp,Y,idxdy,ixsel)
        axis equal
        xlim(ixlim)
        ylim([-13.3 13.3])
        linkaxes([ax1,ax2])
        ax2.Visible = 'off';
        ax2.XTick = [];
        ax2.YTick = [];
        colormap(ax2,'hot')
        ylabel(ax1,'$y$~(m)','interpreter','latex','fontsize',20);
        xlabel(ax1,'$x$~(m)','interpreter','latex','fontsize',20);
        h1=ax1;
        set(h1,'tickdir','out','xminortick','off','yminortick','off');
        set(h1,'ticklength',1*get(h1,'ticklength'));
        set(h1,'fontsize',16);
        set(h1,'ytick',[-10:5:10],'yticklabel',{'-10' '-5' '0' '5' '10'})
        text(ax1,25,14,'(a) Image','interpreter','latex','fontsize',20);
        text(ax1,25,15.5,['Time: ',num2str(i/8,'%.2f'),' s ~~~~~~~~~~$\overline{l_c}=~$',...
            num2str(round(avglen,1)),' m',',~~~$N_{\mathrm{crests}}=~$',num2str(numcrests),...
            ',~~~$N_{\mathrm{crest-ends}}=~$',num2str(numends)],...
            'interpreter','latex','fontsize',20);
        
        ax3 = axes('Position',[0.33 0.095 0.27 0.8]);
        h3 = pcolor(ax3,stereo.orig.x,stereo.orig.y,flipud(squeeze(stereo.orig.z(:,:,i))-squeeze(stereo.orig.msse(:,:,i))));
        set(h3,'EdgeColor','none')
        hold on
        plot([ixsel(1) ixsel(1)],[Y(1,1) Y(end,1)],'Color',[196/265, 118/265, 245/265, 0.7],'LineWidth',2,'LineStyle','-.')
        plot([ixsel(2) ixsel(2)],[Y(1,1) Y(end,1)],'Color',[196/265, 118/265, 245/265, 0.7],'LineWidth',2,'LineStyle','-.')
        axis equal
        xlim(ixlim)
        ylim([-13.3 13.3])
        
        ax4 = axes('Position',[0.33 0.095 0.27 0.8]);
        Stemp = double(Ioz);
        Stemp(Stemp==0)=NaN;
        h4 = pcolor(X,Y,Stemp);
        set(h4,'EdgeColor','none')
        alpha(h4,0.5)
        hold on
%         pltcrestparam(Ioz,Stemp,Y,idxdy,ixsel)
        axis equal
        xlim(ixlim)
        ylim([-13.3 13.3])
        shading flat
        linkaxes([ax3,ax4])
        ax4.Visible = 'off';
        ax4.XTick = [];
        ax4.YTick = [];
        colormap(ax4,'hot')
        colormap(ax3,flipud(cmocean('deep')));
        hc = colorbar(ax3,'Position', [0.6 0.095 0.03 0.74]);
        caxis(ax3,[-0.15001 0.15001])
        box on
        h1=ax3;
        set(h1,'tickdir','out','xminortick','off','yminortick','off');
        set(h1,'ticklength',1*get(h1,'ticklength'));
        set(h1,'fontsize',16);
        set(h1,'ytick',[-10:5:10],'yticklabel',{''});
        text(ax3,34.5,12.3,'$\eta$ (m)','interpreter','latex','fontsize',20);
        xlabel(ax3,'$x$ (m)','interpreter','latex','fontsize',20);
        text(ax3,25,14,'(b) Stereo','interpreter','latex','fontsize',20);
        
        
        ax5 = axes('Position',[0.7 0.095 0.27 0.8]);
        imagesc(ax5,xorig,Y(:,1),imread(imagefile));
        hold on
        plot([ixsel(1) ixsel(1)],[Y(1,1) Y(end,1)],'Color',[196/265, 118/265, 245/265, 0.7],'LineWidth',2,'LineStyle','-.')
        plot([ixsel(2) ixsel(2)],[Y(1,1) Y(end,1)],'Color',[196/265, 118/265, 245/265, 0.7],'LineWidth',2,'LineStyle','-.')
        axis image
        axis equal
        ylim([-13,13]);
        xlim(ixlim);
        set(ax5,'YDir','normal')
        
        ax6 = axes('Position',[0.7 0.095 0.27 0.8]);
        Itemp = double(Icomb);
        Itemp(Itemp==0)=NaN;
        h2 = pcolor(X,Y,Itemp);
        set(h2,'EdgeColor','none')
        alpha(h2,0.5)
        hold on
        pltcrestparam(Icomb,Itemp,Y,idxdy,ixsel)
        axis equal
        xlim(ixlim)
        ylim([-13.3 13.3])
        linkaxes([ax5,ax6])
        ax6.Visible = 'off';
        ax6.XTick = [];
        ax6.YTick = [];
        colormap(ax6,'hot')
        xlabel(ax5,'$x$~(m)','interpreter','latex','fontsize',20);
        h1=ax5;
        set(h1,'tickdir','out','xminortick','off','yminortick','off');
        set(h1,'ticklength',1*get(h1,'ticklength'));
        set(h1,'fontsize',16);
        set(h1,'ytick',[-10:5:10],'yticklabel',{''})
        text(ax5,25,14,'(c) Combined','interpreter','latex','fontsize',20);
        
        
%         Sname = [figfolder,'wavecrest_threshold_',L(i,:)];
%         print(Sname,'-dpng')
%         
        writeVideo(v,getframe(gcf))
    end
    clf
    end
    clear IM Irbw Iob Ioz Icomb
end

close(v)

figure;
histogram(abs(Iprop.yalen(:,1)-Iprop.yalen(:,2)))
sname = 'histogram_combined_alongshore';
print([Tinfo.figfolder,sname],'-dpng')

figure;
histogram(Iprop.yclen)
sname = 'histogram_combined_crest';
print([Tinfo.figfolder,sname],'-dpng')

figure;
histogram(Iprop.yctlen)
sname = 'histogram_combined_crest_transect';
print([Tinfo.figfolder,sname],'-dpng')

figure;
histogram(Iprop.cang)
sname = 'histogram_combined_orient';
print([Tinfo.figfolder,sname],'-dpng')

subname = ['_x',num2str(round(ixsel(1))),'to',num2str(round(ixsel(2))),'_ycrest_res',num2str(idxdy),''];
psname = [Tinfo.savefolder,'ylen_combined_gamma',num2str(gammab),'_thresh',num2str(IMthresh*100),'_minarea',num2str(minarea),'_samprate',num2str(samprate),subname,'.mat'];
eval(['save -v7.3 ',psname,' Iprop',' IMprop',' SRprop']);
% display('Not saving ylen')

psname = [Tinfo.savefolder,'ylen_combined_gamma',num2str(gammab),'_thresh',num2str(IMthresh*100),'_minarea',num2str(minarea),'_samprate',num2str(samprate),subname,'_selectedregions.mat'];
eval(['save -v7.3 ',psname,' X',' Y',' Ithresh']);

% end

function [Icomb,Iprop] = procthreshprop(Icomb,Iprop,minarea,Y,ixsel,iylim,idxdy,imageno,i)
% remove waves touching edges of region

Icomb = bwareaopen(Icomb,minarea,8);
% [L,n] = bwlabel(Icomb);
% exitregflag = zeros(n,1);
% for ip = 1:n
%     temp = zeros(size(L,1),1);
%     temp(L(:,[end-1 end])==ip)=1;
%     if max(temp)>0
% %         Icomb(L==ip)=0;
%         exitregflag(ip)=1;
%     end
%     clear temp
%     temp = zeros(size(L,1),1);
%     temp(L(:,[1 2])==ip)=1;
%     if max(temp)>0
% %         Icomb(L==ip)=0;
%         exitregflag(ip)=1;
%     end
%     clear temp
% end
% Iprop.exitregflag = [Iprop.exitregflag; exitregflag];
Iprop = threshprop(Icomb,Iprop,Y,ixsel,iylim,idxdy,imageno,i);
end

function pltcrestparam(Ithresh,Itemp,Y,idxdy,ixsel)
% add crest properties

Dxplot = [];
Dyplot = [];

[Ln,n] = bwlabel(Ithresh);
for ni = 1:n
%     temp = zeros(size(Ln,1),1);
%     temp(Ln(:,[end-1:end])==ni)=1;
%     temp(Ln(:,[1:2])==ni)=1;
    obtemp = zeros(size(Itemp));
    obtemp(Ln==ni)=1;
    count = 0;
    for ro = 1 : size(obtemp, 1)
        if sum(obtemp(ro,:))>0
            count = count+1;
            B(1) = (find(obtemp(ro, :), 1, 'first')*idxdy)+ixsel(1);
            B(2) = (find(obtemp(ro, :), 1, 'last')*idxdy)+ixsel(1);
            C(count, 1) = mean(B);
            C(count, 2) = (ro*idxdy)+min(Y,[],'all');
        end
        clear B
    end
    D(:,2) = C(:,2);
    D(:,1) = movmean(movmedian(C(:,1),14),14);
%     if sum(temp)==0
%         plot(D(:,1),D(:,2),'y','LineWidth',2);%,'Marker','o','MarkerSize',2)
%     end
    Dx = D([1 end],1);
    Dy = D([1 end],2);
    idx=find(Dx<ixsel(1)+0.2 | Dx>ixsel(2)-0.2);
    Dx(idx) = NaN;
    idy=find(Dy<-13.2 | Dy>13.2);
    Dy(idy) = NaN;
    Dxplot = [Dxplot; Dx'];
    Dyplot = [Dyplot; Dy'];
    crestno = length(Dx)-length(idx);
    if crestno == 2
        plot(D(:,1),D(:,2),'y','LineWidth',2)
    end
    clear C D
end
if size(Dyplot,2)>0
    scatter(Dxplot(:,1),Dyplot(:,1),100,'o','b','LineWidth',2)
    scatter(Dxplot(:,2),Dyplot(:,2),100,'o','r','LineWidth',2)
end

props = regionprops(Ithresh, 'MajorAxisLength', 'Orientation', 'Centroid', 'MinorAxisLength');
t = linspace(0,2*pi,50);
for ni = 1:length(props)
    a = (props(ni).MajorAxisLength*idxdy)/2;
    b = (props(ni).MinorAxisLength*idxdy)/2;
    Xc = props(ni).Centroid(1)*idxdy+ixsel(1);
    Yc = props(ni).Centroid(2)*idxdy+min(Y,[],'all');
    phi = deg2rad(-props(ni).Orientation);
    xe = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    ye = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
    hc = plot(xe,ye,'m','Linewidth',1.2,'LineStyle','--');
    hc.Color = [hc.Color 0.75];
%     [~,iymin]=min(ye);
%     [~,iymax]=max(ye);
%     Xe =  [Xc + a*cos(t)*cos(phi) xe(iymax)];
%     Ye = [ye(iymin) ye(iymax)];
    Xe = Xc+[a*cos(phi) -a*cos(phi)];
    Ye = Yc+[a*sin(phi) -a*sin(phi)];
    idx=find(Xe<ixsel(1) | Xe>ixsel(2));
    Xe(idx)=NaN; Ye(idx)=NaN;
%     scatter(Xe,Ye,100,'x','m','LineWidth',2)
end
end
