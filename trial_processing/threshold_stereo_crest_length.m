function threshold_stereo_crest_length(Tinfo,gamma,minarea,samprate,pltflag,ixlim,iylim,idxdy)
% This code will plot rectified images using code from
% D_gridGenExpampleRect.m in the CIRN-Quantitative-Coastal-Imaging-Toolbox


% Set up paths and clear workspace
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\cameras'))
addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\CIRN-Quantitative-Coastal-Imaging-Toolbox/X_CoreFunctions/'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gamma = 0.4;%0.7%0.6;%0.6 % image threshold
% minarea = 30;%800 % min area of object
% samprate = 2; %Hz
% pltflag = 1;

% ixlim=[25 31.5];
% iylim=[-14.5 14.5];

%% STEP 1: Create paths, files and naming
% general path and names
datapath    = 'E:/';

Tinfo = trial_files(Tinfo);

% Stereo Reconstructions
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

%% Load data

bathymetry = load('E:\data\processed\lidar\Riegl\TRC_bathymetry_estimate_line.mat');
xh = bathymetry.xp;
h = Tinfo.tide-movmean(bathymetry.h,50);
figure; plot(xh,Tinfo.tide-bathymetry.h); hold on; plot(xh,h)

eval(['!mkdir ',Tinfo.figfolder,'stereo'])
subname = '';
Tinfo.cam.regy = iylim;
% F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_resx',num2str((idxdy)*100),'cm_resy',num2str((idxdy)*100),'cm',subname,'.mat'];
if Tinfo.filt == 1
    F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_resx',num2str((Tinfo.cam.dx)*100),'cm_resy',num2str((Tinfo.cam.dx)*100),'cm_filtered.mat'];
elseif Tinfo.filt == 0
    F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_resx',num2str((Tinfo.cam.dx)*100),'cm_resy',num2str((Tinfo.cam.dx)*100),'cm.mat'];
end
    load(F1,'x','y','z');

%% plot

imageno =  Tinfo.cam.imagestart:Tinfo.cam.Hz/samprate:(Tinfo.cam.imagestart+Tinfo.cam.numframes-1);
yalen = [];
yclen = [];
cang = [];
yim = [];
xreg = ixlim;
major_ax = [];
minor_ax = [];
orientation = [];
centroid = [];
yctlen = [];
% iylim = [y(1,1) y(end,1)];

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
% if remap == 1
%     thresh=roundgridfun(X,Y,thresh,x,y,@mean);
% end

imfreq = Tinfo.cam.Hz/samprate;
if pltflag == 1
	figure('units','inches','position',[1 1 6 8],'color','w')
end

for i = 1:size(Z,3)/imfreq
    
    Zim = squeeze(Z(:,:,i*imfreq))-msse;
    Zim(Zim>0.75) = NaN;
    Zim = smoothdata(Zim,1,'movmedian',[3 3],'omitnan');%[4 4]
    Zim = smoothdata(Zim,2,'movmedian',[1 1],'omitnan');% [1 1]
    
%     if remap == 1
%         Zim=roundgridfun(X,Y,Zim,x,y,@mean);
%     end
%     Irbw=roundgridfun(X(1,:),Y(:,1),Irbwo,x,y,@mean);
    
    
    % apply the threshold to 
    Iz   = 0*ones(size(Zim));
    Iz(Zim>thresh) = 1;

    Ioz = bwareaopen(Iz,minarea,8);
    Ithresh(:,:,i) = Ioz;
%     Iob = bwareafilt(Iob,5);
    [L,n] = bwlabel(Ioz);
    
    shapeheight = [];
    for ni = 1:n
        oztemp      = 0.*Ioz;
        oztemp(L==ni)=1;
        props = regionprops(oztemp, 'BoundingBox', 'MajorAxisLength', 'Orientation', 'Centroid', 'MinorAxisLength');
%         shapeheight(ni) = props.BoundingBox(4)*resy;
        aedge(1)    = (props.BoundingBox(2)*resy)+Y(1,1);
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
        for ro = 1 : size(oztemp, 1)
            if sum(oztemp(ro,:))>0
                count = count+1;
                B(1) = (find(oztemp(ro, :), 1, 'first')*resy)+xreg(1);
                B(2) = (find(oztemp(ro, :), 1, 'last')*resy)+xreg(1);
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
        
        % plot
        ax2 = axes('Position',[0.15 0.1 0.45 0.8]);
        if remap == 1
            pcolorjw(x,y,Zim)
        else
            pcolorjw(X,Y,Zim)
        end        
        %     imagesc(nu,nv,Irbw)
        hold on
        shading interp
        colormap(ax2,flipud(cmocean('deep')));
        %     hc = colorbar('Position', [0.88 0.1 0.03 0.73]);
        caxis([-0.2 0.2])
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
        hc = colorbar(ax2,'Position', [0.1 0.1 0.04 0.3]);
        
        ax3 = axes('Position',[0.5 0.1 0.45 0.8]);
        if remap == 1
            pcolor(ax3,x,y,double(Ioz))
        else
            pcolor(ax3,X,Y,double(Ioz))
        end  
        hold on
        shading faceted
        colormap(ax3,'gray');
        caxis([0 1])
        axis equal
        %     grid on
        box on
        shading interp
        xlim(ixlim)
        ylim(iylim)
        h1=gca;
%         set(h1,'ydir','reverse')
        set(h1,'tickdir','out','xminortick','on','yminortick','on');
        set(h1,'ticklength',1*get(h1,'ticklength'));
        set(h1,'fontsize',16);
        set(h1,'ytick',[-10:5:10],'yticklabel',{'' '' '' '' ''});
        text(ax2,13.95,26.25,'$\eta$ (m)','interpreter','latex','fontsize',20);
        xlabel('$x$ (m)','interpreter','latex','fontsize',20);
        text(ax2,25,15.5,['Gamma~ = ',num2str(gamma)],'interpreter','latex','fontsize',15);
        text(ax2,25,14.5,['Min Area~~ = ',num2str(minarea),' pix'],'interpreter','latex','fontsize',15);
        text(ax2,25,13.5,['Frame ~~~~~ = ',num2str(imageno(i))],'interpreter','latex','fontsize',15);

        sname = ['thresh_crest_stereo_',num2str(i,'%04.f')];
        print([Tinfo.figfolder,'stereo/',sname],'-dpng')
        pause(0.1)
        clf
    end
    
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
psname = [Tinfo.savefolder,'ylen_stereo_gamma',num2str(gamma*10),'_minarea',num2str(minarea),'_samprate',num2str(samprate),subname,'.mat'];
eval(['save -v7.3 ',psname,' yalen',' yim',' resx',' resy',' xreg',' yclen',' cang',' yctlen',' elipgeom',' xycent']);
% display('Not saving ylen')

psname = [Tinfo.savefolder,'ylen_stereo_gamma',num2str(gamma*10),'_minarea',num2str(minarea),'_samprate',num2str(samprate),subname,'_selectedregions.mat'];
eval(['save -v7.3 ',psname,' X',' Y',' Ithresh']);

end