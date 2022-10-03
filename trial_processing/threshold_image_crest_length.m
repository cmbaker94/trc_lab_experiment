function threshold_image_crest_length(Tinfo,thresh,minarea,samprate,pltflag,ixlim,iylim,idxdy,ixsel)% This code will plot rectified images using code from
% D_gridGenExpampleRect.m in the CIRN-Quantitative-Coastal-Imaging-Toolbox

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

%% Load data
eval(['!mkdir ',Tinfo.figfolder,'images'])

odir = [Tinfo.savefolder(1:74),'orthos\','x',num2str(ixlim(1)),'to',num2str(round(ixlim(2))),'_y',num2str(round(iylim(2))),'_res',num2str(idxdy*100),'cm'];

imageno =  Tinfo.cam.imagestart:Tinfo.cam.Hz/samprate:(Tinfo.cam.imagestart+Tinfo.cam.numframes-1);
yclen = [];
cang = [];
yalen = [];
yim = [];

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
        props = regionprops(obtemp, 'BoundingBox', 'MajorAxisLength', 'Orientation');
        aedge(1)    = (props.BoundingBox(2)*resy)+Y(1);
        aedge(2)    = aedge(1)+(props.BoundingBox(4)*resy);
        yalen        = [yalen; aedge];
        yclen       = [yclen; props.MajorAxisLength*resy];
        cang        = [cang; props.Orientation];
        yim         = [yim imageno(i)];
    end
  
    if pltflag == 1
%         figure('units','inches','position',[1 1 5 8],'color','w')
        % plot
        ax2 = axes('Position',[0.15 0.1 0.4 0.8]);
        pcolor(X,Y,Irbw)       
        %     imagesc(nu,nv,Irbw)
        hold on
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
histogram(cang)
sname = 'histogram_image_orient';
print([Tinfo.figfolder,sname],'-dpng')

subname = ['_x',num2str(round(xreg(1))),'to',num2str(round(xreg(2))),'_ycrest'];
psname = [Tinfo.savefolder,'ylen_thresh',num2str(thresh*100),'_minarea',num2str(minarea),'_samprate',num2str(samprate),subname,'.mat'];
eval(['save -v7.3 ',psname,' yalen',' yim',' resx',' resy',' xreg',' yclen',' cang']);
% display('Not saving ylen')

psname = [Tinfo.savefolder,'ylen_thresh',num2str(thresh*100),'_minarea',num2str(minarea),'_samprate',num2str(samprate),subname,'_selectedregions.mat'];
eval(['save -v7.3 ',psname,' X',' Y',' Ithresh']);

end
