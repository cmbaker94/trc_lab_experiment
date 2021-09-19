function threshold_image_crest_length(Tinfo,thresh,minarea,samprate,pltflag)% This code will plot rectified images using code from
% D_gridGenExpampleRect.m in the CIRN-Quantitative-Coastal-Imaging-Toolbox


% Set up paths and clear workspace
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\cameras'))
addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\CIRN-Quantitative-Coastal-Imaging-Toolbox/X_CoreFunctions/'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% thresh = 0.45;%0.7%0.6;%0.6 % image threshold
% minarea = 30;%800 % min area of object
% samprate = 2; %Hz
% pltflag = 1;

idxdy=0.02;
ixlim=[25 31.5];
iylim=[-14.5 14.5];

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

subname = '';
load([Tinfo.savefolder,'image_stack_x',num2str(ixlim(1)),'to',num2str(round(ixlim(2))),'_res',num2str(idxdy*100),'cm',subname,'.mat']);

%% plot

imageno =  7200:4:11999;
ylen = [];
yim = [];
xreg = [X(1,1) X(1,end)];

remap = 1;
if remap == 1
    resx = 0.25;
    resy = 0.25;
    xreg = xreg;%[25 30];
    xtemp = xreg(1):resx:xreg(2);
    ytemp = Y(1,1):resy:Y(end,1);
    y = repmat(ytemp',1,length(xtemp));
    x = repmat(xtemp,length(ytemp),1);
end

imfreq = Tinfo.cam.Hz/samprate;

figure('units','inches','position',[1 1 5 8],'color','w')

for i = 1:size(IR,3)/imfreq
    
    Irbw = squeeze(IR(:,:,i*imfreq));
    if remap == 1
        Irbw=roundgridfun(X,Y,Irbw,x,y,@mean);
    end
%     Irbw=roundgridfun(X(1,:),Y(:,1),Irbwo,x,y,@mean);
    
    
    % apply the threshold to 
    Ib   = imbinarize(Irbw,thresh);
   
    Iob = bwareaopen(Ib,minarea,8);
%     Iob = bwareafilt(Iob,5);
    [L,n] = bwlabel(Iob);
    
    shapeheight = [];
    for ni = 1:n
        obtemp = 0.*Iob;
        obtemp(L==ni)=1;
        props = regionprops(obtemp, 'BoundingBox');
        shapeheight(ni) = props.BoundingBox(4)*resy;
    end
    ylen = [ylen shapeheight];
    yim = [yim ones(1,n)*imageno(i)];
    
    if pltflag == 1
        
        % plot
        ax2 = axes('Position',[0.15 0.1 0.4 0.8]);
        if remap == 1
            pcolor(x,y,Irbw)
        else
            pcolor(X,Y,Irbw)
        end        
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
        if remap == 1
            pcolor(x,y,double(Iob))
        else
            pcolor(X,Y,double(Iob))
        end  
        hold on
        shading faceted
        colormap('gray');
        % hc = colorbar('Position', [0.575 0.1 0.02 0.3]);
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
        text(ax2,25,17.5,['Threshold~ = ',num2str(thresh)],'interpreter','latex','fontsize',15);
        text(ax2,25,16.5,['Min Area~~ = ',num2str(minarea),' pix'],'interpreter','latex','fontsize',15);
        text(ax2,25,15.5,['Frame ~~~~~ = ',num2str(imageno(i))],'interpreter','latex','fontsize',15);

        sname = ['thresh_crest_',num2str(i,'%04.f')];
        print([Tinfo.figfolder,'images/',sname],'-dpng')
        pause(0.1)
        clf
    end
    
end

figure;
histogram(ylen)
sname = 'histogram_image';
print([Tinfo.figfolder,sname],'-dpng')

subname = ['_x',num2str(round(xreg(1))),'to',num2str(round(xreg(2)))];
psname = [Tinfo.savefolder,'ylen_thresh',num2str(thresh*100),'_minarea',num2str(minarea),subname,'.mat'];
eval(['save -v7.3 ',psname,' ylen',' yim',' resx',' resy',' xreg']);

end
