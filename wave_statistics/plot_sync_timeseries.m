% This code will extract and save time series from the camera dems and
% lidar mat file to compute the directional spectra as an 'array'.

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\cameras\dem'))
addpath(genpath('E:\code\insitu'))
%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Trial info
    Tinfo.Hs = 0.3;
    Tinfo.Tp = 2;
    Tinfo.tide = 1.07;
    Tinfo.spread = 40;
    
    % for camera 2
    UV06 = [1.6868867 0.93791529]*10^3; 
    UV11 = [1.679672566 1.139658546]*10^3;

    
    Tinfo = trial_files(Tinfo);
    
    % general path and names
    Tinfo.cam = TRC_camera_info(Tinfo.cam);
    
    % Data and figure storage
    [Tinfo] = wc_comp_store(Tinfo); 

    [camera.time] = cam_time(Tinfo);


%% Insitu

% load insitu data
datapath = 'E:\';
F1          = matfile([datapath,'data\processed\insitu\',Tinfo.sz.tdate,'\',Tinfo.sz.tdate,'-insitu.mat']);
pg          = F1.press;
waveg       = F1.wg;
Tinfo_insitu       = F1.Tinfo;
Hz          = pg.Hz;
maxfac = 1.2;
offset  = 0.05;

starttemp   = datenum(Tinfo_insitu.time(1:end-3));
endtemp     = starttemp+((length(waveg.wg1)-1)*datenum(0,0,0,0,0,1/Hz));
timetemp    = starttemp:datenum(0,0,0,0,0,1/Hz):endtemp;

% find overlapping range for insitu
[temp,istart] = nanmin(abs(camera.time(1)-timetemp));
[temp,iend] = nanmin(abs(camera.time(end)-timetemp));

press.time = timetemp(istart:iend)';
press.name  = fieldnames(pg.xyzd)';

% offset for time sync
istart = istart+Tinfo.offsets(1);
iend = iend+Tinfo.offsets(1);

filtramp = 1000;
% compute bulk statistics
for i = 1:length(press.name)
    % load pressure data
    eval(['press.press(:,i)   = pg.',press.name{i},'(istart:iend);'])
    press.eta(:,i) = press2sse_timeseries(press.press(:,i),Hz,offset,maxfac,filtramp);
    eval(['press.xyz(:,i)     = pg.xyzd.',press.name{i},'(1,1:3);']) % z-elevation of sensor in tank coordinates
    press.loc{i} = ['p',sprintf('%02d',str2double(press.name{i}(6:end)))];
end

clear pg *temp T1 T2

% Compute bulk statistics from wave gages
% get wave gage names
wg.names	= fieldnames(waveg.xyzd)';
wg.time     = press.time;
% compute bulk statistics
for i = 1:length(wg.names)
    % load wave gage data
    eval(['wg.z(:,i)                     = waveg.',wg.names{i},'(istart:iend);'])
    eval(['wg.xyz(:,i)              = waveg.xyzd.',wg.names{i},'(1,1:3);']) 
    wg.loc{i} = ['wg',sprintf('%02d',str2double(wg.names{i}(3:end)))];
end

clear i* waveg pg T1 T2

%% STEP 2: Load Stereo Reconstruction Data
camsource = 1;

if camsource == 1
    cam             = load([Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_res',num2str((Tinfo.cam.dx)*100),'cm_filtered.mat']); 
    for i = 1:size(press.xyz,2)
        [temp,ix] = nanmin(abs(cam.x(1,:)-press.xyz(1,i)));
        [temp,iy] = nanmin(abs(cam.y(:,1)-press.xyz(2,i)));
        camera.xy(1,i) = cam.x(1,ix);
        camera.xy(2,i) = round(cam.y(iy,1),2);
        x = squeeze(cam.z(iy,ix,:));
        camera.z(:,i)  = x;
        camera.loc{i} = ['c',sprintf('%02d',str2double(press.name{i}(6:end)))];
    end
elseif camsource == 2
    F1 = matfile([Tinfo.cam.datafolder,'dem_pt_sz_array_400cm2.mat']);
    camera.xy(1,:) = F1.xloc;
    camera.xy(2,:) = F1.yloc;
    camera.z = F1.eta_pt_mean';
    % changing to match rest
%     camera.inst =
%     "press1"    "press2"    "press3"    "press4"    "press5"    "press6"    "press7"    "press8"    "press9"    "press10" "press11"    "press12"
end

clear cam T1 T2 fn*

%% STEP 3: Load VeloLidar
load([datapath,'data\processed\lidar\Velodyne\',Tinfo.lidar.tdate,'\gridded\',Tinfo.lidar.tdate,'_Velodyne-HDL-32-Data_gridded_filtered.mat']);

lidar.Hz = 10;

starttemp   = datenum(Tinfo.lidar.tdate,'yyyy-mm-dd-HH-MM-SS');
endtemp     = datenum(Tinfo.lidar.tdate,'yyyy-mm-dd-HH-MM-SS')+(length(squeeze(griddedData.z(1,1,:)))*datenum(0,0,0,0,0,1/lidar.Hz));
timetemp  = starttemp:datenum(0,0,0,0,0,1/lidar.Hz):endtemp;

[temp,istart] = nanmin(abs(camera.time(1)-timetemp));
[temp,iend] = nanmin(abs(camera.time(end)-timetemp));

lidar.time = squeeze(timetemp(istart:iend));
clear *temp

% offset for time sync
istart = istart+Tinfo.offsets(2);
iend = iend+Tinfo.offsets(2);

for i = 1:size(press.xyz,2)
    [temp,ix] = nanmin(abs(griddedData.xvec-press.xyz(1,i)));
    [temp,iy] = nanmin(abs(griddedData.yvec-press.xyz(2,i)));
    lidar.xy(1,i) = griddedData.xvec(ix);
    lidar.xy(2,i) = griddedData.yvec(iy);
    x = squeeze(griddedData.z(iy,ix,istart:iend));
    lidar.z(:,i)  = x;
    lidar.loc{i} = ['l',sprintf('%02d',str2double(press.name{i}(6:end)))];
end
clear griddedData T1 T2 fn*

%% plot
tiffpath = ['F:\Metashape\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1\',Tinfo.cam.trimname,'dems\'];
imagepath = ['H:\TRC_Fall_Experiment\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1_JPEG\'];
% tiffpath = ['D:\TRC_Fall_Experiment\photoscan\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1\',Tinfo.cam.trimname,'\dems\'];
% imagepath = ['D:\TRC_Fall_Experiment\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1_JPEG\'];
% imagepath = ['D:\TRC_Fall_Experiment\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1_JPEG\c2\'];

imageno = Tinfo.cam.imagestart:Tinfo.cam.imagestart+Tinfo.cam.numframes;
h.presson = nanmean(press.press(:,3));%+0.05;
h.pressoff = nanmean(press.press(:,9));%+0.05;
if camsource == 1
    h.cameraon = nanmean(camera.z(:,3));
    h.cameraoff = nanmean(camera.z(:,9)); 
elseif camsource == 2
    h.cameraon = nanmean(camera.z(:,11));
    h.cameraoff = nanmean(camera.z(:,6)); 
end
    
% h.camera = nanmean(nanmean(camera.z));
h.lidaron = nanmean(lidar.z(:,3));
h.lidaroff = nanmean(lidar.z(:,9)); 
% h.lidar = nanmean(nanmean(lidar.z));
g = 9.81;
rho = 1000;
count = 0;


%%

sname =['timesync_image_vid_25x5cm_mod_imagestereo_timeseries'];
v = VideoWriter([Tinfo.figfolder,'\',sname,'.avi']);
v.FrameRate=8;
open(v)

figure('units','inches','position',[1 1 16 9],'color','w')
reg = 4;
for i = (1+reg)*8:length(camera.time)-(1+reg)*8 
    count = count+1;
    tifffile    = [tiffpath,Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_frames_',Tinfo.cam.imagerange,'_DEM_',num2str(i),'.tif'];
    % read in point cloud
    [xtemp,ytemp,ztemp,resxtemp,resytemp] = read_tiff(tifffile);
    ytemp = fliplr(ytemp);
    ztemp = flipud(ztemp);
    
    % need to make a 2d for finding nans
    ytemp = repmat(ytemp',1,length(xtemp));
    xtemp = repmat(xtemp,length(ytemp),1);

    %read image
    imagefile = [imagepath,getfield(dir([imagepath,'Movie1_Scene1_c2_',sprintf('%05d',imageno(i)),'_*.jpg']),'name')];
    IM = imread(imagefile);
    
    % timeseries 
    [temp,iI] = nanmin(abs(camera.time(i)-press.time));
    [temp,iL] = nanmin(abs(camera.time(i)-lidar.time));
    
    It = press.time(iI-(reg*100):iI+reg*100);
    Ito = (It-camera.time(i))*24*3600;
    Izon = (press.press(iI-reg*100:iI+reg*100,3)-h.presson)/(g*rho);
    Izoff = (press.press(iI-reg*100:iI+reg*100,9)-h.pressoff)/(g*rho);
    Ifzon = press.eta(iI-reg*100:iI+reg*100,3);
    Ifzoff = press.eta(iI-reg*100:iI+reg*100,9);
    
    Ct = camera.time(i-reg*8:i+reg*8);
    Cto = (Ct-camera.time(i))*24*3600;
    Czon = camera.z(i-reg*8:i+reg*8,3)-h.cameraon;
    Czoff = camera.z(i-reg*8:i+reg*8,9)-h.cameraoff;
    
    Lt = lidar.time(iL-reg*10:iL+reg*10);
    Lto = (Lt-camera.time(i))*24*3600;
    Lzon = lidar.z(iL-reg*10:iL+reg*10,3)-h.lidaron;
    Lzoff = lidar.z(iL-reg*10:iL+reg*10,9)-h.lidaroff;
    
    % plot
    ax1 = axes('Position',[0.053 0.56 0.55 0.4]);
    imagesc(IM)
    hold on
    scatter([UV06(1) UV11(1)],[UV06(2) UV11(2)],40,'r','fill')
    text(ax1,UV06(1)-55,UV06(2)-70,'p06','Color','w','interpreter','latex','fontsize',15)
    text(ax1,UV11(1)-55,UV11(2)-70,'p11','Color','w','interpreter','latex','fontsize',15)
    ylim([350,1800])
%     xlim([min(Ct) max(Ct)])
    h1=gca;
    set(h1,'tickdir','out','xminortick','off','yminortick','off');
    set(h1,'YTickLabel',[],'XTickLabel',[]);
    set(h1,'ticklength',0.5*get(h1,'ticklength'));
    set(h1,'fontsize',15);
    text(ax1,25,450,'(a) Image','interpreter','latex','fontsize',15,'Color','w');
    xlabel('$U$ (pix)','interpreter','latex','fontsize',15);
    ylabel('$V$ (pix)','interpreter','latex','fontsize',15);

    ax2 = axes('Position',[0.053 0.08 0.51 0.4]);
    pcolor(ytemp,xtemp,ztemp-Tinfo.tide-.03)
    hold on
    scatter([press.xyz(2,9) press.xyz(2,3)],[press.xyz(1,9) press.xyz(1,3)],40,'r','fill')
    text(ax2,press.xyz(2,9)-0.35,press.xyz(1,9)-0.6,'p06','interpreter','latex','fontsize',15)
    text(ax2,press.xyz(2,3)-0.35,press.xyz(1,3)-0.6,'p11','interpreter','latex','fontsize',15)
    shading flat
    colormap(cmocean('deep'));
    hc = colorbar('Position', [0.575 0.1 0.02 0.3]);
    caxis([-0.2 0.2])
    axis equal
%     grid on
    box on
    shading interp
    ylim([25.5 36])
    xlim([-13.35 13.35])
    h1=gca;
    set(h1,'ydir','reverse')
    set(h1,'tickdir','in','xminortick','off','yminortick','off');
    set(h1,'ticklength',1*get(h1,'ticklength'));
    set(h1,'fontsize',15);
    set(h1,'xtick',[-10:5:10],'xticklabel',{'-10' '-5' '0' '5' '10'});
    text(ax2,13.95,26.25,'$\eta$ (m)','interpreter','latex','fontsize',15);
    xlabel('Alongshore (m)','interpreter','latex','fontsize',15);
    ylabel('Cross-Shore (m)','interpreter','latex','fontsize',15);
    text(ax2,-13,26,'(b) Stereo Reconstruction','interpreter','latex','fontsize',15);
    
    
%     text(ax2,24.7,14.2,['time: ',num2str(round(imageno(i)/8,3),'%.3f'),' sec'],'interpreter','latex','fontsize',15);    
        
    ax3 = axes('Position',[0.68 0.45 0.3 0.3]);
    scatter(Ito,Izoff,4,[0.8 0.8 0.8],'fill')
    hold on
    scatter(Ito,Ifzoff,4,'k','fill')
    scatter(Cto,Czoff,8,'r','fill')
    scatter(Lto,Lzoff,8,'b','fill')
    plot([min(Cto) max(Cto)],[0 0],'LineWidth',0.5,'LineStyle','-','Color',[0.8 0.8 0.8])
    fill([-0.02 -0.02 0.02 0.02],[-0.3 0.5 0.5 -0.3],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
%     scatter(Ito,Izoff,4,[0.8 0.8 0.8],'fill')
    scatter(Ito,Ifzoff,4,'k','fill')
    scatter(Cto,Czoff,8,'r','fill')
    scatter(Lto,Lzoff,8,'b','fill')
    plot(Ito,Izoff,'Color',[0.8 0.8 0.8],'LineWidth',0.5)
    plot(Ito,Ifzoff,'k','LineWidth',0.5)
    plot(Cto,Czoff,'r','LineWidth',0.5)
    plot(Lto,Lzoff,'b','LineWidth',0.5)
%     scatter(Ito((reg*100) + 1),Izoff((reg*100)+1),50,'y','fill')
    scatter(Ito((reg*100) + 1),Ifzoff((reg*100)+1),50,'k','fill')
    scatter(Cto((reg*8)+1),Czoff((reg*8)+1),50,'r','fill')
    scatter(Lto((reg*10)+1),Lzoff((reg*10)+1),50,'b','fill')
%     grid on
    box on
    ylim([-0.3 0.5])
    xlim([min(Cto) max(Cto)])
    h1=gca;
    set(h1,'tickdir','in','xminortick','off','yminortick','off');
    set(h1,'ticklength',1*get(h1,'ticklength'));
    set(h1,'fontsize',15);
    set(h1,'xtick',[-4:1:4],'xticklabel',{''});
    ylabel('$\eta$ (m)','interpreter','latex','fontsize',15);
    h2 = legend('In situ, HR','In situ, TF','Stereo','Lidar','interpreter','latex','fontsize',15)
    set(h2, 'Position', [0.88 0.777 0.10 0.14])
    text(ax3,-3.9, 0.43,['(c) p06, $x=$',num2str(round(press.xyz(1,9),1)),' m'],'interpreter','latex','fontsize',15);
    text(ax3,-4, 0.55,['\textbf{Time Series}'],'interpreter','latex','fontsize',15);
    
    text(ax3,-5,1,['\textbf{Trial}: ',datestr(datenum(Tinfo.cam.tstart(1:end-3),'mm-dd-YYYY-HHMM'),'YYYY-mm-dd HH:MM')],'interpreter','latex','fontsize',15);
    
    text(ax3,-5,0.9,['\textbf{Conditions}: $H_s=',num2str(round(Tinfo.Hs,2)),'$ m, $T_p =',...
        num2str(Tinfo.Tp),'$ s,'],'interpreter','latex','fontsize',15); 
    text(ax3,-2.9,0.8,['$\langle \eta \rangle =',num2str(round(Tinfo.tide,2)),'$ m, $\sigma_{\theta} = ',...
        num2str(Tinfo.spread),' ^{\circ}$'],'interpreter','latex','fontsize',15); 
    text(ax3,-5,0.7,['\textbf{Time Step}: ',datestr(camera.time(i),'HH:MM:SS.FFF')],'interpreter','latex','fontsize',15);

        
    ax4 = axes('Position',[0.68 0.1 0.3 0.3]);
    plot([min(Cto) max(Cto)],[0 0],'LineWidth',1,'LineStyle','-','Color',[.8 .8 .8])
    hold on
    fill([-0.02 -0.02 0.02 0.02],[-0.2 0.2 0.2 -0.2],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
%     scatter(Ito,Izon,4,[0.8 0.8 0.8],'fill')
    scatter(Ito,Ifzon,4,'k','fill')
    scatter(Cto,Czon,8,'r','fill')
    scatter(Lto,Lzon,8,'b','fill')
    plot(Ito,Izon,'Color',[0.8 0.8 0.8],'LineWidth',0.5)
    plot(Ito,Ifzon,'k','LineWidth',0.5)
    plot(Cto,Czon,'r','LineWidth',0.5)
    plot(Lto,Lzon,'b','LineWidth',0.5)
%     scatter(Ito((reg*100) + 1),Izon((reg*100)+1),50,'y','fill')
    scatter(Ito((reg*100) + 1),Ifzon((reg*100)+1),50,'k','fill')
    scatter(Cto((reg*8)+1),Czon((reg*8)+1),50,'r','fill')
    scatter(Lto((reg*10)+1),Lzon((reg*10)+1),50,'b','fill')
%     grid on
    box on
    ylim([-.14 .2])
    xlim([-4 4])
    h1=gca;
    set(h1,'tickdir','in','xminortick','off','yminortick','off');
    set(h1,'ticklength',1*get(h1,'ticklength'));
    set(h1,'fontsize',15);
    set(h1,'xtick',[-4:1:4],'xticklabel',{'-4' '-3' '-2' '-1' '0' '1' '2' '3' '4'});
    xlabel('time (s)','interpreter','latex','fontsize',15);
    ylabel('$\eta$ (m)','interpreter','latex','fontsize',15);
    text(ax4,-3.9, 0.17,['(d) p11, $x=$',num2str(round(press.xyz(1,3),1)),' m'],'interpreter','latex','fontsize',15);
    writeVideo(v,getframe(gcf))
    sname = ['timesync_image_vid_25x5cm_',num2str(count,'%04.f')];
    print([Tinfo.figfolder,sname],'-dpng')
    clf
    clear *temp
end
close(v)

%%

sf = 4;
sname =['timesync_image_vid_25x5cm_mod_imagestereo'];
v = VideoWriter([Tinfo.figfolder,'\',sname,'.avi']);
v.FrameRate=8;
open(v)

figure('units','inches','position',[1 1 16 9],'color','w')
reg = 4;
for i = (1+reg)*8:length(camera.time)-(1+reg)*8 
    count = count+1;
    tifffile    = [tiffpath,Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_frames_',Tinfo.cam.imagerange,'_DEM_',num2str(i),'.tif'];
    % read in point cloud
    [xtemp,ytemp,ztemp,resxtemp,resytemp] = read_tiff(tifffile);
    ytemp = fliplr(ytemp);
    ztemp = flipud(ztemp);
    
    % need to make a 2d for finding nans
    ytemp = repmat(ytemp',1,length(xtemp));
    xtemp = repmat(xtemp,length(ytemp),1);

    %read image
    imagefile = [imagepath,getfield(dir([imagepath,'Movie1_Scene1_c2_',sprintf('%05d',imageno(i)),'_*.jpg']),'name')];
    IM = imread(imagefile);
    
    % plot
    ax1 = axes('Position',[0.07 0.56 0.55 0.4]);
    imagesc(IM)
    hold on
    ylim([350,1800])
%     xlim([min(Ct) max(Ct)])
    h1=gca;
    set(h1,'tickdir','out','xminortick','off','yminortick','off');
    set(h1,'ticklength',0.5*get(h1,'ticklength'));
    set(h1,'fontsize',15);
    text(ax1,25,450,'(a) Image','interpreter','latex','fontsize',15,'Color','w');
    xlabel('$U$ (pix)','interpreter','latex','fontsize',15);
    ylabel('$V$ (pix)','interpreter','latex','fontsize',15);
    text(ax1,0,280,['\textbf{time}: ',num2str(round(count/8)),' s'],'interpreter','latex','fontsize',15);

    ax2 = axes('Position',[0.07 0.08 0.51 0.4]);
    pcolor(ytemp,xtemp,ztemp-Tinfo.tide-.03)
    hold on
    shading flat
    colormap(cmocean('deep'));
    hc = colorbar('Position', [0.595 0.1 0.02 0.3]);
    caxis([-0.2 0.2])
    axis equal
%     grid on
    box on
    shading interp
    ylim([25.5 36])
    xlim([-13.35 13.35])
    h1=gca;
    set(h1,'ydir','reverse')
    set(h1,'tickdir','out','xminortick','on','yminortick','on');
    set(h1,'ticklength',1*get(h1,'ticklength'));
    set(h1,'fontsize',15);
    set(h1,'xtick',[-10:5:10],'xticklabel',{'-10' '-5' '0' '5' '10'});
    text(ax2,13.95,26.25,'$\eta$ (m)','interpreter','latex','fontsize',15);
    xlabel('Alongshore (m)','interpreter','latex','fontsize',15);
    ylabel('Cross-Shore (m)','interpreter','latex','fontsize',15);
    text(ax2,-13,26,'(b) Stereo Reconstruction','interpreter','latex','fontsize',15);       
    writeVideo(v,getframe(gcf))
    sname = ['timesync_image_vid_25x5cm_imagestereo_',num2str(count,'%04.f')];
    print([Tinfo.figfolder,sname],'-dpng')
    clf
    clear *temp
end
close(v)

%% Store mat files
sname = 'szarray_timeseries_allinst_timesync';
eval(['save -v7.3 ',Tinfo.savefolder,sname,' press',' wg',' lidar',' camera'])
clear press wg lidar camera

