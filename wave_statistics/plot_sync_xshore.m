% This code will extract and save time series from the camera dems and
% lidar mat file to compute the directional spectra as an 'array'.

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\cameras\dem'))
addpath(genpath('E:\code\'))
%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Trial info
    Tinfo.Hs = 0.30;
    Tinfo.Tp = 2;
    Tinfo.tide = 1.07;
    Tinfo.spread = 0;
    yloctran    = 0.1;%-0.068;
    yrange      = 0.1;%0.05
    xROI        = [25 37];              % cross-shore region of interest within the tank
    yROI        = [-13 13];             % alongshore region of interest within the tank
    
    % for camera 2
    UV06 = [1.6868867 0.93791529]*10^3; 
    UV11 = [1.679672566 1.139658546]*10^3;

if Tinfo.spread == 0
    % In situ
    sz.Tdate = '09-01-2018-2213UTC';
    is.Tdate = '09-06-2018-1559UTC';
    % Lidar (LI)
    lidarfile = '2018-09-01-22-13-48_Velodyne-HDL-32-Data_gridded';
    % Stereo Reconstructions (SR)
    Tcam.cam.tstart         = '09-01-2018-2214UTC';                   % time starting collection based on spreadsheet
    Tcam.cam.tdate          = '09-01-2018-2155UTC';      % trial date and time - format ex: 09-01-2018-2155UTC
    offsets             = [-487;-69];% [-588; -70]; % index offset relative to camera [in situ, lidar]
elseif Tinfo.spread == 20
    % In situ
    sz.Tdate = '08-30-2018-2222UTC';
    is.Tdate = '09-06-2018-1655UTC';
    % Lidar (LI)
    lidarfile = '2018-08-30-22-22-39_Velodyne-HDL-32-Data_gridded';
    % Stereo Reconstructions (SR)
    Tcam.cam.tstart  = '08-30-2018-2222UTC';                   % time starting collection based on spreadsheet
    Tcam.cam.tdate   = '08-30-2018-2216UTC';      % trial date and time - format ex: 09-01-2018-2155UTC
elseif Tinfo.spread == 30
    % In situ
    sz.Tdate = '08-29-2018-2255UTC'; % 2119+20 min
    is.Tdate = '09-06-2018-1748UTC';
    % Lidar (LI)
    lidarfile = '2018-08-29-22-55-28_Velodyne-HDL-32-Data_gridded';
    % Stereo Reconstructions (SR)
    Tcam.cam.tstart  = '08-29-2018-2255UTC';                   % time starting collection based on spreadsheet
    Tcam.cam.tdate   = '08-29-2018-2236UTC';      % trial date and time - format ex: 09-01-2018-2155UTC
    offsets      = [3655; 341];%[((36.5/8)*100),  ((34/8)*10)];% [-975+130;-98];offset = [3655; 0; 340];
elseif Tinfo.spread == 40
    % In situ
    sz.Tdate = '08-30-2018-2129UTC'; % 2119+20 min
    is.Tdate = '09-06-2018-1841UTC';
    % Lidar (LI)
    lidarfile = '2018-08-30-21-29-26_Velodyne-HDL-32-Data_gridded';
    % Stereo Reconstructions (SR)
    Tcam.cam.tstart  = '08-30-2018-2129UTC';                   % time starting collection based on spreadsheet
    Tcam.cam.tdate   = '08-30-2018-2119UTC';      % trial date and time - format ex: 09-01-2018-2155UTC 
    toff = -89;%68
    offsets      = [((toff\8)*100)+130, ((toff\8)*10)];% [-975+130;-98];
end

%% STEP 1: Create paths, files and naming


% general path and names
datapath    = 'E:\';

% Stereo Reconstructions
Tcam        = TRC_camera_info(Tcam);

% Bathymetry
bathy       = load([datapath,'data\processed\lidar\Riegal\TRC_bathymetry_estimate_line.mat']);
bathy.h = Tinfo.tide-bathy.h;

% % general path and names
% datapath    = 'E:\';
%
% % Stereo Reconstructions
% Tcam.camerasys   = 'TRM';                     % camera setup - e.g., TRM (offshore) or TRC (onshore)
% Tcam.scene       = '1';                       % scene number of trial - typ 1
% Tcam.imagestart  = 7200;                      % images number of first frame on file
% Tcam.numframes   = 4800;                      % number of frames processed
% Tcam.dx          = 0.5; %0.05
% Tcam.dy          = 0.25; %0.05
% Tcam.regx        = [25:Tcam.dx:37];
% Tcam.regy        = [-13:Tcam.dy:13];
% Tcam.Hz           = 8;
%
% Tcam.trialname   = [Tcam.camerasys,'-',Tcam.tdate];
% Tcam.imagerange  = [num2str(Tcam.imagestart,'%05.f'),'-',num2str(Tcam.imagestart+(Tcam.numframes-1),'%05.f')];
% Tcam.trimname    = ['frames_',Tcam.imagerange,'\'];
% Tcam.datafolder  = [datapath,'data\processed\cameras\',Tcam.trialname,'\',Tcam.trimname];

%% STEP 2: Create figure folders

% figure folder
fssubfolder = datestr(date,'yy-mm-dd');
figfolder   = [datapath,'figures\meas_comp\',Tcam.trialname,'\',Tcam.trimname,fssubfolder,'\'];

% make figure folders
eval(['!mkdir ',datapath,'figures\meas_comp\',Tcam.trialname]);
eval(['!mkdir ',datapath,'figures\meas_comp\',Tcam.trialname,'\',Tcam.trimname]);
eval(['!mkdir ',figfolder])

%% Establish timeseries range:

% starttemp       = datenum(Tcam.tstart(1:end-3),'mm-dd-yyyy-HHMM')+datenum(0,0,0,0,0,Tcam.imagestart\Tcam.Hz);
% endtemp         = datenum(Tcam.tstart(1:end-3),'mm-dd-yyyy-HHMM')+datenum(0,0,0,0,0,(Tcam.imagestart+Tcam.numframes-1)\Tcam.Hz);
% camera.time        = starttemp:datenum(0,0,0,0,0,1\Tcam.Hz):endtemp;
% clear *temp

%% Data storage location

display('storing time series')
Tinfo.comp = ['Hs',num2str(Tinfo.Hs*100),'_Tp',num2str(Tinfo.Tp),'_tide',num2str(Tinfo.tide*100),'_spread',num2str(Tinfo.spread)];
datarange = [datestr(Tcam.timevec(1)); datestr(Tcam.timevec(end))];
datarange = ['time_', datarange(1,13:14), datarange(1,16:17), '-' , datarange(2,13:14), datarange(2,16:17)];
procpath = [datapath,'data\processed\conditions\',Tinfo.comp];

eval(['!mkdir ',datapath,'data\processed\conditions\'])
eval(['!mkdir ',procpath])
eval(['!mkdir ',procpath,'\',Tcam.tstart])
eval(['!mkdir ',procpath,'\',Tcam.tstart,'\',datarange])
savefolder = [procpath,'\',Tcam.tstart,'\',datarange,'\'];

%% Insitu

% load insitu data
F1          = matfile([datapath,'data\processed\insitu\',sz.Tdate,'\',sz.Tdate,'-insitu.mat']);
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
[temp,istart] = nanmin(abs(Tcam.timevec(1)-timetemp));
[temp,iend] = nanmin(abs(Tcam.timevec(end)-timetemp));

press.time = timetemp(istart:iend)';
press.name  = fieldnames(pg.xyzd)';

% offset for time sync
istart = istart+offsets(1);
iend = iend+offsets(1);

% compute bulk statistics
for i = 1:length(press.name)
    % load pressure data
    eval(['press.press(:,i)   = pg.',press.name{i},'(istart:iend);'])
    press.eta(:,i) = press2sse_timeseries(press.press(:,i),Hz,offset,maxfac);
    eval(['press.xyz(:,i)     = pg.xyzd.',press.name{i},'(1,1:3);']) % z-elevation of sensor in tank coordinates
    press.loc{i} = ['p',sprintf('%02d',str2double(press.name{i}(6:end)))];
end

clear pg *temp T1 T2

% % Compute bulk statistics from wave gages
% % get wave gage names
% wg.names	= fieldnames(waveg.xyzd)';
% wg.time     = press.time;
% % compute bulk statistics
% for i = 1:length(wg.names)
%     % load wave gage data
%     eval(['wg.z(:,i)                     = waveg.',wg.names{i},'(istart:iend);'])
%     eval(['wg.xyz(:,i)              = waveg.xyzd.',wg.names{i},'(1,1:3);']) 
%     wg.loc{i} = ['wg',sprintf('%02d',str2double(wg.names{i}(3:end)))];
% end

% clear i* waveg pg T1 T2

%% STEP 2: Load Stereo Reconstruction Data
transectvdem = 1;


if transectvdem == 1
    if round(yloctran)>=0
        psname = [Tcam.datafolder,'dem_transect_y',num2str(abs(round(yloctran))),'m_yavg',num2str(yrange*100),'cm.mat'];
    elseif round(yloctran)<0
        psname = [Tcam.datafolder,'dem_transect_yneg',num2str(abs(round(yloctran))),'m_yavg',num2str(yrange*100),'cm.mat'];
    end
    
    cam         = load(psname);
    
    sse = [28, 31.5];
    eval(['cam.x       = cam.x',num2str(round(abs(yloctran))),'(1,:);'])
    eval(['etabar      = nanmean(cam.z',num2str(round(abs(yloctran))),',1);'])
    [temp,istart] = nanmin(abs(sse(1)-cam.x));
    [temp,iend] = nanmin(abs(sse(end)-cam.x));
    SSE = nanmean(etabar(istart:iend)); 
    eval(['cam.eta         = cam.z',num2str(round(abs(yloctran))),'-SSE;'])
    
elseif transectvdem == 2
    subname = '';
    psname = [Tcam.datafolder,'dem_region_x',num2str(xROI(1)),'to',num2str(xROI(end)),'m_yneg',num2str(abs(yROI(1))),'to',num2str(yROI(end)),'m_res',num2str((yrange)*100),'cm',subname,'.mat'];
    camera         = load(psname);
    [temp,iloc] = nanmin(abs(yloctran-camera.y(:,1)));
    cam.y   = camera.y(iloc,1);
    cam.x   = camera.x(1,:)';
    cam.z   = (squeeze(camera.z(iloc,:,:)))';
    clear camera
    
    sse = [28, 31.5];
    etabar      = nanmean(cam.z,1);
    [temp,istart] = nanmin(abs(sse(1)-cam.x));
    [temp,iend] = nanmin(abs(sse(end)-cam.x));
    SSE = nanmean(etabar(istart:iend));  
    cam.eta         = cam.z-SSE;
end


%% STEP 3: Load VeloLidar
load([datapath,'data\processed\lidar\Velodyne\',lidarfile(1:19),'\gridded\',lidarfile,'.mat']);

lidar.Hz = 10;

starttemp   = datenum(lidarfile(1:19),'yyyy-mm-dd-HH-MM-SS');
endtemp     = datenum(lidarfile(1:19),'yyyy-mm-dd-HH-MM-SS')+(length(squeeze(griddedData.z(1,1,:)))*datenum(0,0,0,0,0,1/lidar.Hz));
timetemp  = starttemp:datenum(0,0,0,0,0,1/lidar.Hz):endtemp;

[temp,istart] = nanmin(abs(Tcam.timevec(1)-timetemp));
[temp,iend] = nanmin(abs(Tcam.timevec(end)-timetemp));

lidar.time = squeeze(timetemp(istart:iend));
clear *temp

% offset for time sync
istart = istart+offsets(2);
iend = iend+offsets(2);


% [temp,ix] = nanmin(abs(griddedData.xvec-press.xyz(1,i)));
[temp,iy] = nanmin(abs(griddedData.yvec-yloctran));
lidar.x = griddedData.xvec;
lidar.y = griddedData.yvec(iy);
lidar.z = squeeze(griddedData.z(iy,:,istart:iend))';

etabar      = nanmean(lidar.z,1);

[temp,istart] = nanmin(abs(sse(1)-lidar.x));
[temp,iend] = nanmin(abs(sse(end)-lidar.x));
SSE = nanmean(etabar(istart:iend));

lidar.eta         = lidar.z-SSE;

clear griddedData

%% plot
imagepath = ['D:\TRC_Fall_Experiment\',Tcam.camerasys,'-',Tcam.tdate,'\',Tcam.camerasys,'-',Tcam.tdate,'_Scene1_JPEG\c2\'];
imageno = Tcam.imagestart:Tcam.imagestart+Tcam.numframes;
% presson = nanmean(press.press(:,3));%+0.05;
% pressoff = nanmean(press.press(:,9));%+0.05;
pressx(1)    = press.xyz(1,3);
pressx(2)    = press.xyz(1,9);
presseta(:,1)    = press.eta(:,3);
presseta(:,2)    = press.eta(:,9);

date2sec = (24*3600);
count = 0;

figure('units','inches','position',[1 1 13 9],'color','w')
reg = 4;
for i = 1:length(Tcam.timevec)
    count = count+1;
%     tifffile    = [tiffpath,Tcam.camerasys,'-',Tcam.tdate,'_frames_',Tcam.imagerange,'_DEM_',num2str(i),'.tif'];
%     % read in point cloud
%     [xtemp,ytemp,ztemp,resxtemp,resytemp] = read_tiff(tifffile);
%     ytemp = fliplr(ytemp);
%     ztemp = flipud(ztemp);
%     
%     % need to make a 2d for finding nans
%     ytemp = repmat(ytemp',1,length(xtemp));
%     xtemp = repmat(xtemp,length(ytemp),1);
% 
%     %read image
    imagefile = [imagepath,getfield(dir([imagepath,'Movie1_Scene1_c2_',sprintf('%05d',imageno(i)),'_*.jpg']),'name')];
    IM = imread(imagefile);
    
    % timeseries 
    [temp,iI] = nanmin(abs(Tcam.timevec(i)-press.time));
    [temp,iL] = nanmin(abs(Tcam.timevec(i)-lidar.time));
    
%     It = press.time(iI-(reg*100):iI+reg*100);
%     Ito = (It-camera.time(i))*24*3600;
%     Izon = (press.press(iI-reg*100:iI+reg*100,3)-h.presson)\(g*rho);
%     Izoff = (press.press(iI-reg*100:iI+reg*100,9)-h.pressoff)\(g*rho);
%     Ifzon = press.eta(iI-reg*100:iI+reg*100,3);
%     Ifzoff = press.eta(iI-reg*100:iI+reg*100,9);
    
%     Ct = camera.time(i-reg*8:i+reg*8);
%     Cto = (Ct-camera.time(i))*24*3600;
%     Czon = camera.z(i-reg*8:i+reg*8,3)-h.cameraon;
%     Czoff = camera.z(i-reg*8:i+reg*8,9)-h.cameraoff;
%     
%     Lt = lidar.time(iL-reg*10:iL+reg*10);
%     Lto = (Lt-camera.time(i))*24*3600;
%     Lzon = lidar.z(iL-reg*10:iL+reg*10,3)-h.lidaron;
%     Lzoff = lidar.z(iL-reg*10:iL+reg*10,9)-h.lidaroff;

    % plot
    ax1 = axes('Position',[0.073 0.55 0.72 0.4]);
    imagesc(IM)
    hold on
    scatter([UV06(1) UV11(1)],[UV06(2) UV11(2)],40,'r','fill')
    text(ax1,UV06(1)-55,UV06(2)-70,'p06','Color','w','interpreter','latex','fontsize',15)
    text(ax1,UV11(1)-55,UV11(2)-70,'p11','Color','w','interpreter','latex','fontsize',15)
    ylim([350,1800])
%     xlim([min(Ct) max(Ct)])
    h1=gca;
    set(h1,'tickdir','out','xminortick','off','yminortick','off');
    set(h1,'ticklength',0.5*get(h1,'ticklength'));
    set(h1,'fontsize',15);
    text(ax1,25,270,'(a) Image','interpreter','latex','fontsize',15,'Color','k');
    xlabel('$U$ (pix)','interpreter','latex','fontsize',15);
    ylabel('$V$ (pix)','interpreter','latex','fontsize',15);
%      text(ax1,25,450,['\textbf{Time Series}'],'interpreter','latex','fontsize',15);
    text('Parent',ax1,'FontSize',15,'Interpreter','latex','String',['\textbf{Trial}:',sprintf('\n'),...
        datestr(datenum(Tcam.tstart(1:end-3),'mm-dd-YYYY-HHMM'),'YYYY-mm-dd HH:MM')],'Position',[3610 454 0])
%     text(ax1,25,450,['\textbf{Trial}: ',datestr(datenum(Tcam.tstart(1:end-3),'mm-dd-YYYY-HHMM'),'YYYY-mm-dd HH:MM')],'interpreter','latex','fontsize',15);
    text('Parent',ax1,'FontSize',15,'Interpreter','latex','String',['\textbf{Conditions}:',sprintf('\n'),...
        '$H_s=',num2str(round(Tinfo.Hs,2)),'$ m,',sprintf('\n'),...
        '$T_p =',num2str(Tinfo.Tp),'2$ s,',sprintf('\n'),...
        '$\langle \eta \rangle =',num2str(round(Tinfo.tide,2)),'$ m,',sprintf('\n'),...
        '$\sigma_{\theta} = ',num2str(Tinfo.spread),' ^{\circ}$'],'Position',[3610 890 0]);
    text('Parent',ax1,'FontSize',15,'Interpreter','latex','String',['\textbf{Time Step}:',sprintf('\n'),...
        datestr(Tcam.timevec(i),'HH:MM:SS.FFF')],'Position',[3610 1370 0]);
%     text(ax1,25,450,['\textbf{Conditions}: $H_s=',num2str(round(Tinfo.Hs,2)),'$ m, $T_p =',...
%         num2str(Tinfo.Tp),'$ s,'],'interpreter','latex','fontsize',15); 
%     text(ax1,25,450,['$\langle \eta \rangle =',num2str(round(Tinfo.tide,2)),'$ m, $\sigma_{\theta} = ',...
%         num2str(Tinfo.spread),' ^{\circ}$'],'interpreter','latex','fontsize',15); 
%     text(ax1,25,450,['\textbf{Time Step}: ',datestr(Tcam.timevec(i),'HH:MM:SS.FFF')],'interpreter','latex','fontsize',15);

%     ax2 = axes('Position',[0.053 0.08 0.51 0.4]);
%     pcolor(ytemp,xtemp,ztemp-Tinfo.tide-.03)
%     hold on
%     scatter([press.xyz(2,9) press.xyz(2,3)],[press.xyz(1,9) press.xyz(1,3)],40,'r','fill')
%     text(ax2,press.xyz(2,9)-0.35,press.xyz(1,9)-0.6,'p06','interpreter','latex','fontsize',15)
%     text(ax2,press.xyz(2,3)-0.35,press.xyz(1,3)-0.6,'p11','interpreter','latex','fontsize',15)
%     shading flat
%     colormap(cmocean('deep'));
%     hc = colorbar('Position', [0.575 0.1 0.02 0.3]);
%     caxis([-0.2 0.2])
%     axis equal
% %     grid on
%     box on
%     shading interp
%     ylim([25.5 36])
%     xlim([-13.35 13.35])
%     h1=gca;
%     set(h1,'ydir','reverse')
%     set(h1,'tickdir','out','xminortick','on','yminortick','on');
%     set(h1,'ticklength',1*get(h1,'ticklength'));
%     set(h1,'fontsize',15);
%     set(h1,'xtick',[-10:5:10],'xticklabel',{'-10' '-5' '0' '5' '10'});
%     text(ax2,13.95,26.25,'$\eta$ (m)','interpreter','latex','fontsize',15);
%     xlabel('Alongshore (m)','interpreter','latex','fontsize',15);
%     ylabel('Cross-Shore (m)','interpreter','latex','fontsize',15);
%     text(ax2,-13,26,'(b) Stereo Reconstruction','interpreter','latex','fontsize',15);
    
    
%     text(ax2,24.7,14.2,['0.078173611111111 0.110379584681769 0.9 0.322916666666667time: ',num2str(round(imageno(i)\8,3),'%.3f'),' sec'],'interpreter','latex','fontsize',15);    

    ax3 = axes('Position',[0.073 0.08 0.87 0.32]);
    scatter(pressx,presseta(iI,:),10,'k','fill')
    hold on
    scatter(cam.x,cam.eta(i,:),8,'r','fill')
    scatter(lidar.x,lidar.eta(iL,:),8,'b','fill')
    plot(bathy.xp,-bathy.h,'Color',[0.3 0.3 0.3],'LineWidth',1.5);
%     plot([min(Cto) max(Cto)],[0 0],'LineWidth',0.5,'LineStyle','-','Color',[0.8 0.8 0.8])
%     fill([-0.02 -0.02 0.02 0.02],[-0.3 0.5 0.5 -0.3],[239 242 157]\256,'LineStyle','none','EdgeColor','none')
    plot([min(cam.x) max(cam.x)],[0 0],'LineWidth',0.5,'LineStyle','-','Color',[0.8 0.8 0.8])
    scatter(cam.x,cam.eta(i,:),8,'r','fill')
    scatter(lidar.x,lidar.eta(iL,:),8,'b','fill')
    plot(cam.x,cam.eta(i,:),'r','LineWidth',0.5)
    plot(lidar.x,lidar.eta(iL,:),'b','LineWidth',0.5)
    scatter(pressx,presseta(iI,:),50,'k','fill')
    
    bathyfill = fill([bathy.xp; flipud(bathy.xp)],[-bathy.h; -2*ones(size(bathy.h))],[0.7 0.7 0.7])
%     plot([18 40],[0 0],'Color','k','LineStyle','--','LineWidth',1);
    alpha(bathyfill,0.6)
%     grid on
    box on
    ylim([-0.5 0.5])
    xlim([24 36])
    h1=gca;
    set(h1,'tickdir','out','xminortick','on','yminortick','on');
    set(h1,'ticklength',0.7*get(h1,'ticklength'));
    set(h1,'fontsize',15);
    set(h1,'ytick',[-0.6:0.2:0.6],'yticklabel',{'-0.6' '-0.4' '-0.2' '0' '0.2' '0.4' '0.6'});
%     set(h1,'xtick',[-4:1:4],'xticklabel',{''});
    xlabel('$x$ (m)','interpreter','latex','fontsize',15);
    ylabel('$\eta$ (m)','interpreter','latex','fontsize',15);
    h2 = legend('In situ','Stereo','Lidar','h','interpreter','latex','fontsize',15)
    set(h2, 'Position', [0.84 0.45 0.1 0.16])
    text(ax3,24.1, 0.6,['(b) Cross-Shore Transect'],'interpreter','latex','fontsize',15);
    tvecplot = ((Tcam.timevec(i)-lidar.time(iL))*date2sec)*1000;
    display(tvecplot)
    text(ax3,24.1, 0.43,['$t_{\mathrm{cam}} -  t_{\mathrm{lidar}}$ = ',num2str(round(tvecplot)),' ms'],'interpreter','latex','fontsize',13);
       
    sname = ['timesync_image_vid_xshore_',num2str(count,'%04.f')];
    print([figfolder,sname],'-dpng')
    clf
    clear *temp
end

%% Store mat files
% sname = 'szarray_timeseries_allinst_timesync';
% eval(['save -v7.3 ',savefolder,sname,' press',' wg',' lidar',' camera'])
% clear press wg lidar camera

