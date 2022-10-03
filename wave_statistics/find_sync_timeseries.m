% This code will extract and save time series from the camera dems and
% lidar mat file to compute the directional spectra as an 'array'.

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('/Users/cmbaker9/Documents/MTOOLS'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Trial info
    Tinfo.Hs = 0.30;
    Tinfo.Tp = 2;
    Tinfo.tide = 1.07;
    Tinfo.spread = 40;
    
    % for camera 2
    UV06 = [1.6868867 0.93791529]*10^3; 
    UV11 = [1.679672566 1.139658546]*10^3;
    
    Tinfo = trial_files(Tinfo);

% Stereo Reconstructions
Tinfo.cam = TRC_camera_info(Tinfo.cam);
camera.time        = cam_time(Tinfo);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

% if Tinfo.spread == 0
%     % In situ
%     sz.Tdate = '09-01-2018-2213UTC';
%     is.Tdate = '09-06-2018-1559UTC';
%     % Lidar (LI)
%     lidarfile = '2018-09-01-22-13-48_Velodyne-HDL-32-Data_gridded';
%     % Stereo Reconstructions (SR)
%     Tcam.tstart         = '09-01-2018-2214UTC';                   % time starting collection based on spreadsheet
%     Tcam.tdate          = '09-01-2018-2155UTC';      % trial date and time - format ex: 09-01-2018-2155UTC
%     offsets             = [-487;-69];% [-588; -70]; % index offset relative to camera [in situ, lidar]
% elseif Tinfo.spread == 20
%     % In situ
%     sz.Tdate = '08-30-2018-2222UTC';
%     is.Tdate = '09-06-2018-1655UTC';
%     % Lidar (LI)
%     lidarfile = '2018-08-30-22-22-39_Velodyne-HDL-32-Data_gridded';
%     % Stereo Reconstructions (SR)
%     Tcam.tstart  = '08-30-2018-2222UTC';                   % time starting collection based on spreadsheet
%     Tcam.tdate   = '08-30-2018-2216UTC';      % trial date and time - format ex: 09-01-2018-2155UTC 
%     offsets      = [291, 0];
% elseif Tinfo.spread == 40
%     % In situ
%     sz.Tdate = '08-30-2018-2129UTC'; % 2119+20 min
%     is.Tdate = '09-06-2018-1841UTC';
%     % Lidar (LI)
%     lidarfile = '2018-08-30-21-29-26_Velodyne-HDL-32-Data_gridded';
%     % Stereo Reconstructions (SR)
%     Tcam.tstart  = '08-30-2018-2129UTC';                   % time starting collection based on spreadsheet
%     Tcam.tdate   = '08-30-2018-2119UTC';      % trial date and time - format ex: 09-01-2018-2155UTC 
%     offsets      = [130, 0];
% end
% 
% %% STEP 1: Create paths, files and naming
% 
% % general path and names
% datapath    = '/Users/cmbaker9/Documents/Research/Lab_Experiments/';
% 
% % Stereo Reconstructions
% Tcam.camerasys   = 'TRM';                     % camera setup - e.g., TRM (offshore) or TRC (onshore)
% Tcam.scene       = '1';                       % scene number of trial - typ 1
% Tcam.imagestart  = 7200;                      % images number of first frame on file 
% Tcam.numframes   = 4800;                      % number of frames processed
% Tcam.dx          = 0.05;
% Tcam.dy          = 0.05;
% Tcam.regx        = [25:Tcam.dx:37];
% Tcam.regy        = [-13:Tcam.dy:13];
% Tcam.Hz           = 8;
% 
% Tcam.trialname   = [Tcam.camerasys,'-',Tcam.tdate,'_Scene1'];
% Tcam.imagerange  = [num2str(Tcam.imagestart,'%05.f'),'-',num2str(Tcam.imagestart+(Tcam.numframes-1),'%05.f')];
% Tcam.trimname    = ['frames_',Tcam.imagerange,'/'];  
% Tcam.datafolder  = [datapath,'data/processed/cameras/',Tcam.trialname,'/',Tcam.trimname];
% 
% %% STEP 2: Create figure folders
% 
% % figure folder
% fssubfolder = datestr(date,'yy-mm-dd');
% figfolder   = [datapath,'figures/meas_comp/',Tcam.trialname,'/',Tcam.trimname,fssubfolder,'/'];
% 
% % make figure folders
% eval(['!mkdir ',datapath,'figures/meas_comp/',Tcam.trialname]);
% eval(['!mkdir ',datapath,'figures/meas_comp/',Tcam.trialname,'/',Tcam.trimname]);
% eval(['!mkdir ',figfolder])
% 
% %% Establish timeseries range:
% 
% starttemp       = datenum(Tcam.tstart(1:end-3),'mm-dd-yyyy-HHMM')+datenum(0,0,0,0,0,Tcam.imagestart/Tcam.Hz);
% endtemp         = datenum(Tcam.tstart(1:end-3),'mm-dd-yyyy-HHMM')+datenum(0,0,0,0,0,(Tcam.imagestart+Tcam.numframes-1)/Tcam.Hz);
% camera.time        = starttemp:datenum(0,0,0,0,0,1/Tcam.Hz):endtemp;
% clear *temp
% 
% %% Data storage location
% 
% display('storing time series')
% Tinfo.comp = ['Hs',num2str(Tinfo.Hs*100),'_Tp',num2str(Tinfo.Tp),'_tide',num2str(Tinfo.tide*100),'_spread',num2str(Tinfo.spread)];
% datarange = [datestr(camera.time(1)); datestr(camera.time(end))];
% datarange = ['time_', datarange(1,13:14), datarange(1,16:17), '-' , datarange(2,13:14), datarange(2,16:17)];
% procpath = [datapath,'data/processed/conditions/',Tinfo.comp];
% 
% eval(['!mkdir ',datapath,'data/processed/conditions/'])
% eval(['!mkdir ',procpath])
% eval(['!mkdir ',procpath,'/',Tcam.tstart])
% eval(['!mkdir ',procpath,'/',Tcam.tstart,'/',datarange])
% savefolder = [procpath,'/',Tcam.tstart,'/',datarange,'/'];

%% Insitu

% load insitu data
F1          = matfile([datapath,'data/processed/insitu/',sz.Tdate,'/',sz.Tdate,'-insitu.mat']);
pg          = F1.press;
waveg       = F1.wg;
Tinfo_insitu       = F1.Tinfo;
Hz          = pg.Hz;
maxfac = 1.2;
offset  = 0.05;
rho = 1000;
g = 9.81;

starttemp   = datenum(Tinfo_insitu.time(1:end-3));
endtemp     = starttemp+((length(waveg.wg1)-1)*datenum(0,0,0,0,0,1/Hz));
timetemp    = starttemp:datenum(0,0,0,0,0,1/Hz):endtemp;
PRESSTIME = timetemp;

% find overlapping range for insitu
% [temp,istart] = nanmin(abs(camera.time(1)-timetemp));
% [temp,iend] = nanmin(abs(camera.time(end)-timetemp));

% press.time = timetemp(istart:iend)';
istart = offsets(1);
press.time = timetemp(1:end-istart+1)';
press.name  = fieldnames(pg.xyzd)';

% offset for time sync
% istart = istart+offsets(1);
% iend = iend+offsets(1);


istart = 1; % TEMPPPPPPPPPP


% compute bulk statistics
for i = 1:length(press.name)
    % load pressure data
    eval(['press.press(:,i)   = pg.',press.name{i},'(istart:end);'])
    %     press.eta(:,i) = press2sse_timeseries(press.press(:,i),Hz,offset,maxfac);
    press.eta(:,i) = press.press(:,i)/(rho*g);
    press.eta(:,i) = press.eta(:,i)-nanmean(press.eta(:,i));
    eval(['press.xyz(:,i)     = pg.xyzd.',press.name{i},'(1,1:3);']) % z-elevation of sensor in tank coordinates
    press.loc{i} = ['p',sprintf('%02d',str2double(press.name{i}(6:end)))];
end
plot(press.eta(:,5)); xlim([0,10000])
clear pg *temp T1 T2

% Compute bulk statistics from wave gages
% get wave gage names
wg.names	= fieldnames(waveg.xyzd)';
wg.time     = press.time;
% compute bulk statistics
for i = 1:length(wg.names)
    % load wave gage data
    eval(['wg.z(:,i)                     = waveg.',wg.names{i},'(istart:end);'])
    eval(['wg.xyz(:,i)              = waveg.xyzd.',wg.names{i},'(1,1:3);']) 
    wg.loc{i} = ['wg',sprintf('%02d',str2double(wg.names{i}(3:end)))];
end


% timestr = cellstr(datestr(press.time, 'hh:MM:ss.FFF'));
% T1 = array2table(press.press);
% T1.Properties.VariableNames = press.loc;
% T1.Properties.RowNames = timestr;
% writetable(T1,[savefolder,'pressure_array_timeseries_full.csv'])
% 
% timestr = cellstr(datestr(wg.time, 'hh:MM:ss.FFF'));
% T1 = array2table(wg.z);
% T1.Properties.VariableNames = wg.loc;
% T1.Properties.RowNames = timestr;
% writetable(T1,[savefolder,'wg_array_timeseries_full.csv'])

clear i* waveg pg T1 T2

%% STEP 2: Load Stereo Reconstruction Data

cam             = load([Tcam.datafolder,'dem_region_x',num2str(Tcam.regx(1)),'to',num2str(Tcam.regx(end)),'m_yneg',num2str(abs(Tcam.regy(1))),'to',num2str(Tcam.regy(end)),'m_res',num2str((Tcam.dx)*100),'cm.mat']);

for i = 1:size(press.xyz,2)
    [temp,ix] = nanmin(abs(cam.x(1,:)-press.xyz(1,i)));
    [temp,iy] = nanmin(abs(cam.y(:,1)-press.xyz(2,i)));
    camera.xy(1,i) = cam.x(1,ix);
    camera.xy(2,i) = round(cam.y(iy,1),2);
    x = squeeze(cam.z(iy,ix,:));
    camera.z(:,i)  = x;
    camera.loc{i} = ['c',sprintf('%02d',str2double(press.name{i}(6:end)))];
end

clear cam T1 T2 fn* cam

%% STEP 3: Load VeloLidar
load([datapath,'data/processed/lidar/Velodyne/',lidarfile,'.mat']);

lidar.Hz = 10;

starttemp   = datenum(lidarfile(1:19),'yyyy-mm-dd-HH-MM-SS');
endtemp     = datenum(lidarfile(1:19),'yyyy-mm-dd-HH-MM-SS')+(length(squeeze(griddedData.z(1,1,:)))*datenum(0,0,0,0,0,1/lidar.Hz));
timetemp  = starttemp:datenum(0,0,0,0,0,1/lidar.Hz):endtemp;
LIDARTIME = timetemp;

[temp,istart] = nanmin(abs(press.time(1)-timetemp));
% [temp,iend] = nanmin(abs(camera.time(end)-timetemp));
istart = istart+offsets(2);
lidar.time = squeeze(timetemp(istart:end-1));
clear *temp

% offset for time sync
% istart = offsets(2);


istart = 1; % TEMPPPPPPPPPP


for i = 1:size(press.xyz,2)
    [temp,ix] = nanmin(abs(griddedData.xvec-press.xyz(1,i)));
    [temp,iy] = nanmin(abs(griddedData.yvec-press.xyz(2,i)));
    lidar.xy(1,i) = griddedData.xvec(ix);
    lidar.xy(2,i) = griddedData.yvec(iy);
    x = squeeze(griddedData.z(iy,ix,istart:end));
    lidar.z(:,i)  = x;
    lidar.loc{i} = ['l',sprintf('%02d',str2double(press.name{i}(6:end)))];
end
clear griddedData T1 T2 fn*

% fntime = [savefolder,'lidar_array_timeseries_full.csv'];
% 
% timestr = cellstr(datestr(lidar.time, 'hh:MM:ss.FFF'));
% T1 = array2table(lidar.z);
% T1.Properties.VariableNames = lidar.loc;
% T1.Properties.RowNames = timestr;
% writetable(T1,fntime)

%% testing time sync
lidarseries = detrend(fillmissing(lidar.z(:,5),'linear'));
figure
plot(press.time(2000:3500),press.eta(2000:3500,5))
hold on
% plot(lidar.time(200:350),lidar.z(200:350,5)-nanmean(lidar.z(:,5)))
plot(lidar.time(200:350),lidarseries(200:350))

%% cross-correlation
lidarseries = detrend(fillmissing(lidar.z(:,5),'linear'));
cameraseries = detrend(fillmissing(camera.z(:,5),'linear'));
ltime  = lidar.time(1):datenum(0,0,0,0,0,1/8):lidar.time(end);

lin = timeseries(lidarseries,lidar.time(1:end));
lout = resample(lin,ltime);
ldata = squeeze(lout.data);


%% plot
h.presson = nanmean(press.press(:,3));%+0.05;
h.pressoff = nanmean(press.press(:,9));%+0.05;
h.cameraon = nanmean(camera.z(:,3));
h.cameraoff = nanmean(camera.z(:,9)); 
% h.camera = nanmean(nanmean(camera.z));
h.lidaron = nanmean(lidar.z(:,3));
h.lidaroff = nanmean(lidar.z(:,9)); 
% h.lidar = nanmean(nanmean(lidar.z));
g = 9.81;
rho = 1000;
count = 0;

figure('units','inches','position',[1 1 20 9],'color','w')
reg = 30;
% MVCAM = -24.38-3+30.4250;
% itime = camera.time(1500);
% [temp,iC] = nanmin(abs(camera.time(1500)-camera.time));
iC = 2400;
Ct = camera.time(iC-reg*8:iC+reg*8);
Cto = (Ct-camera.time(iC))*24*3600;
Czon = camera.z(iC-reg*8:iC+reg*8,3)-h.cameraon;
Czoff = camera.z(iC-reg*8:iC+reg*8,9)-h.cameraoff;



for i = (1+reg)*8:length(camera.time)-(1+reg)*8;
    count = count+1;
    
    itime = camera.time(i);
    % timeseries 
    [temp,iI] = nanmin(abs(itime-press.time));
    [temp,iL] = nanmin(abs(itime-lidar.time));
    
    It = press.time(iI-(reg*100):iI+reg*100);
    Ito = (It-itime)*24*3600;
    Izon = (press.press(iI-reg*100:iI+reg*100,3)-h.presson)/(g*rho);
    Izoff = (press.press(iI-reg*100:iI+reg*100,9)-h.pressoff)/(g*rho);
    Ifzon = press.eta(iI-reg*100:iI+reg*100,3);
    Ifzoff = press.eta(iI-reg*100:iI+reg*100,9);

    Lt = lidar.time(iL-reg*10:iL+reg*10);
    Lto = (Lt-itime)*24*3600;
    Lzon = lidar.z(iL-reg*10:iL+reg*10,3)-h.lidaron;
    Lzoff = lidar.z(iL-reg*10:iL+reg*10,9)-h.lidaroff;
    
    
%     text(ax2,24.7,14.2,['time: ',num2str(round(imageno(i)/8,3),'%.3f'),' sec'],'interpreter','latex','fontsize',15);    
        
    ax3 = subplot(2,1,1); %axes('Position',[0.68 0.45 0.3 0.3]);
    scatter(Ito,Izoff,4,[0.8 0.8 0.8],'fill')
    hold on
    scatter(Ito,Ifzoff,4,'k','fill')
    scatter(Cto,Czoff,8,'r','fill')
    scatter(Lto,Lzoff,8,'b','fill')
    plot([min(Cto) max(Cto)],[0 0],'LineWidth',0.5,'LineStyle','-','Color',[0.8 0.8 0.8])
    fill([-0.02 -0.02 0.02 0.02],[-0.3 0.5 0.5 -0.3],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
    scatter(Ito,Izoff,4,[0.8 0.8 0.8],'fill')
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
    ylim([-0.2 0.4])
    xlim([min(Cto) max(Cto)])
    h1=gca;
    set(h1,'tickdir','out','xminortick','off','yminortick','off');
    set(h1,'ticklength',1*get(h1,'ticklength'));
    set(h1,'fontsize',15);
%     set(h1,'xtick',[-8:1:8],'xticklabel',{''});
    ylabel('$/eta$ (m)','interpreter','latex','fontsize',15);
    h2 = legend('In situ, HR','In situ, TF','Stereo','Lidar','interpreter','latex','fontsize',15)
    set(h2, 'Position', [0.88 0.777 0.10 0.14])
    text(ax3,-3.9, 0.43,['(c) p06, $x=$',num2str(round(press.xyz(1,9),1)),' m'],'interpreter','latex','fontsize',15);
%     text(ax3,-4, 0.55,['/textbf{Time Series}'],'interpreter','latex','fontsize',15);
    
    text(ax3,-30,0.55,['/textbf{Trial}: ',datestr(datenum(Tcam.tstart(1:end-3),'mm-dd-YYYY-HHMM'),'YYYY-mm-dd HH:MM')],'interpreter','latex','fontsize',15);
    
    text(ax3,-20,0.55,['/textbf{Conditions}: $H_s=',num2str(round(Tinfo.Hs,2)),'$ m, $T_p =',...
        num2str(Tinfo.Tp),'$ s,','$/langle /eta /rangle =',num2str(round(Tinfo.tide,2)),'$ m, $/sigma_{/theta} = ',...
        num2str(Tinfo.spread),' ^{/circ}$'],'interpreter','latex','fontsize',15); 
    text(ax3,5,0.55,['/textbf{Time Step}: ',datestr(camera.time(iC),'HH:MM:SS.FFF')],'interpreter','latex','fontsize',15);
    text(ax3,20,0.55,['/textbf{Offset}: ',num2str((iC-i)/8), ' sec'],'interpreter','latex','fontsize',15);

        
    ax4 = subplot(2,1,2); %axes('Position',[0.68 0.1 0.3 0.3]);
    plot([min(Cto) max(Cto)],[0 0],'LineWidth',1,'LineStyle','-','Color',[.8 .8 .8])
    hold on
    fill([-0.02 -0.02 0.02 0.02],[-0.2 0.2 0.2 -0.2],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
    scatter(Ito,Izon,4,[0.8 0.8 0.8],'fill')
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
    xlim([min(Cto) max(Cto)])
    h1=gca;
    set(h1,'tickdir','out','xminortick','off','yminortick','off');
    set(h1,'ticklength',1*get(h1,'ticklength'));
    set(h1,'fontsize',15);
%     set(h1,'xtick',[-8:2:8],'xticklabel',{'-8' '-6' '-4' '-2' '0' '2' '4' '6' '8'});
    xlabel('time (s)','interpreter','latex','fontsize',15);
    ylabel('$/eta$ (m)','interpreter','latex','fontsize',15);
    text(ax4,-3.9, 0.17,['(d) p11, $x=$',num2str(round(press.xyz(1,3),1)),' m'],'interpreter','latex','fontsize',15);
    
    sname = ['timesync_trial_',num2str(count,'%04.f')];
    print([figfolder,sname],'-dpng')
    clf
    clear *temp
end

%% Store mat files
sname = 'szarray_timeseries_allinst_timesync';
eval(['save -v7.3 ',savefolder,sname,' press',' wg',' lidar',' camera'])
clear press wg lidar camera

