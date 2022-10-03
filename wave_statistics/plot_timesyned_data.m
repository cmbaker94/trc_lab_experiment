% This code will extract and save time series from the camera dems and
% lidar mat file to compute the directional spectra as an 'array'.

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('/Users/cmbaker9/Documents/MTOOLS'))
addpath(genpath('/Users/cmbaker9/Documents/Research/Lab_Experiments/codes/insitu'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

sprd = 30;%[0, 20, 40];

for op = 1:length(sprd)
    % Trial info
    Tinfo.Hs = 0.25;
    Tinfo.Tp = 2;
    Tinfo.tide = 1.07;
    Tinfo.spread = sprd(op);
    
     Tinfo = trial_files(Tinfo);
    
    % general path and names
    Tinfo.cam = TRC_camera_info(Tinfo.cam);
    
    % Data and figure storage
    [Tinfo] = wc_comp_store(Tinfo); 

    [camera.time] = cam_time(Tinfo);

% if Tinfo.spread == 0
%     % In situ
%     sz.Tdate = '09-01-2018-2213UTC';
%     is.Tdate = '09-06-2018-1559UTC';
%     % Lidar (LI)
%     lidarfile = '2018-09-01-22-13-48_Velodyne-HDL-32-Data_gridded';
%     % Stereo Reconstructions (SR)
%     Tcam.tstart         = '09-01-2018-2214UTC';                   % time starting collection based on spreadsheet
%     Tcam.tdate          = '09-01-2018-2155UTC';      % trial date and time - format ex: 09-01-2018-2155UTC
%     offset = [-487;0;-70];%[217; 56; 0];
%     % pressure is 2.125 sec ahead of the lidar
%     % camera is 7 sec ahead of the lidar
%     % so relative to the camera, the lidar is 7 sec behind the camera and
%     % the pressure gages are 4.875 sec behind the camera
% elseif Tinfo.spread == 20
%     % In situ
%     sz.Tdate = '08-30-2018-2222UTC';
%     is.Tdate = '09-06-2018-1655UTC';
%     % Lidar (LI)
%     lidarfile = '2018-08-30-22-22-39_Velodyne-HDL-32-Data_gridded';
%     % Stereo Reconstructions (SR)
%     Tcam.tstart  = '08-30-2018-2226UTC';                   % time starting collection based on spreadsheet
%     Tcam.tdate   = '08-30-2018-2216UTC';      % trial date and time - format ex: 09-01-2018-2155UTC 
% elseif Tinfo.spread == 40
%     % In situ
%     sz.Tdate = '08-30-2018-2129UTC'; % 2119+20 min
%     is.Tdate = '09-06-2018-1841UTC';
%     % Lidar (LI)
%     lidarfile = '2018-08-30-21-29-26_Velodyne-HDL-32-Data_gridded';
%     % Stereo Reconstructions (SR)
%     Tcam.tstart  = '08-30-2018-2129UTC';                   % time starting collection based on spreadsheet
%     Tcam.tdate   = '08-30-2018-2119UTC';      % trial date and time - format ex: 09-01-2018-2155UTC 
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
% Tcam.trialname   = [Tcam.camerasys,'-',Tcam.tdate];
% Tcam.imagerange  = [num2str(Tcam.imagestart,'%05.f'),'-',num2str(Tcam.imagestart+(Tcam.numframes-1),'%05.f')];
% Tcam.trimname    = ['frames_',Tcam.imagerange,'/'];  
% Tcam.datafolder  = [datapath,'data/processed/DEM/',Tcam.trialname,'_Scene',Tcam.scene,'/',Tcam.trimname];
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

% % load surf zone insitu data
% tab = readtable([datapath,'data/processed/insitu/surfzone_pressure_gages_xyz.csv']);
% VarNames = tab.Properties.VariableNames;
% for i = 2:length(VarNames)
%     gage_sz(i-1) = str2double(VarNames{i}(6:end));
%     xy_sz(:,i-1)   = tab{1:2,i};
% end
% 
% % load inner shelf insitu data
% tab = readtable([datapath,'data/processed/insitu/innershelf_pressure_gages_xyz.csv']);
% VarNames = tab.Properties.VariableNames;
% for i = 2:length(VarNames)
%     gage_is(i-1) = str2double(VarNames{i}(6:end));
%     xy_is(:,i-1)   = tab{1:2,i};
% end


% load insitu datt
F1          = matfile([Tinfo.datapath,'data/processed/insitu/',Tinfo.sz.tdate,'/',Tinfo.sz.tdate,'-insitu.mat']);
pg          = F1.press;
waveg       = F1.wg;
Tinfo_insitu       = F1.Tinfo;
Hz          = pg.Hz;

starttemp   = datenum(Tinfo_insitu.time(1:end-3));
endtemp     = starttemp+((length(waveg.wg1)-1)*datenum(0,0,0,0,0,1/Hz));
timetemp    = starttemp:datenum(0,0,0,0,0,1/Hz):endtemp;

% find overlapping range for insitu
[temp,istart] = nanmin(abs(camera.time(1)-timetemp));
[temp,iend] = nanmin(abs(camera.time(end)-timetemp));

press.time = timetemp(istart:iend)';
press.name  = fieldnames(pg.xyzd)';

istart = istart+Tinfo.offsets(1);
iend = iend+Tinfo.offsets(1);

% compute bulk statistics
for i = 1:length(press.name)
    % load pressure data
    eval(['press.press(:,i)   = pg.',press.name{i},'(istart:iend);'])
    eval(['press.xyz(:,i)     = pg.xyzd.',press.name{i},'(1,1:3);']) % z-elevation of sensor in tank coordinates
    press.loc{i} = ['p',sprintf('%02d',str2double(press.name{i}(6:end)))];
end

% timestr = cellstr(datestr(press.time, 'hh:MM:ss.FFF'));
% T1 = array2table(press.press);
% T1.Properties.VariableNames = press.loc;
% T1.Properties.RowNames = timestr;
% writetable(T1,[savefolder,'pressure_array_timeseries.csv'])
% 
% T2 = array2table(press.xyz);
% T2.Properties.VariableNames = press.loc;
% T2.Properties.RowNames = {'x'; 'y'; 'z'};
% writetable(T2,[savefolder,'pressure_array_xyz.csv'])

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

% timestr = cellstr(datestr(wg.time, 'hh:MM:ss.FFF'));
% T1 = array2table(wg.z);
% T1.Properties.VariableNames = wg.loc;
% T1.Properties.RowNames = timestr;
% writetable(T1,[savefolder,'wg_array_timeseries.csv'])
% 
% T2 = array2table(wg.xyz);
% T2.Properties.VariableNames = wg.loc;
% T2.Properties.RowNames = {'x'; 'y'; 'z'};
% writetable(T2,[savefolder,'wg_array_xyz.csv'])

clear i* waveg pg T1 T2

%% STEP 2: Load Stereo Reconstruction Data

cam             = load([Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_res',num2str((Tinfo.cam.dx)*100),'cm.mat']);

interptseries = 0;
stepx = 0;%[0, 1, -1, 0, 0];
stepy = 0;%[0, 0, 0, 1, -1];

istart = 1;
camera.time = camera.time(istart:end);

for j = 1:length(stepx)
    for i = 1:size(press.xyz,2)
        [temp,ix] = nanmin(abs(cam.x(1,:)-press.xyz(1,i)));
        [temp,iy] = nanmin(abs(cam.y(:,1)-press.xyz(2,i)));
        ix = ix + stepx(j);
        iy = iy + stepy(j);
        camera.xy(1,i) = cam.x(1,ix);
        camera.xy(2,i) = round(cam.y(iy,1),2);
        x = squeeze(cam.z(iy,ix,istart:end));
%         if interptseries == 1
%             nanx = isnan(x);
%             t    = 1:numel(x);
%             x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
%             if isnan(x(1))
%                 x(1) = x(2);
%             end
%             if isnan(x(end))
%                 x(end) = x(end-1);
%                 endstr
%             end
%         end
        camera.z(:,i)  = x;
        camera.loc{i} = ['c',sprintf('%02d',str2double(press.name{i}(6:end)))];
    end
%     if length(stepx)>1
%         fntime = [savefolder,'cam_array_timeseries_stepx',num2str(stepx(j)),'_stepy',num2str(stepy(j)),'.csv'];
%         fnxyz = [savefolder,'cam_array_xy_stepx',num2str(stepx(j)),'_stepy',num2str(stepy(j)),'.csv'];
%     else
%         fntime = [savefolder,'cam_array_timeseries.csv'];
%         fnxyz = [savefolder,'cam_array_xy.csv'];
%     end
%     timestr = cellstr(datestr(camera.time, 'hh:MM:ss.FFF'));
%     T1 = array2table(camera.z);
%     T1.Properties.VariableNames = camera.loc;
%     T1.Properties.RowNames = timestr;
%     writetable(T1,fntime)
%     
%     T2 = array2table(camera.xy);
%     T2.Properties.VariableNames = camera.loc;
%     T2.Properties.RowNames = {'x'; 'y'};
%     writetable(T2,fnxyz)
end
clear cam T1 T2 fn*

    %% STEP 3: Load VeloLidar
    load([Tinfo.datapath,'data/processed/lidar/Velodyne/',Tinfo.lidar.tdate,'/gridded/',Tinfo.lidar.tdate,'_Velodyne-HDL-32-Data_gridded_filtered.mat']);
    
    lidar.Hz = 10;
    
    starttemp   = datenum(Tinfo.lidar.tdate(1:19),'yyyy-mm-dd-HH-MM-SS');
    endtemp     = datenum(Tinfo.lidar.tdate(1:19),'yyyy-mm-dd-HH-MM-SS')+(length(squeeze(griddedData.z(1,1,:)))*datenum(0,0,0,0,0,1/lidar.Hz));
    timetemp  = starttemp:datenum(0,0,0,0,0,1/lidar.Hz):endtemp;

[temp,istart] = nanmin(abs(camera.time(1)-timetemp));
[temp,iend] = nanmin(abs(camera.time(end)-timetemp));


lidar.time = squeeze(timetemp(istart:iend));
clear *temp

istart = istart+Tinfo.offsets(2);
iend = iend+Tinfo.offsets(2);

interptseries = 0;
stepx = 0;%[0, 1, -1, 0, 0];
stepy = 0;%[0, 0, 0, 1, -1];

for j = 1:length(stepx)
    for i = 1:size(press.xyz,2)
        [temp,ix] = nanmin(abs(griddedData.xvec-press.xyz(1,i)));
        [temp,iy] = nanmin(abs(griddedData.yvec-press.xyz(2,i)));
        ix = ix + stepx(j);
        iy = iy + stepy(j);
        lidar.xy(1,i) = griddedData.xvec(ix);
        lidar.xy(2,i) = griddedData.yvec(iy);
        x = squeeze(griddedData.z(iy,ix,istart:iend));
%         if interptseries == 1
%             nanx = isnan(x);
%             t    = 1:numel(x);
%             x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
%             if isnan(x(1))
%                 x(1) = x(2);
%             end
%             if isnan(x(end))
%                 x(end) = x(end-1);
%             end
%         end
        lidar.z(:,i)  = x;
        lidar.loc{i} = ['l',sprintf('%02d',str2double(press.name{i}(6:end)))];
    end
%     if length(stepx)>1
%         fntime = [savefolder,'lidar_array_timeseries_stepx',num2str(stepx(j)),'_stepy',num2str(stepy(j)),'.csv'];
%         fnxyz = [savefolder,'lidar_array_xy_stepx',num2str(stepx(j)),'_stepy',num2str(stepy(j)),'.csv'];
%     else
%         fntime = [savefolder,'lidar_array_timeseries.csv'];
%         fnxyz = [savefolder,'lidar_array_xy.csv'];
%     end
%     timestr = cellstr(datestr(lidar.time, 'hh:MM:ss.FFF'));
%     T1 = array2table(lidar.z);
%     T1.Properties.VariableNames = lidar.loc;
%     T1.Properties.RowNames = timestr;
%     writetable(T1,fntime)
%     
%     T2 = array2table(lidar.xy);
%     T2.Properties.VariableNames = lidar.loc;
%     T2.Properties.RowNames = {'x'; 'y'};
%     writetable(T2,fnxyz)
end
clear griddedData T1 T2 fn*

%%
 
ptemp = press.press(1:6000,9)/(9.81*1000);
ptemp = ptemp-mean(ptemp);
ctemp = camera.z(1:480,9)-nanmean(camera.z(1:480,9));
ltemp = lidar.z(1:600,9)-nanmean(lidar.z(1:600,9));
 
figure 
plot(press.time(1:6000),ptemp,'k')
hold on
plot(camera.time(1:480),ctemp,'r')
plot(lidar.time(1:600),ltemp,'b')


%% Store mat files
sname = 'szarray_timeseries_allinst';
eval(['save -v7.3 ',savefolder,sname,' press',' wg',' lidar',' camera'])
clear press wg lidar camera
end
