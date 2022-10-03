% This code will extract and save time series from the camera dems and
% lidar mat file to compute the directional spectra as an 'array'.
% The timeseries are extracted at the gridded locations.

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\trc_lab_experiment'))
addpath(genpath('E:\codes\insitu'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

sprd    = [20 40 40 40];
Hs = [0.3 0.3 0.25 0.2];
Tp = [3, 2, 2, 2];
tide = [1.07 1 1.07 1.07];
dx = 0.5;
dy = 0.5;
xcam    = [27.5:dx:33];
ycam    = [-10:dy:4];%10
xlid    = [27.5:dx:33];
ylid    = [-10:dy:4];

for op = 1:length(sprd)
    % Trial info
    Tinfo.Hs = Hs(op);
    Tinfo.Tp = Tp(op);
    Tinfo.tide = tide(op);
    Tinfo.spread = sprd(op);

%% STEP 1: Create paths, files and naming

datapath    = 'E:\';
% create file structure
% Tcam        = TRC_camera_info(Tcam);
% [camera.time] = cam_time(Tcam);

Tinfo = trial_files(Tinfo);

% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);
camera.time        = cam_time(Tinfo);

%% STEP 2: Create figure folders

% figure folder
fssubfolder = datestr(date,'yy-mm-dd');
figfolder   = [datapath,'figures\meas_comp\',Tinfo.cam.trialname,'\',Tinfo.cam.trimname,fssubfolder,'\'];

% make figure folders
eval(['!mkdir ',datapath,'figures\meas_comp\',Tinfo.cam.trialname]);
eval(['!mkdir ',datapath,'figures\meas_comp\',Tinfo.cam.trialname,'\',Tinfo.cam.trimname]);
eval(['!mkdir ',figfolder])

%% Data storage location

eval(['!mkdir ',datapath,'data\processed\conditions\'])
eval(['!mkdir ',Tinfo.procpath])
eval(['!mkdir ',Tinfo.procpath,'\',Tinfo.cam.tstart])
eval(['!mkdir ',Tinfo.procpath,'\',Tinfo.cam.tstart,'\',Tinfo.datarange])
savefolder = [Tinfo.procpath,'\',Tinfo.cam.tstart,'\',Tinfo.datarange,'\'];

%% Create grid to extract timeseries

[X,Y] = meshgrid(xcam,ycam);

%% STEP 2: Load Stereo Reconstruction Data

cam             = load([Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_res',num2str((Tinfo.cam.dx)*100),'cm.mat']);
cam.time        = cam_time(Tinfo);

interptseries = 0;
ix = knnsearch(cam.x(1,:)',xcam');
iy = knnsearch(cam.y(:,1),ycam');
z = squeeze(cam.z(iy,ix,:));

if interptseries == 1
    for i = 1:length(ycam)
        for j = 1:length(xcam)
            zint = z(i,j,:);
            nanx = isnan(zint);
            t    = 1:numel(zint);
            zint(nanx) = interp1(t(~nanx), zint(~nanx), t(nanx));
            if isnan(zint(1))
                zint(1) = zintx(2);
            end
            if isnan(zint(end))
                zint(end) = zint(end-1);
                endstr
            end
        end
        z(i,j,:) = zint;
    end
end

% reorganize data for csv
XY = [X(:)'; Y(:)'];
Z = reshape(z,size(XY,2),[])';

fntime = [savefolder,'cam_grid_timeseries_regx',num2str(round(xcam(1))),'-',num2str(round(xcam(end))),'_regy',num2str(round(ycam(1))),'-',num2str(round(ycam(end))),'_dx',num2str(dx*100),'_dy',num2str(dy*100),'.csv'];
fnxyz = [savefolder,'cam_grid_xy_regx',num2str(round(xcam(1))),'-',num2str(round(xcam(end))),'_regy',num2str(round(ycam(1))),'-',num2str(round(ycam(end))),'_dx',num2str(dx*100),'_dy',num2str(dy*100),'.csv'];

timestr = cellstr(datestr(cam.time, 'hh:MM:ss.FFF'));
T1 = array2table(Z);
T1.Properties.RowNames = timestr;
writetable(T1,fntime)

T2 = array2table(XY);
T2.Properties.RowNames = {'x'; 'y'};
writetable(T2,fnxyz)

sname = [savefolder,'cam_grid_timeseries_regx',num2str(round(xcam(1))),'_',num2str(round(xcam(end))),'_regyneg',num2str(abs(round(ycam(1)))),'_',num2str(round(ycam(end))),'_dx',num2str(dx*100),'_dy',num2str(dy*100)];
x = X;
y = Y;
eval(['save -v7.3 ',sname,' x',' y',' z'])
% % end
clear T1 T2 fn* x y z X Y Z

%% STEP 3: Load VeloLidar

[X,Y] = meshgrid(xlid,ylid);
load([datapath,'data\processed\lidar\Velodyne\',Tinfo.lidar.tdate,'\gridded\',Tinfo.lidar.tdate,'_Velodyne-HDL-32-Data_gridded.mat']);

lidar.Hz = 10;

starttemp   = datenum(Tinfo.lidar.tdate,'yyyy-mm-dd-HH-MM-SS');
endtemp     = datenum(Tinfo.lidar.tdate,'yyyy-mm-dd-HH-MM-SS')+(length(squeeze(griddedData.z(1,1,:)))*datenum(0,0,0,0,0,1/lidar.Hz));
timetemp  = starttemp:datenum(0,0,0,0,0,1/lidar.Hz):endtemp;

[temp,istart] = nanmin(abs(camera.time(1)-timetemp));
[temp,iend] = nanmin(abs(camera.time(end)-timetemp));

lidar.time = squeeze(timetemp(istart:iend));
lidar.z    = squeeze(griddedData.z(:,:,istart:iend));
lidar.x     = griddedData.xvec;
lidar.y     = griddedData.yvec;
clear *temp griddedData

interptseries = 0;
ix = knnsearch(lidar.x',xlid');
iy = knnsearch(lidar.y',ylid');
z = squeeze(lidar.z(iy,ix,:));

if interptseries == 1
    for i = 1:length(ylid)
        for j = 1:length(xlid)
            zint = z(i,j,:);
            nanx = isnan(zint);
            t    = 1:numel(zint);
            zint(nanx) = interp1(t(~nanx), zint(~nanx), t(nanx));
            if isnan(zint(1))
                zint(1) = zintx(2);
            end
            if isnan(zint(end))
                zint(end) = zint(end-1);
                endstr
            end
        end
        z(i,j,:) = zint;
    end
end

% reorganize data for csv
XY = [X(:)'; Y(:)'];
Z = reshape(z,size(XY,2),[])';

fntime = [savefolder,'lidar_grid_timeseries_regx',num2str(round(xlid(1))),'-',num2str(round(xlid(end))),'_regy',num2str(round(ylid(1))),'-',num2str(round(ylid(end))),'_dx',num2str(dx*100),'_dy',num2str(dy*100),'.csv'];
fnxyz = [savefolder,'lidar_grid_xy_regx',num2str(round(xlid(1))),'-',num2str(round(xlid(end))),'_regy',num2str(round(ylid(1))),'-',num2str(round(ylid(end))),'_dx',num2str(dx*100),'_dy',num2str(dy*100),'.csv'];

timestr = cellstr(datestr(lidar.time, 'hh:MM:ss.FFF'));
T1 = array2table(Z);
T1.Properties.RowNames = timestr;
writetable(T1,fntime)

T2 = array2table(XY);
T2.Properties.RowNames = {'x'; 'y'};
writetable(T2,fnxyz)

sname = [savefolder,'lidar_grid_timeseries_regx',num2str(round(xlid(1))),'_',num2str(round(xlid(end))),'_regyneg',num2str(abs(round(ylid(1)))),'_',num2str(round(ylid(end))),'_dx',num2str(dx*100),'_dy',num2str(dy*100)];

x = X;
y = Y;
eval(['save -v7.3 ',sname,' x',' y',' z'])

end
