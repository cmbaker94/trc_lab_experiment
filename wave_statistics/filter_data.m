% Compute sea-surface elevation spectra, Hs, and Tp for in situ, stereo,
% and lidar

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\cameras'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

numpts = 3;

% Trial info
Tinfo.Hs        = 0.30;
Tinfo.Tp        = 2;
Tinfo.tide      = 1.07;
Tinfo.spread    = 40;


%% STEP 1: Create paths, files and naming

% % general path and names
datapath    = 'E:\';
% 
Tinfo = trial_files(Tinfo);

% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

%% camera

camfile = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_res',num2str((Tinfo.cam.dx)*100),'cm.mat'];
load(camfile);

z = medfilt3(z);

sname = [camfile(1:end-4),'_filtered.mat'];
eval(['save -v7.3 ',sname,' x',' y',' z'])
clear x y z

%% lidar

lidarfile = [datapath,'data\processed\lidar\Velodyne\',Tinfo.lidar.tdate,'\gridded\',Tinfo.lidar.tdate,'_Velodyne-HDL-32-Data_gridded.mat'];
load(lidarfile);

for i = 1:length(griddedData.xvec)
    z(:,i,:) = medfilt2(squeeze(griddedData.z(:,i,:)));
end

xvec = griddedData.xvec;
yvec = griddedData.yvec;

clear griddedData
griddedData.xvec = xvec;
griddedData.yvec = yvec;
griddedData.z = z;

sname = [lidarfile(1:end-4),'_filtered.mat'];
eval(['save -v7.3 ',sname,' griddedData'])
