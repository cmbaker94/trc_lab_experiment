% This code will extract and save time series from the camera dems and
% lidar mat file to compute the directional spectra as an 'array'.

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('/Users/cmbaker9/Documents/MTOOLS'))
addpath(genpath('/Users/cmbaker9/Documents/Research/Lab_Experiments/codes/insitu'))
addpath(genpath('/Users/cmbaker9/Documents/Research/Lab_Experiments/codes/trc_lab_experiment/toolbox'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trial info
Tinfo.Hs = 0.30;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
Tinfo.spread = 40;

%%%%%%%%%%%%%%%%%% END OF USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%

%% STEP 0: Prepare file structures, folders, etc.

Tinfo = trial_files(Tinfo);

% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

%% STEP 1: Load insitu data

[press,wg] = load_insitu(Tinfo);

sensors = {'11','6'};
for i = 1:length(sensors)
    sensor = sensors(i);
    [x,y,eta] = grab_press_sensor(press,sensor);
end

sensor = {'4'};
[x,y,eta] = grab_wg_sensor(wg,sensor);

%% STEP 2: Load stereo reconstruction data

% full grid
cam         = load([Tinfo.cam.datafolder,'dem_region_x',num2str(Tino.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_res',num2str((Tinfo.cam.dx)*100),'cm.mat']);
cam.time    = Tinfo.cam.timevec;
 
%% STEP 3: Load velodyne lidar data

lidar = lidar_data_prep(Tinfo);
% load([Tinfo.datapath,'data/processed/lidar/Velodyne/',Tinfo.lidar.Tdate,'.mat']);
% [Tinfo,lidar] = TRC_lidar_info(Tinfo,griddedData);
% clear griddedData % this would be after extract z elevation




%% process all

% given a structure of the trial data, T

tgrab = 30;
[iI,iC,iL] = grab_time_allinst(T,tgrab)




