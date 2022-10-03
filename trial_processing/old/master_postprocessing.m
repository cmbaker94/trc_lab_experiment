% Master post-processing run for trials
% This code can be used after dems are created in metashape. Before running
% this code, enter the trial information in trial_files.m

% Set up paths and clear workspace
clc 
clear all
close all
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\CIRN-Quantitative-Coastal-Imaging-Toolbox/X_CoreFunctions/'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))
cd E:\code\trc_lab_experiment\trial_processing

%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sprd = [0, 20, 30, 40, 0, 10, 20, 30, 40]%, 40, 40, 20, 40];%[30,  40];%
Hs = [0.3, 0.3, 0.3, 0.3, 0.25, 0.25, 0.25, 0.25, 0.25]%, 0.2, 0.3, 0.3, 0.3]; %[0.25, 0.25];%
Tp = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3];
h = [1.07, 1.07, 1.07, 1.07, 1.07, 1.07, 1.07, 1.07, 1.07, 1.07, 1.00, 1.07, 1.07];

% sprd = [0, 20, 30, 40, 40, 40, 20, 40];
% Hs = [0.3, 0.3, 0.3, 0.3, 0.2, 0.3, 0.3, 0.3]; 
% Tp =[2, 2, 2, 2, 2, 2, 3, 3];
% h = [1.07, 1.07, 1.07, 1.07, 1.07, 1.00, 1.07, 1.07];
res = [0.05];%[1 0.5 0.2 0.05 0.02];

for op = 1:length(sprd)
    for ires = 1:length(res)
% trial info
% Trial info
    Tinfo.Hs = Hs(op);
    Tinfo.Tp = Tp(op);
    Tinfo.tide = h(op);
    Tinfo.spread = sprd(op);
    Tinfo = trial_files(Tinfo);

% filter data
Tinfo.filt = 1; % 1 if true

% % insitu details
% % surf zone
% sz.day     = '7';
% sz.trial   = '10';
% 
% % inner shelf
% is.day     = '10';
% is.trial   = '03';

%%%%%%%%%%%%%%%%%%% BEGIN RAW DATA PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%

% % Process in situ raw data, lidar raw data, and dems
% 
% % in situ data
% % surf zone
% szname = read_insitu_raw_data(sz.day,sz.trial)
% % inner shelf
% isname = read_insitu_raw_data(is.day,is.trial)
% 
% % lidar data
% extract_velo_lidar(Tinfo)

% % dems
trialloc = [0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; % trialloc: if old loction = 0, if new location = 1
% extraction: point = 0, y-transect = 1, x-transect = 2, region = 3
% for extraction = 3%1:3
%     extract_tif_4D(Tinfo,trialloc(op),extraction)
%     display('change the region back to 13 to 13!!!!! currently 14 to 14')
% end
% changed the region to be dx and dy saving name. the ones with just dx =
% 0.05 cm are actually 0.25 m resolution in cross-shore. This is what the
% code is reading in currently for most of the analysis.
% 
% % read and store intrinsics and extrinsics
% extract_camera_intrisics_extrinsics(Tinfo,trialloc(op))
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%d%%%% Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if Tinfo.filt == 1
%     filter_data(Tinfo)
%     display('change the region back to 13 to 13!!!!! currently 14 to 14')
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%% Extract timeseries %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% stereo and lidar grid
% medfilt = 1; % filter camera spatially with median filter
% extract_timeseries_grid(Tinfo,medfilt)



% % stereo and lidar in situ loc
% extract_timeseries_insitu_loc(Tinfo)
% 
% %%%%%%%%%%%%%%%%%%%%%%%% Wave height, period %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% compute
% master_calc_Hs(Tinfo)
% % 
% extract Hs Tp at locations
% extract_Hs_stats(Tinfo)
% 
% % plot
% plot_Hs_spectra(Tinfo)
% 
% %%%%%%%%%%%%%%%%%%%%%% Directional Spread %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % puv
% calc_insitu_dirspectra_PUV(Tinfo)
% 
% % possibly python?
% hyd = 1;
% calc_directional_spec_iter_press(Tinfo,hyd)
% calc_directional_spec_iter_cam(Tinfo)
% calc_directional_spec_iter_lidar(Tinfo)

% 
% remote sensing and pressure gage instruments and interations
rsinst = 14;
rsiter = 120;
pginst = 5;
pgiter = 120;
frange = [0.3 0.8]; % range over which to compute the S(theta) for the cos2s dist fitting
matorpy = 1; % 0 for python, 1 for mat
% estimate_spread(Tinfo,rsinst,rsiter,pginst,pgiter,frange,matorpy)
% display('redo this, especially for 0 deg spread, with less oints withere max is selected')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% Crest length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % imagery
% calc_wave_crest_length_image(Tinfo)

% % sea-surface elevation
% % THIS CODE ISN'T USING FILTERED DATA
% calc_wave_crest_length_stereo(Tinfo)
% 
% THRESHOLD STEREO
gamma = 0.27;%0.7%0.6;%0.6 % image threshold
minarea = 150;%800 % min area of object
samprate = 8; %Hz
pltflag = 0;
ixlim=[26.5 31.5];%25
iylim = [-14 14];
idxdy = 0.05;
% threshold_stereo_crest_length(Tinfo,gamma,minarea,samprate,pltflag,ixlim,iylim,idxdy)
% 
% rectify images larger version
% idxdy=res(ires);%0.02;
% ixlim=[25 35];%[18 35];
% iylim=[-14 14];
% imlen=0;
% rectify_images_met(Tinfo,idxdy,ixlim,iylim,imlen)
% 
% RECTIFY IMAGES FOR CREST LENGTH INDENTIFICATION
idxdy=0.05;
ixlim=[25 33.5];%25
iylim=[-14 14];
imlen=0;
%%%% rectify_images_met(Tinfo,idxdy,ixlim,iylim,imlen)
rectify_images_stereo(Tinfo,idxdy,ixlim,iylim,imlen)

% 
% % % threshold images
thresh = 0.35;%0.7%0.6;%0.6 % image threshold
minarea = 250;%30;%800 % min area of object
samprate =8; %Hz
pltflag = 1;
ixsel = [26.5 31.5];
% threshold_image_crest_length(Tinfo,thresh,minarea,samprate,pltflag,ixlim,iylim,idxdy,ixsel)
%
% THRESHOLD COMBINED
idxdy=0.05;
ixlim=[25 33.5];%25
iylim=[-14 14];
imlen=0;
gammab = 0.2;
% IMthresh = 0.38; % image threshold
IMthresh =  0.2; % adjust the mean cross-shore pixel intensity for thresholding imagery
minareaIMST = 40;% min area of crest for stereo and imagery
minarea = 100; % min area of crest combined
distconnect = 6;
samprate =8; %Hz
ixsel = [27 31.5];
% threshold_combined_crest_length(Tinfo,gammab,IMthresh,minareaIMST,minarea,distconnect,samprate,pltflag,ixlim,iylim,idxdy,ixsel)

% %% Copy files
% Tinfo = trial_files(Tinfo);
% Tinfo.cam = TRC_camera_info(Tinfo.cam);
% Tinfo = wc_comp_store(Tinfo);
% clustfn = '/data2/baker/lab_experiments/data/processed/conditions/';
% clear Tinfo

%% surface currents
% rectify images larger version
idxdy=0.01;%0.02;
ixlim=[18 35];%[18 35];
iylim=[-14 14];
imlen=1; %1 if rectify all trial
rectify_images_stereo(Tinfo,idxdy,ixlim,iylim,imlen)

a = 4;%[2 4];%4 6 6 8 8 10];
o = 2;% [1 2];% 2 2 3 2 4 5];
for iao = 1:length(a)
idxdy=res(ires);%0.02;
ixlim=[25 35];%[18 35];
iylim=[-14 14];
awin = a(iao); % Averaging Window (seconds)
ovr = o(iao); % Overlap (seconds)
% rectifiedWAM(Tinfo,ixlim,iylim,idxdy,awin,ovr)
end
    end
end
%didn't finish for  deg deg
addpath('E:\code\trc_lab_experiment\wave_statistics\crest\')
Hs = calc_comb_thresh_crest();

%% vbar

% % rectify images
% idxdy=0.02;
% ixlim=[24 36];
% iylim=[-14.5 14.5];
% rectify_images(Tinfo,idxdy,ixlim,iylim)