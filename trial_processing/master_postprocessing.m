% Master post-processing run for trials
% This code can be used after dems are created in metashape. Before running
% this code, enter the trial information in trial_files.m

clc 
clear all
close all
cd E:\code\trc_lab_experiment\trial_processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sprd = [0, 20, 30, 40, 0, 10, 20, 30, 40, 40, 40, 20, 40];
Hs = [0.3, 0.3, 0.3, 0.3, 0.25, 0.25, 0.25, 0.25, 0.25, 0.2, 0.3, 0.3, 0.3];
Tp = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3];
h = [1.07, 1.07, 1.07, 1.07, 1.07, 1.07, 1.07, 1.07, 1.07, 1.07, 1.00, 1.07, 1.07];

for op = 1:length(sprd)
% trial info
% Trial info
    Tinfo.Hs = Hs(op);
    Tinfo.Tp = Tp(op);
    Tinfo.tide = h(op);
    Tinfo.spread = sprd(op);

% filter data
Tinfo.filt = 0; % 1 if true

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

% % % dems
% trialloc = 1; % trialloc: if old loction = 0, if new location = 1
% % extraction: point = 0, y-transect = 1, x-transect = 2, region = 3
% for extraction = 1:3
%     extract_tif_4D(Tinfo,trialloc,extraction)
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%d%%%% Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if Tinfo.filt == 1
%     filter_data(Tinfo)
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%% Extract timeseries %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % stereo and lidar grid
% extract_timeseries_grid(Tinfo)
% 
% % stereo and lidar in situ loc
% extract_timeseries_insitu_loc(Tinfo)
% 
% %%%%%%%%%%%%%%%%%%%%%%%% Wave height, period %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % compute
% master_calc_Hs(Tinfo)
% % 
% % extract Hs Tp at locations
% extract_Hs_stats(Tinfo)
% 
% plot
plot_Hs_spectra(Tinfo)
% 
% %%%%%%%%%%%%%%%%%%%%%% Directional Spread %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % puv
% calc_insitu_dirspectra_PUV(Tinfo)
% 
% % possibly python?
% 
% % remote sensing and pressure gage instruments and interations
% rsinst = 14;
% rsiter = 200;
% pginst = 5;
% pgiter = 200;
% frange = [0.3 0.8]; % range over which to compute the S(theta) for the cos2s dist fitting
% estimate_spread(Tinfo,rsinst,rsiter,pginst,pgiter,frange)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% Crest length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % imagery
% calc_wave_crest_length_image(Tinfo)
% 
% % sea-surface elevation
% % THIS CODE ISN'T USING FILTERED DATA
% calc_wave_crest_length_stereo(Tinfo)
% 
% % threshold stereo
% gamma = 0.4;%0.7%0.6;%0.6 % image threshold
% minarea = 20;%800 % min area of object
% samprate = 2; %Hz
% pltflag = 0;
% ixlim=[25 31.5];
% threshold_stereo_crest_length(Tinfo,gamma,minarea,samprate,pltflag,ixlim)

% % rectify images
% idxdy=0.02;
% ixlim=[25 31.5];
% iylim=[-14.5 14.5];
% rectify_images(Tinfo,idxdy,ixlim,iylim)
% 
% % threshold images
% thresh = 0.45;%0.7%0.6;%0.6 % image threshold
% minarea = 20;%30;%800 % min area of object
% samprate = 2; %Hz
% pltflag = 0;
% threshold_image_crest_length(Tinfo,thresh,minarea,samprate,pltflag)
% 
% %% Copy files
% Tinfo = trial_files(Tinfo);
% Tinfo.cam = TRC_camera_info(Tinfo.cam);
% Tinfo = wc_comp_store(Tinfo);
% clustfn = '/data2/baker/lab_experiments/data/processed/conditions/';
% clear Tinfo
end
%didn't finish for  deg deg

%% vbar

% % rectify images
% idxdy=0.02;
% ixlim=[24 36];
% iylim=[-14.5 14.5];
% rectify_images(Tinfo,idxdy,ixlim,iylim)