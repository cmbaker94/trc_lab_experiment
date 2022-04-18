% function rectifiedvideo(Tinfo,ixlim,iylim,idxdy,awin,ovr)
% function to create a wave-averaged movie. credits to B. Bruder for code
% INPUT:
% Tinfo             trial information
% ixlim, iylim      x,y limitis (m)
% idxdy             resolution (m)
% awin              averaging window (s)
% ovr               overlap window (s)
% OUTPUT:
% saved mat and video vile of wave averaged movie
clc 
clear all
close all
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\CIRN-Quantitative-Coastal-Imaging-Toolbox/X_CoreFunctions/'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))
cd E:\code\trc_lab_experiment\trial_processing

Tinfo.Hs = 0.25;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
Tinfo.spread = 40;
idxdy=0.02;
ixlim=[18 35];
iylim=[-14 14];

%% STEP 1: Create paths, files and naming

Tinfo = trial_files(Tinfo);

% Stereo Reconstructions
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

%% Load File Names + Times
subname = '';
geoname = ['x',num2str(ixlim(1)),'to',num2str(round(ixlim(2))),'_y',num2str(round(iylim(2))),'_res',num2str(idxdy*100),'cm',subname];
odir = [Tinfo.savefolder(1:74),'orthos\',geoname];

cnames = ls([odir,'\c2_*.tiff']);
L = cnames(:,4:8);
t = str2num(L)/Tinfo.cam.Hz;
dt=mode(diff(t));
ts=(t-t(1));

%% Create Structure to Hold Images
Io=imread(fullfile(odir,['c2_',L(1,:),'.tiff']));
[r c co]=size(Io);
% Iwam=zeros(r,c,length(Tbin));
clear('Io')


%% Make A Movie

v = VideoWriter([Tinfo.savefolder(1:74),'dirsprd_rectified_',geoname,'.avi']);
v.FrameRate=8;
open(v)

for k=1:60*8
    
    I=imread(fullfile(odir,['c2_',L(k,:),'.tiff']));
    writeVideo(v,I)

k
end  
   close(v) 
   disp('saved movie')

