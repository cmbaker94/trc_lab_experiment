% function rectifiedWAM(Tinfo,ixlim,iylim,idxdy,awin,ovr)
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
    Tinfo.spread = 30;


% filter data
Tinfo.filt = 1; % 1 if true

idxdy=0.05;
ixlim=[25 35];
iylim=[-14 14];
imlen=0;
imsel = [0:359]+7200;% [0:8*15*60-1];%[7200:7200+4800-1]';
fps = 1;%8;

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
% L = cnames(:,4:8);
for i = 1:length(imsel)
    L(i,:) = num2str(imsel(i),'%05.f');
end
t = str2num(L)/Tinfo.cam.Hz;
dt=mode(diff(t));
ts=(t-t(1));

%% Create Structure to Hold Images
Io=imread(fullfile(odir,['c2_',L(1,:),'.tiff']));
[r c co]=size(Io);
clear('Io')

oname=['_c2_',num2str(imsel(1)),'_',num2str(imsel(end))];
sname = fullfile(Tinfo.savefolder(1:74),'WAMs',geoname,oname);
eval(['!mkdir ',sname])

%% Load Images
v = VideoWriter([odir,oname,'R_framerate',num2str(fps),'.avi']);
v.FrameRate=fps; %?????
open(v)

for k=1:length(L)
    
% Load Image (Gray- Matlab cannot handle RGB array size)
% I=double(rgb2gray(imread(fullfile(odir,['c2_',L(k,:),'.tiff']))));
% I=imread(fullfile(odir,['c2_',L(k,:),'.tiff']));
I=rgb2gray(imread(fullfile(odir,['c2_',L(k,:),'.tiff'])));

writeVideo(v,I)

k
end  
close(v)
disp('saved movie')

