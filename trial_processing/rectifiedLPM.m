% function rectifiedLPM(Tinfo,ixlim,iylim,idxdy,fcut)
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

Tinfo.Hs = 0.3;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
Tinfo.spread = 30;
idxdy   = 0.02;
ixlim   = [18 35];
iylim   = [-14 14];
fps     = 8;
cutoff  = 4; %s

subname = 'trimmed';
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

display('TEMP TRIMING')
% L = L(1:(60/cutoff)*fps*3,:);
% L = L(7201:12000,:);

t = str2num(L)/Tinfo.cam.Hz;
dt=mode(diff(t));
ts=(t-t(1));

%% Create Structure to Hold Images
Io=imread(fullfile(odir,['c2_',L(1,:),'.tiff']));
[r c co]=size(Io);
clear('Io')

%% Load Images

for k=1:length(L)
    
% Load Image (Gray- Matlab cannot handle RGB array size)
I(:,:,k)=double(rgb2gray(imread(fullfile(odir,['c2_',L(k,:),'.tiff']))));

k
end

disp('Loaded Images')

%% Selected time

Tbin= (cutoff/2):cutoff/2: (ts(end)-cutoff/2); % Centered On
Tbin2 = (cutoff/2):cutoff/4: (ts(end)-cutoff/2); % Centered On
% Tlow=Tbin-cutoff/2; % Upper Limit
% Thigh=Tbin+cutoff/2; % Lower Limit

%% Filter video

for i = 1:r
    for j = 1:c
        Ilp(i,j,:) = pl64ta(squeeze(I(i,j,:)),cutoff*fps);
    end
    i
end


%% Take Average

Ilpsub = uint8(squeeze(Ilp(:,:,Tbin*fps)));
Ilpsub2 = uint8(squeeze(Ilp(:,:,Tbin2*fps)));
Ilp = uint8(Ilp);

TbinS=Tbin;
Tbin=Tbin+t(1);
disp('Averaged Images')

%% Save File (MAT)

oname=['c2_',L(1,:),'_',L(end,:),'_cutoff',num2str(cutoff),subname];
sname = fullfile(Tinfo.savefolder(1:74),'WAMs',geoname,oname);
eval(['!mkdir ',sname])
save(fullfile(sname,[oname 'R.mat']),'Ilp','Tbin','TbinS','cutoff')
disp('saved Mat File')
save(fullfile(sname,['c2_',L(1,:),'_',L(end,:),'_nofiltR.mat']),'I','t')
disp('saved Mat File')

%% Make A Movie
    
v = VideoWriter(fullfile(sname,[oname,'R.avi']));
v.FrameRate=24; %?????
open(v)

for k=1:length(t)
writeVideo(v,Ilp(:,:,k))

k
end  
   close(v) 
   
%% plot instantenous 
inname = ['c2_',L(1,:),'_',L(end,:),'_instaneous'];
   
v = VideoWriter(fullfile(sname,[inname,'R.avi']));
v.FrameRate=24; %?????
open(v)

for k=1:length(t)
writeVideo(v,uint8(I(:,:,k)))

k
end  
   close(v) 
   
   %% Make A Movie

eval(['!mkdir ',sname,'\freq',num2str(cutoff/2)])
v = VideoWriter(fullfile(sname,['freq',num2str(cutoff/2)],[oname,'_subsamp',num2str(cutoff/2),'_R.avi']));
v.FrameRate=1; %?????
open(v)

for k=1:length(Tbin)
writeVideo(v,Ilpsub(:,:,k))

k
end  
   close(v) 
   disp('saved movie')
   
   
   %% Make A Movie

eval(['!mkdir ',sname,'\freq',num2str(cutoff/4)])
v = VideoWriter(fullfile(sname,['freq',num2str(cutoff/4)],[oname,'_subsamp',num2str(cutoff/4),'_R.avi']));
v.FrameRate=1; %?????
open(v)

for k=1:length(Tbin2)
writeVideo(v,Ilpsub2(:,:,k))

k
end  
   close(v) 
   disp('saved movie')

