% This code will plot rectified images using code from
% D_gridGenExpampleRect.m in the CIRN-Quantitative-Coastal-Imaging-Toolbox


% Set up paths and clear workspace
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\cameras'))
addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\CIRN-Quantitative-Coastal-Imaging-Toolbox/X_CoreFunctions/'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

idxdy=0.02;
ixlim=[24 35];
iylim=[-14.5 14.5];

% Trial info
Tinfo.Hs = 0.30;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
Tinfo.spread = 40;

%% STEP 1: Create paths, files and naming
% general path and names
datapath    = 'E:/';

Tinfo = trial_files(Tinfo);
Tinfo.cam.tstart = '08-30-2018-2129UTC';
Tinfo.cam.tdate = '08-30-2018-2119UTC';
% Stereo Reconstructions
Tinfo.cam = VBAR_GENERAL_TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

imagepath = ['D:\TRC_Fall_Experiment\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1_JPEG\c2\'];

[camera.time] = cam_time(Tinfo);

%% Insitu

load([datapath,'data/processed/cameras/c2_intrinsics_extrinsics.mat']);

% extrinsics = [39.35 0.02 11.03 267.7*pi/180 32.82*pi/180 0.13*pi/180];
extrinsics(1) = extrinsics(1)+0.75;

%% Create matrix to rectify images
% see code: D_gridGenExampleRect.m

localOrigin = [0, 0]; % [ x y]
localAngle =[0]; % Degrees +CCW from Original World X
localFlagInput=1;

iz=0;

%  World Extrinsics, need to make into sell
Extrinsics{1}=extrinsics;
Intrinsics{1}=intrinsics;

% %  Local Extrinsics
% localExtrinsics{k} = localTransformExtrinsics(localOrigin,localAngle,1,xtrinsics{1});

%  Create Equidistant Input Grid
[iX iY]=meshgrid([ixlim(1):idxdy:ixlim(2)],[iylim(1):idxdy:iylim(2)]);

%  Make Elevation Input Grid
iZ=iX*0+iz;

% If entered as Local
if localFlagInput==1
    % Assign local Grid as Input Grid
    localX=iX;
    localY=iY;
    localZ=iZ;
    
    % Assign world Grid as Rotated local Grid
    [ X Y]=localTransformEquiGrid(localOrigin,localAngle,0,iX,iY);
    Z=X*.0+iz;
end

teachingMode = 0;

%% plot
imageno =  10000:1:10999;
% IMrec = nan(size(X));

for i = 1:length(imageno)
    %read image
    imagefile = [imagepath,getfield(dir([imagepath,'Movie1_Scene1_c2_',sprintf('%05d',imageno(i)),'_*.jpg']),'name')];
    IM{1} = imread(imagefile);
    
    % World Rectification
    [Ir]= imageRectifier(IM,Intrinsics,Extrinsics,X,Y,Z,teachingMode);
    Irbw = double(rgb2gray(Ir))/255;
    
    IR(:,:,i) = Irbw;
    clear IM Ir Irbw
end

subname = '';%'onewindow_11movmean_detrended';
res = idxdy;
psname = [Tinfo.savefolder,'image_stack_x',num2str(ixlim(1)),'to',num2str(round(ixlim(2))),'_res',num2str(res*100),'cm',subname,'.mat'];
eval(['save -v7.3 ',psname,' IR',' X',' Y',' res']);
