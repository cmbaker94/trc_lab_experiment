function rectify_images(Tinfo,idxdy,ixlim,iylim,imlen)

% This code will plot rectified images using code from
% D_gridGenExpampleRect.m in the CIRN-Quantitative-Coastal-Imaging-Toolbox


% Set up paths and clear workspace
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\cameras'))
addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\CIRN-Quantitative-Coastal-Imaging-Toolbox/X_CoreFunctions/'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% idxdy=0.02;
% ixlim=[25 31.5];
% iylim=[-14.5 14.5];

% % Trial info
% Tinfo.Hs = 0.30;
% Tinfo.Tp = 2;
% Tinfo.tide = 1.07;
% Tinfo.spread = 0;

%% STEP 1: Create paths, files and naming
% general path and names
datapath    = 'E:/';

Tinfo = trial_files(Tinfo);

% Stereo Reconstructions
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

if Tinfo.spread == 0
    if Tinfo.Hs == 0.3
        imagepath = ['D:\TRC_Fall_Experiment\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1_JPEG\c2\'];
    else 
        imagepath = ['H:\TRC_Fall_Experiment\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1_JPEG\'];
    end
elseif Tinfo.spread == 40
    if Tinfo.Hs == 0.25
        imagepath = ['D:\TRC_Fall_Experiment\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1_JPEG\'];
%     elseif Tinfo.Hs == 0.3
%         imagepath = ['D:\TRC_Fall_Experiment\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1_JPEG\c2\'];
%         display('revert to H!!!')
    else
        imagepath = ['H:\TRC_Fall_Experiment\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1_JPEG\'];
    end
else
    imagepath = ['H:\TRC_Fall_Experiment\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1_JPEG\'];
end

[camera.time] = cam_time(Tinfo);

%% Insitu

cam = load([datapath,'data/processed/cameras/c2_intrinsics_extrinsics.mat']);

% extrinsics = [39.35 0.02 11.03 267.7*pi/180 32.82*pi/180 0.13*pi/180];
% extrinsics(1) = extrinsics(1);

%% Create matrix to rectify images
% see code: D_gridGenExampleRect.m

localOrigin = [0, 0]; % [ x y]
localAngle =[0]; % Degrees +CCW from Original World X
localFlagInput=1;

iz=Tinfo.tide;

%  World Extrinsics, need to make into sell
Ex{1}=cam.extrinsics;
In{1}=cam.intrinsics;

% %  Local Extrinsics
% localExtrinsics{k} = localTransformExtrinsics(localOrigin,localAngle,1,xtrinsics{1});

%  Create Equidistant Input Grid
[X Y]=meshgrid([ixlim(1):idxdy:ixlim(2)],[iylim(1):idxdy:iylim(2)]);

%  Make Elevation Input Grid
Z=X*0+iz;

% % If entered as Local
% if localFlagInput==1
%     % Assign local Grid as Input Grid
%     localX=iX;
%     localY=iY;
%     localZ=iZ;
%     
%     % Assign world Grid as Rotated local Grid
%     [ X Y]=localTransformEquiGrid(localOrigin,localAngle,0,iX,iY);
%     Z=X*.0+iz;
% end

teachingMode = 0;

%% plot
imageno =  str2double(Tinfo.cam.imagerange(1:5)):1:str2double(Tinfo.cam.imagerange(7:end));
% IMrec = nan(size(X));
eval(['!mkdir ',Tinfo.savefolder(1:74),'orthos'])
subname = '';%'onewindow_11movmean_detrended';
res = idxdy;
odir = [Tinfo.savefolder(1:74),'orthos\','x',num2str(ixlim(1)),'to',num2str(round(ixlim(2))),'_y',num2str(round(iylim(2))),'_res',num2str(res*100),'cm',subname];
eval(['!mkdir ',odir])

if imlen == 0
    for i = 1:length(imageno)%266:8*60%length(imageno)
        %read image
        imagefile = [imagepath,getfield(dir([imagepath,'Movie1_Scene1_c2_',sprintf('%05d',imageno(i)),'_*.jpg']),'name')];
        IM{1} = imread(imagefile);
        
        % World Rectification
        [Ir]= imageRectifier(IM,In,Ex,X,Y,Z,teachingMode);
        %     Irbw = double(rgb2gray(Ir))/255;
        
        %     IR(:,:,i) = Irbw;
        
        % Save Image
        imwrite(flipud(Ir),string(fullfile(odir,['c2_',sprintf('%05d',imageno(i)),'.tiff'])))
        display(imageno(i))
        clear IM Ir Irbw
    end
elseif imlen == 1
    imageno = dir([imagepath,'Movie1_Scene1_c2_*']);
    for i = 1:length(imageno)%266:8*60%length(imageno)
        %read image
        imagefile = [imagepath,getfield(dir([imagepath,'Movie1_Scene1_c2_',sprintf('%05d',i-1),'_*.jpg']),'name')];
        IM{1} = imread(imagefile);
        
        % World Rectification
        [Ir]= imageRectifier(IM,In,Ex,X,Y,Z,teachingMode);
        %     Irbw = double(rgb2gray(Ir))/255;
        
        %     IR(:,:,i) = Irbw;
        
        % Save Image
        imwrite(flipud(Ir),string(fullfile(odir,['c2_',sprintf('%05d',i-1),'.tiff'])))
        display(i-1)
        clear IM Ir Irbw
    end
end

% subname = '';%'onewindow_11movmean_detrended';
% res = idxdy;
% psname = [Tinfo.savefolder,'image_stack_x',num2str(ixlim(1)),'to',num2str(round(ixlim(2))),'_res',num2str(res*100),'cm',subname,'.mat'];
% eval(['save -v7.3 ',psname,' IR',' X',' Y',' res']);
