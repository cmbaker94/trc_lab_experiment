function [Ir,X,Y] = rectify_image_single(Tinfo,imagefile,idxdy,ixlim,iylim,imageno)

% This code will plot rectified images using code from
% D_gridGenExpampleRect.m in the CIRN-Quantitative-Coastal-Imaging-Toolbox

% Set up paths and clear workspace
% addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
% addpath(genpath('E:\code\cameras'))
% addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\CIRN-Quantitative-Coastal-Imaging-Toolbox/X_CoreFunctions/'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

% camprop = load('E:/data/processed/cameras/c2_intrinsics_extrinsics.mat');
cam = load([Tinfo.cam.datafolder,'c2_prop_guiBathy.mat']);
cam.extrinsics = [cam.X cam.Y cam.Z cam.az cam.tilt cam.roll];
cam.intrinsics = [cam.NU cam.NV cam.cUo cam.cVo cam.Fx cam.Fy cam.d1 cam.d2 cam.d3 cam.t1 cam.t2];

%% load  stere

if Tinfo.filt == 1
    F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(iylim(1))),'to',num2str(iylim(end)),'m_resx',num2str((Tinfo.cam.dx)*100),'cm_resy',num2str((Tinfo.cam.dx)*100),'cm_filtered.mat'];
elseif Tinfo.filt == 0
    F1 = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(iylim(1))),'to',num2str(iylim(end)),'m_resx',num2str((Tinfo.cam.dx)*100),'cm_resy',num2str((Tinfo.cam.dx)*100),'cm.mat'];
end
stereo =  load(F1,'x','y');
F1 = matfile(F1);
stereo.z = squeeze(F1.z(:,:,imageno-Tinfo.cam.imagestart));
%% Create matrix to rectify images
% see code: D_gridGenExampleRect.m

localOrigin = [0, 0]; % [ x y]
localAngle =[0]; % Degrees +CCW from Original World X
localFlagInput=1;

iz=Tinfo.tide;

%  World Extrinsics, need to make into sell
Extrinsics{1}=cam.extrinsics;
Intrinsics{1}=cam.intrinsics;

%  Create Equidistant Input Grid
[X Y]=meshgrid([ixlim(1):idxdy:ixlim(2)],[iylim(1):idxdy:iylim(2)]);

%  Make Elevation Input Grid
% Z=X*0+iz;
Z  = prepZrect(Tinfo,stereo.x,stereo.y,stereo.z,X,Y);

% znotnan = ~isnan(Z);
% Z = interp2(X(~isnan(Z)),Y(~isnan(Z)),Z(~isnan(Z)),X,Y)%,'spline');
% zxavg(10:end-10) = movmedian(zxavg(10:end-10),20,'omitnan');

teachingMode = 0;

%read image
IM{1} = imread(imagefile);

% World Rectification
[Ir]= imageRectifier(IM,Intrinsics,Extrinsics,X,Y,Z,teachingMode);
Ir = flipud(Ir);

