function [Ir,X,Y] = rectify_image_single(Tinfo,imagefile,idxdy,ixlim,iylim)

% This code will plot rectified images using code from
% D_gridGenExpampleRect.m in the CIRN-Quantitative-Coastal-Imaging-Toolbox

% Set up paths and clear workspace
% addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
% addpath(genpath('E:\code\cameras'))
% addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\CIRN-Quantitative-Coastal-Imaging-Toolbox/X_CoreFunctions/'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

camprop = load('E:/data/processed/cameras/c2_intrinsics_extrinsics.mat');

%% Create matrix to rectify images
% see code: D_gridGenExampleRect.m

localOrigin = [0, 0]; % [ x y]
localAngle =[0]; % Degrees +CCW from Original World X
localFlagInput=1;

iz=Tinfo.tide;

%  World Extrinsics, need to make into sell
Extrinsics{1}=camprop.extrinsics;
Intrinsics{1}=camprop.intrinsics;

%  Create Equidistant Input Grid
[X Y]=meshgrid([ixlim(1):idxdy:ixlim(2)],[iylim(1):idxdy:iylim(2)]);

%  Make Elevation Input Grid
Z=X*0+iz;

teachingMode = 0;

%read image
IM{1} = imread(imagefile);

% World Rectification
[Ir]= imageRectifier(IM,Intrinsics,Extrinsics,X,Y,Z,teachingMode);
Ir = flipud(Ir);

