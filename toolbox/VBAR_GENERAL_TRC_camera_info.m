function cam = VBAR_GENERAL_TRC_camera_info(cam)
% This code is used to import the default information about the camera
% setup as well as specify folder paths.
% INPUT:
% Tinfo.tstart:   date\time trial was actually started
% Tinfo.tdate:    date\time of trial file
% OUTPUT:
% Tinfo:     structure with camera information
%
% CB created on Feb 7, 2020

% general path and names
datapath    = 'E:\';

% Stereo Reconstructions Information
cam.camerasys   = 'TRM';                     % camera setup - e.g., TRM (offshore) or TRC (onshore)
cam.scene       = '1';                       % scene number of trial - typ 1
cam.imagestart  = 10000;                      % images number of first frame on file 
cam.numframes   = 1000;                      % number of frames processed
cam.dx          = 0.05;
cam.dy          = 0.05;
cam.regx        = [25:cam.dx:37];
cam.regy        = [-13:cam.dy:13];
cam.Hz           = 8;

% create file paths
% cam.trialname   = [cam.camerasys,'-',cam.tdate];
cam.trialname   = [cam.camerasys,'-',cam.tdate];
cam.imagerange  = [num2str(cam.imagestart,'%05.f'),'-',num2str(cam.imagestart+(cam.numframes-1),'%05.f')];
cam.trimname    = ['frames_',cam.imagerange,'\'];
cam.datafolder  = [datapath,'data\processed\cameras\',cam.trialname,'\',cam.trimname];

starttemp       = datenum(cam.tstart(1:end-3),'mm-dd-yyyy-HHMM')+datenum(0,0,0,0,0,cam.imagestart/cam.Hz);
endtemp         = datenum(cam.tstart(1:end-3),'mm-dd-yyyy-HHMM')+datenum(0,0,0,0,0,(cam.imagestart+cam.numframes-1)/cam.Hz);
cam.timevec     = starttemp:datenum(0,0,0,0,0,1/cam.Hz):endtemp;
