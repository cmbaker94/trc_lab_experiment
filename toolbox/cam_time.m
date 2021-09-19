function [time] = cam_time(Tinfo)
% Creates a time vector for the camera system
% INPUT: Tcam = trial information structure
% OUTPUT: time = vector of time for cameras in matlab format

starttemp       = datenum(Tinfo.cam.tstart(1:end-3),'mm-dd-yyyy-HHMM')+datenum(0,0,0,0,0,Tinfo.cam.imagestart/Tinfo.cam.Hz);
endtemp         = datenum(Tinfo.cam.tstart(1:end-3),'mm-dd-yyyy-HHMM')+datenum(0,0,0,0,0,(Tinfo.cam.imagestart+Tinfo.cam.numframes-1)/Tinfo.cam.Hz);
time        = starttemp:datenum(0,0,0,0,0,1/Tinfo.cam.Hz):endtemp;