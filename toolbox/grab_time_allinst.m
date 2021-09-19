function [iI,iC,iL] = grab_time_allinst(T,tgrab)
% grab index of the of a specificed time from the camera, based on seconds
% after sterio results started

iC = round(tgrab*8);
tcam = T.camera.time(iC);
[temp,iI] = nanmin(abs(tcam-T.press.time));
[temp,iL] = nanmin(abs(tcam-T.lidar.time));