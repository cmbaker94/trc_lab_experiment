function [Hs,Tp] = grab_Hs_Tp(spec,gage)
% function to extract Hs and tp from mat files.

eval(['Hs(1,1) = spec.press',gage,'.press.Hs;'])
eval(['Tp(1,1) = spec.press',gage,'.press.Tp;'])

eval(['Hs(1,2) = spec.press',gage,'.cam.Hs;'])
eval(['Tp(1,2) = spec.press',gage,'.cam.Tp;'])

eval(['Hs(1,3) = spec.press',gage,'.lidar.Hs;'])
eval(['Tp(1,3) = spec.press',gage,'.lidar.Tp;'])