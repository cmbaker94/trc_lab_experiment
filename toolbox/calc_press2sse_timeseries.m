% this is just a testing code for now
clear all
close all
clc
addpath(genpath('/Users/cmbaker9/Documents/MTOOLS'))


data = '/Users/cmbaker9/Documents/Research/Lab_Experiments/data/processed/DEM/TRM-09-01-2018-2155UTC_Scene1/frames_07200-11999/';
file = 'pressure_gage_timeseries.csv';
A = csvread([data,file])
press = A(:,1);
Hz = 100;
maxfac = 1.2;
offset = 0.05;

[SSE] = press2sse_timeseries(press,Hz,offset,maxfac);