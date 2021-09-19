function [lidar] = calc_Hs_lidar(Tinfo,sdet)
% This function will compute the sea-surface elevation spectra and compute
% the significant wave height from lidar measurments over a specified time
% range. 
% INPUT:
% lidarfile     = name of lidar output file
% trange        = time wave statistics will be computed over, 2x1 vector:
%                 [start, end]
% sdet          = structure with defined parameters for computing the
%                 sea-surface elevation spectra and Hs
%       WL      = window length, default = []
%       OL      = overlap length, default = []
%       Hz      = measurement frequency, default = 10 Hz
%       frange  = frequency range Hs integrated over, 2x1 vector: 
%                 [start, end], default = [0.25, 1] Hz
%       nancutoff = percent of nan's allowable in timeseries when computing
%                   the spectra, default = 10%
%       plotspec = flag, plot each spectra as computed, where 1 = plot
%                  and 0 = no plots, default = 0 
% OUTPUT:
% lidar         = structure with spectra (m^2/Hz), Hs (m), frequency (Hz)

if Tinfo.filt == 1
    load(['E:\data\processed\lidar\Velodyne\',Tinfo.lidar.tdate,'\gridded\',Tinfo.lidar.tdate,'_Velodyne-HDL-32-Data_gridded_filtered.mat']);
else
    load(['E:\data\processed\lidar\Velodyne\',Tinfo.lidar.tdate,'\gridded\',Tinfo.lidar.tdate,'_Velodyne-HDL-32-Data_gridded.mat']);
end
% set defaults if not specified in input
if ~isfield(sdet,'WL')
    sdet.WL = [];
    sdet.OL = [];
end
if ~isfield(sdet,'Hz')
    sdet.Hz = 10;
end
if ~isfield(sdet,'frange')
    sdet.frange = [0.25, 1.25];
end
if ~isfield(sdet,'nancutoff')
    sdet.nancutoff = 0.1;
end
if ~isfield(sdet,'plotspec')
    sdet.plotspec = 0;
end


% lidar record start and end time
starttemp           = datenum(Tinfo.lidar.tdate,'yyyy-mm-dd-HH-MM-SS');
endtemp             = datenum(Tinfo.lidar.tdate,'yyyy-mm-dd-HH-MM-SS')+(length(squeeze(griddedData.z(1,1,:)))*datenum(0,0,0,0,0,1/sdet.Hz));
griddedData.time    = starttemp:datenum(0,0,0,0,0,1/sdet.Hz):endtemp;

% index of analysis start and end 
[temp,istart]       = nanmin(abs(Tinfo.cam.timevec(1)-griddedData.time));
[temp,iend]         = nanmin(abs(Tinfo.cam.timevec(end)-griddedData.time));

% extract time series of interest
z       = squeeze(griddedData.z(:,:,istart:iend));
% tFrame  = squeeze(griddedData.tFrame(istart:iend));
% numPts  = squeeze(griddedData.numPts(:,:,istart:iend));
% tAll    = squeeze(griddedData.tAll(:,:,istart:iend));
x       = squeeze(griddedData.xvec);
y       = squeeze(griddedData.yvec);
% time    = squeeze(griddedData.time(istart:iend));
clear griddedData *temp

% LiDAR Bulk Stats
% stillwat = nanmean(nanmean(nanmean(squeeze(z(6:50,12:20,1:30),3))));
% lidar.Hs_4std = 4*nanstd(lidar.z,[],3);
% lidar.zmean = nanmean(lidar.z,3);

% LiDAR spectral stats
[lidar.Hs,lidar.Tp,lidar.See,lidar.f,lidar.Seec,lidar.nanratio] = calc_Hs_spectra_remotesensing(z,x,y,sdet.WL,sdet.OL,sdet.frange,sdet.Hz,sdet.nancutoff,sdet.plotspec);
lidar.x = x;
lidar.y = y;
