function [lidar] = lidar_data_prep(Tinfo)
% This code will read and trip the lidar data file
tlidar = Tinfo.lidar.tdate;

if Tinfo.filt == 1
    load([Tinfo.datapath,'data/processed/lidar/Velodyne/',tlidar,'/gridded/',tlidar,'_Velodyne-HDL-32-Data_gridded_filtered.mat']);
else
    load([Tinfo.datapath,'data/processed/lidar/Velodyne/',tlidar,'/gridded/',tlidar,'_Velodyne-HDL-32-Data_gridded.mat']);
end

lidar.Hz = 10;
starttemp   = datenum(tlidar,'yyyy-mm-dd-HH-MM-SS');
endtemp     = datenum(tlidar,'yyyy-mm-dd-HH-MM-SS')+(length(squeeze(griddedData.z(1,1,:)))*datenum(0,0,0,0,0,1/lidar.Hz));
timetemp  = starttemp:datenum(0,0,0,0,0,1/lidar.Hz):endtemp;

[temp,istart] = nanmin(abs(Tinfo.cam.timevec(1)-timetemp));
[temp,iend] = nanmin(abs(Tinfo.cam.timevec(end)-timetemp));
lidar.time = squeeze(timetemp(istart:iend));

% offset for time sync
istart = istart+Tinfo.offsets(2);
iend = iend+Tinfo.offsets(2);

% extract time series of interest
lidar.z       = squeeze(griddedData.z(:,:,istart:iend));
% lidar.tFrame  = squeeze(griddedData.tFrame(istart:iend));
% lidar.numPts  = squeeze(griddedData.numPts(:,:,istart:iend));
% lidar.tAll    = squeeze(griddedData.tAll(:,:,istart:iend));
lidar.x       = squeeze(griddedData.xvec);
lidar.y       = squeeze(griddedData.yvec);