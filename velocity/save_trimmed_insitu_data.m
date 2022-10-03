% This code will extract and save time series from the camera dems and
% lidar mat file to compute the directional spectra as an 'array'.

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\cameras\dem'))
addpath(genpath('E:\code\insitu'))
%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trial info
Tinfo.Hs = 0.3;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
Tinfo.spread = 30;

Tinfo = trial_files(Tinfo);

% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

[camera.time] = cam_time(Tinfo);

saveinsitu = [Tinfo.savefolder,'insitu_trimmed'];
eval(['!mkdir ',saveinsitu])

%% Insitu

% load insitu data
datapath = 'E:\';
F1          = matfile([datapath,'data\processed\insitu\',Tinfo.sz.tdate,'\',Tinfo.sz.tdate,'-insitu.mat']);
pg          = F1.press;
waveg       = F1.wg;
velocities  = F1.vel;
Tinfo_insitu       = F1.Tinfo;
Hz          = pg.Hz;
maxfac = 1.2;
offset  = 0.05;

starttemp   = datenum(Tinfo_insitu.time(1:end-3));
endtemp     = starttemp+((length(waveg.wg1)-1)*datenum(0,0,0,0,0,1/Hz));
timetemp    = starttemp:datenum(0,0,0,0,0,1/Hz):endtemp;

% find overlapping range for insitu
[temp,istart] = nanmin(abs(camera.time(1)-timetemp));
[temp,iend] = nanmin(abs(camera.time(end)-timetemp));
iend = iend+1;

press.time = timetemp(istart:iend)';
press.name  = fieldnames(pg.xyzd)';

% offset for time sync
istart = istart+Tinfo.offsets(1);
iend = iend+Tinfo.offsets(1);

filtramp = 1000;
% compute bulk statistics
for i = 1:length(press.name)
    % load pressure data
    eval(['press.press(:,i)   = pg.',press.name{i},'(istart:iend);'])
    press.eta(:,i) = press2sse_timeseries(press.press(:,i),Hz,offset,maxfac,filtramp);
    eval(['press.xyz(:,i)     = pg.xyzd.',press.name{i},'(1,1:3);']) % z-elevation of sensor in tank coordinates
    press.loc{i} = ['p',sprintf('%02d',str2double(press.name{i}(6:end)))];
end

clear pg *temp T1 T2

% Compute bulk statistics from wave gages
% get wave gage names
wg.names	= fieldnames(waveg.xyzd)';
wg.time     = press.time;
% compute bulk statistics
for i = 1:length(wg.names)
    % load wave gage data
    eval(['wg.z(:,i)                     = waveg.',wg.names{i},'(istart:iend);'])
    eval(['wg.xyz(:,i)              = waveg.xyzd.',wg.names{i},'(1,1:3);']) 
    wg.loc{i} = ['wg',sprintf('%02d',str2double(wg.names{i}(3:end)))];
end

%%
% Compute bulk statistics from wave gages
% get wave gage names
vel.names	= fieldnames(velocities.xyzd)';
vel.time     = press.time;
% compute bulk statistics
for i = 1:length(vel.names)
    % load wave gage data
%     eval(['vel.u(:,i)                     = velocities.',vel.names{i}(4:end),'(istart:iend);'])
    eval(['adv',num2str(i,'%02.f'),'(:,2)    = velocities.u',vel.names{i}(4:end),'(istart:iend);'])
    eval(['adv',num2str(i,'%02.f'),'(:,3)    = velocities.v',vel.names{i}(4:end),'(istart:iend);'])
    eval(['adv',num2str(i,'%02.f'),'(:,4)    = velocities.w',vel.names{i}(4:end),'(istart:iend);'])
    eval(['T = array2table(adv',num2str(i,'%02.f'),');'])
    T.Properties.VariableNames(1:3) = {'u (m/s)','v (m/s)','w (m/s)'};
    writetable(T,[saveinsitu,'\adv',num2str(i,'%02.f'),'.csv'])
    
    eval(['vel.xyz(:,i)              = waveg.xyzd.',wg.names{i},'(1,1:3);'])
    vel.xyz(3,i)              = 1.07-press.xyz(3,i);
%     wg.loc{i} = ['wg',sprintf('%02d',str2double(wg.names{i}(3:end)))];
end
clear T
T = array2table(vel.xyz);
T.Properties.VariableNames(1:12) = {'adv1','adv2','adv3','adv4','adv5','adv6','adv7','adv8','adv9','adv10','adv11','adv12'};
writetable(T,[saveinsitu,'\adv_xyh.csv'])

clear T
timetemp = [0:0.01:(length(adv08(:,1))-1)*0.01]'+15*60;
T = array2table(timetemp);
T.Properties.VariableNames(1) = {'time'};
writetable(T,[saveinsitu,'\adv_time.csv'])

clear i* waveg pg T1 T2
