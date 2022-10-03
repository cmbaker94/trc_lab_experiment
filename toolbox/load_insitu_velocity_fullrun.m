function [vel] = load_insitu_velocity_fullrun(Tinfo,tsel)
% load and extract the right time from the timeseries.
% input: Tinfo structure, tsel (slected time)
% output: press, wg data

% load insitu data
F1          = matfile([Tinfo.datapath,'data/processed/insitu/',Tinfo.sz.tdate,'/',Tinfo.sz.tdate,'-insitu.mat']);
pg          = F1.press;
waveg       = F1.wg;
velocity    = F1.vel;
tstart      = F1.Tinfo;
Hz          = pg.Hz;
maxfac = 1.2;
offset  = 0.05;
 
starttemp   = datenum(tstart.time(1:end-3));
endtemp     = starttemp+((length(waveg.wg1)-1)*datenum(0,0,0,0,0,1/Hz));
timetemp    = starttemp:datenum(0,0,0,0,0,1/Hz):endtemp;
 
% find overlapping range for insitu
[temp,istart] = nanmin(abs(tsel(1)-timetemp));
[temp,iend] = nanmin(abs(tsel(end)-timetemp));

if istart==iend
    istart = 100*60*15;
    iend = (100*60*25)-1;
    display('This trial is for inner shelf so picking the 10 min range for in situ sz.')
end

press.time = timetemp(istart:iend)';
press.name  = fieldnames(pg.xyzd)';
 
% offset for time sync
istart = istart+Tinfo.offsets(1);
iend = iend+Tinfo.offsets(1);
 
filtramp = 1000;
g       = 9.81;
rho     = 1000;
% compute bulk statistics
% for i = 1:length(press.name)
%     % load pressure data
%     eval(['press.press(:,i)   = pg.',press.name{i},'(istart:iend);'])
%     press.eta(:,i) = press2sse_timeseries(press.press(:,i),Hz,offset,maxfac,filtramp);
%     eval(['press.xyz(:,i)     = pg.xyzd.',press.name{i},'(1,1:3);']) % z-elevation of sensor in tank coordinates
%     press.loc{i} = ['p',sprintf('%02d',str2double(press.name{i}(6:end)))];
%     press.etahyd(:,i) = (press.press(:,i)./(rho*g));
% end
 
clear pg *temp T1 T2

% wave gages
% get wave gage names
wg.names    = fieldnames(waveg.xyzd)';
wg.time     = press.time;
% compute bulk statistics
% for i = 1:length(wg.names)
%     % load wave gage data
%     eval(['wg.z(:,i)                     = waveg.',wg.names{i},'(istart:iend);'])
%     eval(['wg.xyz(:,i)              = waveg.xyzd.',wg.names{i},'(1,1:3);']) 
%     wg.loc{i} = ['wg',sprintf('%02d',str2double(wg.names{i}(3:end)))];
% end
 
clear waveg pg T1 T2

% velocity
% get vel names
vel.names    = fieldnames(velocity.xyzd)';
vel.time     = press.time;
% compute bulk statistics
for i = 1:length(vel.names)
    % load wave gage data
    eval(['vel.u(:,i)                     = velocity.u',vel.names{i}(4:end),'(istart:iend);'])
    eval(['vel.v(:,i)                     = velocity.v',vel.names{i}(4:end),'(istart:iend);'])
    eval(['vel.xyz(:,i)              = velocity.xyzd.',vel.names{i},'(1,1:3);']) 
    vel.loc{i} = ['vel',sprintf('%02d',str2double(vel.names{i}(4:end)))];
end

clear i* velocity T1 T2