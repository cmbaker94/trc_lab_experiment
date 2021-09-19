function [Hs,Tp,See,Seec,f,inst,xyz,Tinfo,data] = calc_Hs_spectra_insitu(Tinfo,time,WL,OL,range,maxfac,bathy,szvsis)
% Calc Bulk Stats Insitu
% INPUT:
% Tdate = time/date of the trial
% time  = time range to compute over (could be a long string or just 2x1) 
%           Note: should be in datenum format
% WL    = window length in seconds for pwelch
% OL    = overlap length in seconds for pwelch
% range = range of spectras to compute over
% TdateSZtrial = time of equivalent surfzone trial if this is an innershelf
%       insitu data collection
% OUTPUT:
% Hs    = vector of Hs from all inst
% Tp    = vector of Tp from all inst
% xy    = locations of instruments
% See   = matrix of spectra for each instrument
% f     = matrix of frequencies for each instrument (should be the same...)
% inst  = instrument names (just for reference

rho     = 1000;

%% STEP 1: Create paths, files and naming

% path to insitu data
datapath    = 'E:\data\processed\insitu\';
% timecam     = [Tinfo.cam.timevec(1) Tinfo.cam.timevec(end)+datenum(0,0,0,0,0,1/Tinfo.cam.Hz)];

% load insitu data
if szvsis == 0
    Tdate = Tinfo.sz.tdate;
elseif szvsis == 1
    Tdate = Tinfo.is.tdate;
end

F1          = matfile([datapath,Tdate,'\',Tdate,'-insitu.mat']);
press       = F1.press;
wg          = F1.wg;
Tinst       = F1.Tinfo;
Hz          = press.Hz;

%% STEP 2: Get time range

% if this is a innershelf trial it won't have the same time range so use
% the TdateSZtrial

starttemp   = datenum(Tinst.time(1:end-3));
endtemp     = starttemp+((length(wg.wg1)-1)*datenum(0,0,0,0,0,1/Hz));
t_inst       = starttemp:datenum(0,0,0,0,0,1/Hz):endtemp;

starttime  = datenum(Tinst.time(1:end-3))+time(1);
endtime    = datenum(Tinst.time(1:end-3))+time(2);
% t_inst 	= starttemp:datenum(0,0,0,0,0,1/Hz):endtemp;

% if szvsis == 0%isempty(TdateSZtrial)
%     starttemp   = datenum(Tinst.time(1:end-3));
%     endtemp     = starttemp+((length(wg.wg1)-1)*datenum(0,0,0,0,0,1/Hz));
%     t_inst       = starttemp:datenum(0,0,0,0,0,1/Hz):endtemp;
% elseif szvsis == 1
%     F2          = matfile([datapath,TdateSZtrial,'/',TdateSZtrial,'-insitu.mat']);
%     szTinst       = F2.Tinfo;
%     display('Using time of surf zone trial because this is an insitu trial')
%     starttemp   = datenum(szTinst.time(1:end-3));
%     endtemp     = starttemp+((length(wg.wg1)-1)*datenum(0,0,0,0,0,1/Hz));
%     t_inst       = starttemp:datenum(0,0,0,0,0,1/Hz):endtemp;
%     clear szT*
% elseif szvsis == 2
%     starttemp   = datenum(Tinst.time(1:end-3));
%     endtemp     = starttemp+((length(wg.wg1)-1)*datenum(0,0,0,0,0,1/Hz));
%     t_inst       = starttemp:datenum(0,0,0,0,0,1/Hz):endtemp;
% end

if isempty(time)
    % for now, use the outer surfzone as the same region
    istart  = 100000;%1;%
    iend    = 199999;%%length(wg.wg1)-1
    display(['Choosing time from ',num2str(istart),' to ',num2str(iend)])
else
        % find overlapping range for insitu
        [temp,istart] = nanmin(abs(starttime-t_inst));
        [temp,iend] = nanmin(abs(endtime-t_inst));
%         [temp,istart] = nanmin(abs(time(1)-t_inst));
%         [temp,iend] = nanmin(abs(time(end)-t_inst));
end

data.time = t_inst(istart:iend);
% data.time = t_inst;
%% Compute bulk stat with pressure gages

% depth of instruments
% NOTE STILL NEED TO FIND DEPTH AT THE XYZ LOCATIONS OF PRESSURE GAGES
% This can be based on the x-location of the instrument d.
% JK FOR NOW...
% Just giving an offset of the instrument below

% get file pressure gauge names
pnames  = fieldnames(press.xyzd);

% % fix the different elevation for stacked - ONLY FOR IS
% if szvsis == 1
%     pstack(1,:) = press.xyzd.press1;
%     pstack(2,:) = press.xyzd.press2;
%     pstack(3,:) = press.xyzd.press3;
%     stackoffset(1) = 0;
%     stackoffset(2) = pstack(2,3)-pstack(1,3);
%     stackoffset(3) = pstack(3,3)-pstack(1,3);
% end
    
% compute bulk statistics
for i = 1:length(pnames)
    % load pressure data
    eval(['ptemp = press.',pnames{i},'(istart:iend);'])
    % define geometry
%     offset = 0.05;
    eval(['x = press.xyzd.',pnames{i},'(1,1);']) 
    eval(['z = press.xyzd.',pnames{i},'(1,3);']) % z-elevation of sensor in tank coordinates
%     z = -(Tinfo.stillwat-z); % depth of sensor relative to the surface (positive upwards)
    % define depth as X amount below the elevation of the sensor
%     h       = -z+offset;
    [temp,ibathy] = nanmin(abs(x-bathy.xp));
    offset = z-bathy.h(ibathy);
    offset
    
%     % fix for stack - ONLY FOR PRESS on IS
%     if szvsis == 1
%         if length(pnames{i}) == 6
%             if pnames{i} == 'press2'
%                 offset = offset+stackoffset(2);
%             elseif pnames{i} == 'press3'
%                 offset = offset+stackoffset(3);
%             end
%         end
%     end
    % Compute spectra and bulk statistics
    [See(i,:),f(i,:),Seec(i,:,:)]   = press2spectra(ptemp,Hz,WL*Hz,OL*Hz,rho,offset,maxfac);
    
    [Hs(i),Tp(i)]       = spec2HsTp(f(i,:),See(i,:),range);
    inst{i}             = pnames{i};
    eval(['xyz(i,:) = press.xyzd.',pnames{i},'(1,1:3);']) 
    eval(['data.',pnames{i},' = ptemp;']) 
    clear ptemp z
end
clear press ptemp pnames z 

count = i;

%% Compute bulk statistics from wave gages

% get wave gage names
wgnames  = fieldnames(wg.xyzd);
% compute bulk statistics
for i = 1:length(wgnames)
    count = count+1;
    % load wave gage data
    eval(['wgtemp               = wg.',wgnames{i},'(istart:iend);'])
    % Compute spectra and bulk statistics
    [See(count,:),f(count,:),Seec(count,:,:)]   = pwelch(wgtemp,WL*Hz,OL*Hz,[],Hz,'ConfidenceLevel',0.95);
    % compute bulk statistics
    [Hs(count),Tp(count)]       = spec2HsTp(f(count,:),See(count,:),range);
    inst{count}                 = wgnames{i};
    eval(['xyz(count,:)             = wg.xyzd.',wgnames{i},'(1,1:3);']) 
    eval(['data.',wgnames{i},' = wgtemp;']) 
    clear wgtemp 
end




