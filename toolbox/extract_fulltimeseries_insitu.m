function [press,wg,camera,lidar] = extract_fulltimeseries_insitu(Tinfo,camera)
% load data and extract the remote sensing at the location of the insitu
% data

[press,wg] = load_insitu(Tinfo);

% % load insitu data
% F1          = matfile([Tinfo.datapath,'data/processed/insitu/',Tinfo.sz.tdate,'/',Tinfo.sz.tdate,'-insitu.mat']);
% pg          = F1.press;
% waveg       = F1.wg;
% Tinfo_insitu       = F1.Tinfo;
% Hz          = pg.Hz;
% maxfac = 1.2;
% offset  = 0.05;
%  
% starttemp   = datenum(Tinfo_insitu.time(1:end-3));
% endtemp     = starttemp+((length(waveg.wg1)-1)*datenum(0,0,0,0,0,1/Hz));
% timetemp    = starttemp:datenum(0,0,0,0,0,1/Hz):endtemp;
%  
% % find overlapping range for insitu
% [temp,istart] = nanmin(abs(camera.time(1)-timetemp));
% [temp,iend] = nanmin(abs(camera.time(end)-timetemp));
%  
% press.time = timetemp(istart:iend)';
% press.name  = fieldnames(pg.xyzd)';
%  
% % offset for time sync
% istart = istart+Tinfo.offsets(1);
% iend = iend+Tinfo.offsets(1);
%  
% filtramp = 1000;
% % compute bulk statistics
% for i = 1:length(press.name)
%     % load pressure data
%     eval(['press.press(:,i)   = pg.',press.name{i},'(istart:iend);'])
%     press.eta(:,i) = press2sse_timeseries(press.press(:,i),Hz,offset,maxfac,filtramp);
%     eval(['press.xyz(:,i)     = pg.xyzd.',press.name{i},'(1,1:3);']) % z-elevation of sensor in tank coordinates
%     press.loc{i} = ['p',sprintf('%02d',str2double(press.name{i}(6:end)))];
% end
%  
% clear pg *temp T1 T2
%  
% % Compute bulk statistics from wave gages
% % get wave gage names
% wg.names    = fieldnames(waveg.xyzd)';
% wg.time     = press.time;
% % compute bulk statistics
% for i = 1:length(wg.names)
%     % load wave gage data
%     eval(['wg.z(:,i)                     = waveg.',wg.names{i},'(istart:iend);'])
%     eval(['wg.xyz(:,i)              = waveg.xyzd.',wg.names{i},'(1,1:3);']) 
%     wg.loc{i} = ['wg',sprintf('%02d',str2double(wg.names{i}(3:end)))];
% end
%  
% clear i* waveg pg T1 T2
 
%% STEP 2: Load Stereo Reconstruction Data
camsource = 1;
 
if camsource == 1
    if Tinfo.filt == 1
        cam             = load([Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_res',num2str((Tinfo.cam.dx)*100),'cm_filtered.mat']); 
    else 
        cam             = load([Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_res',num2str((Tinfo.cam.dx)*100),'cm.mat']);
    end
    for i = 1:size(press.xyz,2)
        [temp,ix] = nanmin(abs(cam.x(1,:)-press.xyz(1,i)));
        [temp,iy] = nanmin(abs(cam.y(:,1)-press.xyz(2,i)));
        camera.xy(1,i) = cam.x(1,ix);
        camera.xy(2,i) = round(cam.y(iy,1),2);
        x = squeeze(cam.z(iy,ix,:));
        camera.z(:,i)  = x;
        camera.loc{i} = ['c',sprintf('%02d',str2double(press.name{i}(6:end)))];
    end
elseif camsource == 2
    F1 = matfile([Tcam.datafolder,'dem_pt_sz_array_400cm2.mat']);
    camera.xy(1,:) = F1.xloc;
    camera.xy(2,:) = F1.yloc;
    camera.z = F1.eta_pt_mean';
    % changing to match rest
end
 
clear cam T1 T2 fn*

%% STEP 3: Load VeloLidar
% load([datapath,'data\processed\lidar\Velodyne\',lidarfile(1:19),'\',lidarfile,'.mat']);
lid = lidar_data_prep(Tinfo);
 
lidar.Hz = 10;
lidar.time = lid.time;

for i = 1:size(press.xyz,2)
    [temp,ix] = nanmin(abs(lid.x-press.xyz(1,i)));
    [temp,iy] = nanmin(abs(lid.y-press.xyz(2,i)));
    lidar.xy(1,i) = lid.x(ix);
    lidar.xy(2,i) = lid.y(iy);
    lidar.z(:,i)  = squeeze(lid.z(iy,ix,:));
    lidar.loc{i} = ['l',sprintf('%02d',str2double(press.name{i}(6:end)))];
end
clear lid
 