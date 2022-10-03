function filter_data(Tinfo)

% Set up paths and clear workspace
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\cameras'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

%% STEP 1: Create paths, files and naming

% % general path and names
datapath    = 'E:\';
% 
Tinfo = trial_files(Tinfo);

% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

%% camera
Tinfo.cam.regy = [-14  14];
camfile = [Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_resx',num2str((Tinfo.cam.dx)*100),'cm_resy',num2str((Tinfo.cam.dx)*100),'cm.mat'];
load(camfile);

% z = medfilt3(z);
for i = 1:size(z,1)
    for j = 1:size(z,2)
        zfilt(i,j,:) = movmedian(squeeze(z(i,j,:)),3,'omitnan');
    end
end

clear z
z = zfilt;
sname = [camfile(1:end-4),'_filtered.mat'];
eval(['save -v7.3 ',sname,' x',' y',' z'])
clear x y zfilt

%% lidar

lidarfile = [datapath,'data\processed\lidar\Velodyne\',Tinfo.lidar.tdate,'\gridded\',Tinfo.lidar.tdate,'_Velodyne-HDL-32-Data_gridded.mat'];
load(lidarfile);

% for i = 1:length(griddedData.xvec)
%     z(:,i,:) = medfilt2(squeeze(griddedData.z(:,i,:)));
% end
for i = 1:size(griddedData.z,1)
    for j = 1:size(griddedData.z,2)
        zfilt(i,j,:) = movmedian(squeeze(griddedData.z(i,j,:)),3,'omitnan');
    end
end

xvec = griddedData.xvec;
yvec = griddedData.yvec;
tFrame = griddedData.tFrame;

clear griddedData
griddedData.xvec = xvec;
griddedData.yvec = yvec;
griddedData.z = zfilt;
griddedData.tFrame = tFrame;

sname = [lidarfile(1:end-4),'_filtered.mat'];
eval(['save -v7.3 ',sname,' griddedData'])
end
