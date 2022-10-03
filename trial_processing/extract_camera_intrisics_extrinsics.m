function extract_camera_intrisics_extrinsics(Tinfo,trialloc)
% Plot XYZ coordinates of point cloud exported from Photoscan.
% This will make a series of several plots of sea surface elevation and
% compute the sea surface elevation as a function of time. 

% extraction: point = 0, y-transect = 1, x-transect = 2, region = 3
% trialloc: if old loction = 0, if new location = 1

% trial information based on wave conditions
Tinfo = trial_files(Tinfo);
% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);
% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);
% STEP 1: Load Stereo Reconstruction Data
[camera.time] = cam_time(Tinfo);


%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

camerasys   = 'TRM';                % camera setup - e.g., TRM (offshore) or TRC (onshore)
tdate       = Tinfo.cam.tdate;%'08-30-2018-2119UTC'; % trial date and time - format ex: 09-01-2018-2155UTC
scene       = '1';                  % scene number of trial - typ 1
imagestart  = 07200;                  % images number of first frame on file 
numframes   = 4800;                   % number of frames processed
xROI        = [25 37];              % cross-shore region of interest within the tank
yROI        = [-14 14];             % alongshore region of interest within the tank
% xloc        = [31;31;31;31;31;31;31];    % x-location of timeseries
% yloc        = [-12;-8;-4;0;4;8;12];      % y-locatino of timeseries
xloctran    = [29; 30; 31];                   % cross-shore location of transect
yloctran    = 0.1;%-0.068;
tide        = 1.07;                 % still water level in tank (m)
range       = 0.05;%0.25;%[0.05 0.05];                 % range in x and y where points are averaged and std is taken
plotFlag    = 0;                    % 1 means make DEM plots, 0 means not to make DEM plots, 2 plots the figure flipped to compare with LiDAR

% extraction: point = 0, y-transect = 1, x-transect = 2, region = 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% STEP 1: Create paths, files and naming

% general path and names
if trialloc == 0
    datapath    = 'D:\TRC_Fall_Experiment\photoscan\';
elseif trialloc == 1
    datapath    = 'F:\Metashape\';%
end

trialname   = [camerasys,'-',tdate];
imagerange  = [num2str(imagestart,'%05.f'),'-',num2str(imagestart+(numframes-1),'%05.f')];
trimname    = ['frames_',imagerange,'\']; 
trimnamein    = ['frames_',num2str(imagestart,'%05.f'),'-',num2str(imagestart+9,'%05.f'),'\'];  

% intrinsics path
inpath    = [datapath,trialname,'_Scene',scene,'\',trimnamein,''];
expath    = [datapath,trialname,'_Scene',scene,'\',trimname,''];
savefolder = ['E:\data\processed\cameras\',trialname,'\',trimname];
% psname = [savefolder,'camera2_intrinsics_extrinsics.mat'];


%% load intrinsic text file
        
% Load text file
fileID      = fopen([expath,'c2_calibration','.txt'],'r');
filetemp    = textscan(fileID,'%s','Delimiter','\n');

% grab introduction information
intro = filetemp{1,1}(1:3,1);

width = filetemp{1,1}(4,1);
width = str2double(width{1,1}(8:11));

height = filetemp{1,1}(5,1);
height = str2double(height{1,1}(9:12));

f = filetemp{1,1}(6,1);
f = str2double(f{1,1}(4:end-4));

cx = filetemp{1,1}(7,1);
cx = str2double(cx{1,1}(5:end-5));

cy = filetemp{1,1}(8,1);
cy = str2double(cy{1,1}(5:end-5));

for ik = 1:3
    temp = filetemp{1,1}(8+ik,1);
    k(ik) = str2double(temp{1,1}(5:end-5));
end

for ip = 1:2
    temp = filetemp{1,1}(11+ip,1);
    p(ip) = str2double(temp{1,1}(5:end-5));
end    

b(1:2) = 0;

%% load extrinsic text file
        
% Load text file
fileID      = fopen([expath,'camera_extrinsics.txt'],'r');
filetemp    = textscan(fileID,'%s','Delimiter','\t');

X = str2double(filetemp{1,1}(20,1));
Y = str2double(filetemp{1,1}(21,1));
Z = str2double(filetemp{1,1}(22,1));
Omega = str2double(filetemp{1,1}(23,1));
Phi = str2double(filetemp{1,1}(24,1));
Kappa = str2double(filetemp{1,1}(25,1));
% r11, r12, r13, r21, r22, r23, r31, r32, r33

count = 0;
for ir1 = 1:3
    for ir2 = 1:3
        count = count+1;
        R(ir2,ir1) = str2double(filetemp{1,1}(25+count,1));
    end
end

eval(['save -v7.3 ',[savefolder,'c2_prop_Metashape.mat'],' X',' Y',' Z',' Omega',' Phi',' Kappa',' R',...
    ' width',' height',' f',' cx',' cy',' k',' p',' b'])

%% Convert to guiBathy

NU = width;
NV = height
cUo = width/2+cx;
cVo = height/2+cy;
Fx = b(1)+f;
Fy = f;
d1 = k(1);
d2 = k(2);
d3 = k(3);
t1 = p(2);
t2 = p(1);
[az,tilt,roll] = R2Angles(R');
tilt =  pi-tilt;
az   = pi+az;
roll = roll;

eval(['save -v7.3 ',[savefolder,'c2_prop_guiBathy.mat'],' X',' Y',' Z',' az',' tilt',' roll',' NU',' NV',' cUo',' cVo',' Fx',' Fy',...
    ' d1',' d2',' d3',' t1',' t2',' R'])
end

