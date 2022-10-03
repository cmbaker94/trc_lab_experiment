% This code will plot rectified images using code from
% D_gridGenExpampleRect.m in the CIRN-Quantitative-Coastal-Imaging-Toolbox


% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\cameras'))
addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\CIRN-Quantitative-Coastal-Imaging-Toolbox\X_CoreFunctions/'))
%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Trial info
    Tinfo.Hs = 0.30;
    Tinfo.Tp = 2;
    Tinfo.tide = 1.07;
    Tinfo.spread = 40;
    time_range = [datenum(0,0,0,0,0,11364/8) datenum(0,0,0,0,0,12960/8)];
    
if Tinfo.spread == 0
    % In situ
    sz.Tdate = '09-01-2018-2213UTC';
    is.Tdate = '09-06-2018-1559UTC';
    % Lidar (LI)
    lidarfile = '2018-09-01-22-13-48_Velodyne-HDL-32-Data_gridded';
    % Stereo Reconstructions (SR)
    Tcam.tstart         = '09-01-2018-2214UTC';                   % time starting collection based on spreadsheet
    Tcam.tdate          = '09-01-2018-2155UTC';      % trial date and time - format ex: 09-01-2018-2155UTC
    offsets             = [-487;-69];% [-588; -70]; % index offset relative to camera [in situ, lidar]
elseif Tinfo.spread == 20
    % In situ
    sz.Tdate = '08-30-2018-2222UTC';
    is.Tdate = '09-06-2018-1655UTC';
    % Lidar (LI)
    lidarfile = '2018-08-30-22-22-39_Velodyne-HDL-32-Data_gridded';
    % Stereo Reconstructions (SR)
    Tcam.tstart  = '08-30-2018-2226UTC';                   % time starting collection based on spreadsheet
    Tcam.tdate   = '08-30-2018-2216UTC';      % trial date and time - format ex: 09-01-2018-2155UTC 
elseif Tinfo.spread == 40
    % In situ
    sz.Tdate = '08-30-2018-2129UTC'; % 2119+20 min
    is.Tdate = '09-06-2018-1841UTC';
    isoffset = 481;
    % Lidar (LI)
    lidarfile = '2018-08-30-21-29-26_Velodyne-HDL-32-Data_gridded';
    % Stereo Reconstructions (SR)
    Tcam.tstart  = '09-06-2018-1842UTC';                   % time starting collection based on spreadsheet
    Tcam.tdate   = '09-06-2018-1836UTC';      % trial date and time - format ex: 09-01-2018-2155UTC 
end

%% STEP 1: Create paths, files and naming

Tcam = TRC_camera_info(Tcam);
imagepath = ['H:\TRC_Fall_Experiment\',Tcam.camerasys,'-',Tcam.tdate,'\',Tcam.camerasys,'-',Tcam.tdate,'_Scene1_JPEG\'];

%% STEP 2: Create figure folders

datapath = 'E:\';

% figure folder
fssubfolder = datestr(date,'yy-mm-dd');
figfolder   = [datapath,'figures\meas_comp\',Tcam.trialname,'\innershelfarray_rectimage\',fssubfolder,'\'];

% make figure folders
eval(['!mkdir ',datapath,'figures\meas_comp\',Tcam.trialname]);
eval(['!mkdir ',datapath,'figures\meas_comp\',Tcam.trialname,'\','\innershelfarray_rectimage\']);
eval(['!mkdir ',figfolder])

%% Camera details

load([datapath,'data/processed/cameras/c2_intrinsics_extrinsics.mat']);

% extrinsics = [39.35 0.02 11.03 267.7*pi/180 32.82*pi/180 0.13*pi/180];
extrinsics(1) = extrinsics(1)+0.75;

% starttemp       = datenum(Tcam.tstart(1:end-3),'mm-dd-yyyy-HHMM')+datenum(0,0,0,0,0,Tcam.imagestart/Tcam.Hz);
% endtemp         = datenum(Tcam.tstart(1:end-3),'mm-dd-yyyy-HHMM')+datenum(0,0,0,0,0,(Tcam.imagestart+Tcam.numframes-1)/Tcam.Hz);
% camera.time        = starttemp:datenum(0,0,0,0,0,1/Tcam.Hz):endtemp;

cam.time = time_range(1):datenum(0,0,0,0,0,1/Tcam.Hz):time_range(end);
cam.istart = round(time_range(1)*24*3600*8);
cam.iend = round(time_range(2)*24*3600*8);
clear *temp

%% Insitu

% load insitu data
F1          = matfile([datapath,'data/processed/insitu/',is.Tdate,'/',is.Tdate,'-insitu.mat']);
% press       = F1.press;
% wg          = F1.wg;
Tcam.is_insitu    = F1.Tinfo;
vel       = F1.vel;
ptemp   	= F1.press;
Hz          = ptemp.Hz;
clear ptemp

t_inst = 0:datenum(0,0,0,0,0,1/Hz):((length(vel.u1)-1)*datenum(0,0,0,0,0,1/Hz));
% find overlapping range for insitu
[temp,istart] = nanmin(abs(time_range(1)-t_inst));
[temp,iend] = nanmin(abs(time_range(2)-t_inst));
inst.time = t_inst(istart:iend)';
istart = istart - isoffset;
iend = iend - isoffset;


velnames  = fieldnames(vel.xyzd);

% % fix the different elevation for stacked - ONLY FOR IS
% if szvsis == 1
%     pstack(1,:) = press.xyzd.press1;
%     pstack(2,:) = press.xyzd.press2;
%     pstack(3,:) = press.xyzd.press3;
%     stackoffset(1) = 0;
%     stackoffset(2) = pstack(2,3)-pstack(1,3);
%     stackoffset(3) = pstack(3,3)-pstack(1,3);
% end

mvavgt= 8;% sec

% compute bulk statistics
for i = 1:length(velnames)
    % load pressure data
    instnum = str2double(velnames{i}(4:end));
%     display(instnum)
    eval(['u(:,instnum) = vel.u',velnames{i}(4:end),'(istart:iend);'])
    eval(['v(:,instnum) = vel.v',velnames{i}(4:end),'(istart:iend);'])
    ulf(:,instnum) = movmean(u(:,instnum),100*mvavgt)-nanmean(u(:,instnum));
    vlf(:,instnum) = movmean(v(:,instnum),100*mvavgt)-nanmean(v(:,instnum));
    uLF2(:,instnum) = pl64ta(ulf(:,instnum),mvavgt*30);
    vLF2(:,instnum) = pl64ta(vlf(:,instnum),mvavgt*30);
    uLF(:,instnum) = pl64ta(u(:,instnum),mvavgt*100);
    vLF(:,instnum) = pl64ta(v(:,instnum),mvavgt*100);
    instname{instnum}             = velnames{i};
    eval(['xyz(instnum,:) = vel.xyzd.',velnames{i},'(1,1:3);']) 
end
% clear vel

figure;
subplot(4,1,1)
plot(u)
ylim([-.5 .5])
subplot(4,1,2)
plot(uLF)
ylim([-.5 .5])
subplot(4,1,3)
plot(ulf)
ylim([-.5 .5])
subplot(4,1,4)
plot(uLF2)
ylim([-.5 .5])

figure;
subplot(4,1,1)
plot(v)
ylim([-.2 .2])
subplot(4,1,2)
plot(vLF)
ylim([-.2 .2])
subplot(4,1,3)
plot(vlf)
ylim([-.2 .2])
subplot(4,1,4)
plot(vLF2)
ylim([-.2 .2])

%% Create matrix to rectify images
% see code: D_gridGenExampleRect.m

localOrigin = [0, 0]; % [ x y]
localAngle =[0]; % Degrees +CCW from Original World X
localFlagInput=1;

ixlim=[19 36];
iylim=[-14.5 14.5];
idxdy=0.01;

iz=0;

%  World Extrinsics, need to make into sell
Extrinsics{1}=extrinsics;
Intrinsics{1}=intrinsics;

% %  Local Extrinsics
% localExtrinsics{k} = localTransformExtrinsics(localOrigin,localAngle,1,xtrinsics{1});

%  Create Equidistant Input Grid
[iX iY]=meshgrid([ixlim(1):idxdy:ixlim(2)],[iylim(1):idxdy:iylim(2)]);

%  Make Elevation Input Grid
iZ=iX*0+iz;

% If entered as Local
if localFlagInput==1
    % Assign local Grid as Input Grid
    localX=iX;
    localY=iY;
    localZ=iZ;
    
    % Assign world Grid as Rotated local Grid
    [ X Y]=localTransformEquiGrid(localOrigin,localAngle,0,iX,iY);
    Z=X*.0+iz;
end

teachingMode = 0;


%% plot
imageno =  0:1:length(dir([imagepath, 'Movie1_Scene1_c2_*']))-1;
sf = 5;

figure('units','inches','position',[1 1 12 8],'color','w')

count = 0;
for i = cam.istart:cam.iend
    count = count+1;
    
    % find insitu vale
    tcamtemp = cam.time(count);
    [temp,iI] = nanmin(abs(cam.time(count)-inst.time));
    utemp = u(iI,:);
    ulftemp = uLF(iI,:);
    vtemp = v(iI,:);
    vlftemp = vLF(iI,:);
    
    %read image
    imagefile = [imagepath,getfield(dir([imagepath,'Movie1_Scene1_c2_',sprintf('%05d',imageno(i)),'_*.jpg']),'name')];
    IM{1} = imread(imagefile);
    
    % World Rectification
    [Ir]= imageRectifier(IM,Intrinsics,Extrinsics,X,Y,Z,teachingMode);
    
    
    % plot
    ax1 = axes('Position',[0.07 0.1 0.85 0.85]);
    imagesc(Y(:,1)',X(1,:)',flipud(rot90(Ir)))
    hold on
    scatter(xyz(:,2),xyz(:,1),40,[186, 94, 242]/256,'fill','MarkerEdgeColor','k')
    quiver(xyz(:,2),xyz(:,1),ulftemp'*sf,vlftemp'*sf,'AutoScale','off','LineWidth',2,'MaxHeadSize',4,'Color',[186, 94, 242]/256);
%     text(425,950,'(a)','interpreter','latex','fontsize',20,'Color','w');
    text(11.1,21,'$0.2~\mathrm{m/s}$','interpreter','latex','fontsize',16,'Color','w');
    text(11.2,22,'$\langle \vec{u} \rangle$','interpreter','latex','fontsize',25,'Color','w');
    quiver(11.5,20.5,0,-.2*sf,'AutoScale','off','LineWidth',1.2,'MaxHeadSize',0.7,'Color','w');
    quiver(11.5,20.5,.2*sf,0,'AutoScale','off','LineWidth',1.2,'MaxHeadSize',0.7,'Color','w');

%     imagesc(Ir',X,Y,1)
%     scatter(sz.UVd(1,:),sz.UVd(2,:),40,'r','fill','MarkerEdgeColor','k')
%     scatter(is.UVd(1,:),is.UVd(2,:),40,[186, 94, 242]/256,'fill','MarkerEdgeColor','k')
    axis equal
%     grid on
    box on
    ylim([19 36])
    xlim([-13.35 13.35])
    h1=gca;
    set(h1,'ydir','reverse')
    set(h1,'tickdir','out','xminortick','on','yminortick','on');
    set(h1,'ticklength',1*get(h1,'ticklength'));
    set(h1,'fontsize',15);
    set(h1,'xtick',[-10:5:10],'xticklabel',{'-10' '-5' '0' '5' '10'});
    xlabel('Alongshore (m)','interpreter','latex','fontsize',15);
    ylabel('Cross-Shore (m)','interpreter','latex','fontsize',15);
%     text(ax1,-13.3,24.8,'(b) Stereo Reconstruction','interpreter','latex','fontsize',15);
%     h2 = legend('Surf Zone Array','Inner Shelf Array','interpreter','latex','fontsize',15);
%     set(h2, 'Position', [0.72 0.85 0.19 0.064]);
    text(ax1,7.9,18.4,['\textbf{Time Step}: ',datestr((imageno(i)/8)/(24*3600),'MM:SS.FFF')],'interpreter','latex','fontsize',15);
     
    sname = ['innershelfarray_rectimage_vid_',num2str(imageno(i),'%05.f')];
    print([figfolder,sname],'-dpng')
    clf
    clear *temp
end
