function extract_velo_lidar(Tinfo)
% This code will extract the xyz information of the velodyne lidar
% collections from the TRC experiment.

addpath('E:\code\lidar\Velodyne')

Tinfo = trial_files(Tinfo);

%% STEP 0: Define Parameters 
vfile = [Tinfo.lidar.tdate,'_Velodyne-HDL-32-Data'];                     % File to Process
eval(['!mkdir E:\data\processed\lidar\Velodyne\',vfile(1:19)])

p.pcapFile = ['E:\data\raw\lidar\Velodyne\',vfile,'.pcap'];                          % location of file to read
p.batchFileLocation = ['E:\data\processed\lidar\Velodyne\',vfile(1:19)];   % location to store batch files
p.batchByteSize = 500000000;                                 % size of batches: 500 MB (bytes); use Inf for all data       
p.batchLimit = Inf;                                                  % number of batches to extract; use Inf for all data

%% STEP 1: Define LiDAR Rotational Matrix
% Define Rotation Parameters & Lidar Location

% NOTES:
% Guesses need to be close enough to put the basin roughly where it's
% supposed to be. roll & tilt will be refined later. Code right now assumes
% yaw is 180 and does not refine - could be added later if have at least two
% known objects in the field of view or using a fit to lab walls.

p.tilt =47.9;              % guess rotation around y
p.roll = 2.6;            	% guess rotation around x
p.yaw = 180;              	% guess rotation around z


p.x = 21.954;   %22.115;  	% x (m) location of lidar in WRL coordinates
p.y = -3.530;   %-3.757; 	% y (m) location of lidar in WRL coordinates
p.z = 7.292;    %7.25;  	% z (m) location of lidar in WRL coordinates

% Define Region of good data
p.boxExtentsX = [19,35];        % spatial extents in rotated coordinate system relative to lidar
p.boxExtentsY = [-13.5 13.5];
% p.boxExtentsZ = [-8,-5]; 
p.boxExtentsZ = [-0.75,2.25];

% Define Grid Parameters (assume Grid region is the good data region)
p.dX = 0.5;
p.dY = 0.25;

%% STEP 2: Begin processing
% read in file and save batch files as mat files in specified folder
% location
veloBatch(p)

% Define rough rotation matrix
RT = veloMakeRT(p.roll,p.tilt,p.yaw,p.x,p.y,p.z);

% Organize into individual frames, trim data to ROI, and save to disk
% NOTE:right now this saves over the original raw data read from the pcap
% file, so if you want data before rectification & trimming you need to
% change that. % CB to stop the overwrite
veloOrganize(p,RT);
%%%%%% NOTE THAT THE FILES ARE GROWING IN SIZE AT THIS STEP. SOMEWHERE IN
%%%%%% HERE IT IS ADDING THE FILES TOGETHER


% Make surfaces from each frame (DOES NO SPRAY FILTERING RIGHT NOW; VERY SIMPLE GRIDDING METHOD)
griddedData = veloMakeZGrid(p);

%% STEP 3: Make some figures
[M,N,P]=size(griddedData.z);
% Hs = ;

figure('Units','Normalized','Position',[0.1 0.3 0.8 0.5])
ax1 = subplot(1,4,1);
pcolor(griddedData.xvec,griddedData.yvec',griddedData.z(:,:,round(P/2)));shading flat;
colormap(ax1,'jet');caxis auto
c = colorbar('EastOutside');
xlabel('Cross-Shore (m)');
ylabel('Alongshore (m)');
title(['Frame #' num2str(round(P/2))]);
caxis([0.8 1.6])

ax2 = subplot(1,4,2);
pcolor(griddedData.xvec,griddedData.yvec',nanmean(griddedData.z,3));shading flat;
colormap(ax2,'jet');caxis auto
c = colorbar('EastOutside');
xlabel('Cross-Shore (m)');
title('Mean Elevation');
caxis([0.1 1.4])

ax3 = subplot(1,4,3);
pcolor(griddedData.xvec,griddedData.yvec',4*nanstd(griddedData.z,[],3));shading flat;
colormap(ax3,'jet');caxis auto
c = colorbar('EastOutside');
xlabel('Cross-Shore (m)');
title('H_s');
caxis([0 0.4])

ax4 = subplot(1,4,4);
pcolor(griddedData.xvec,griddedData.yvec',nanmean(griddedData.numPts,3));shading flat;colormap(ax4,'parula')
c = colorbar('EastOutside');
cmap=colormap(gcf);cmap(1,:)=[1 1 1];colormap(ax4,cmap);
xlabel('Cross-Shore (m)');
title('Number of Data Points per Bin');

figure('Units','Normalized','Position',[0.1 0.1 0.8 0.8]);
ax1=subplot(3,1,1);
indCenter = find(griddedData.yvec>=0,1,'first');
colors = jet(8);
for i=1:8
    plot(griddedData.xvec,griddedData.z(indCenter,:,round(P/2)+3*i),'color',colors(i,:));
    hold all
end
xlabel('X (m)')
ylabel('Z (m below scanner)')
title(['Example Free-surface Profile at Y = ' num2str(griddedData.yvec(indCenter))]);

ax2=subplot(3,1,2);
pcolor(griddedData.tFrame',griddedData.xvec,squeeze(griddedData.z(indCenter,:,:)));shading flat
colormap(ax2,'jet');caxis([-6.4 6]);c = colorbar('EastOutside');caxis auto;
xlim([griddedData.tFrame(round(P/2))-45 griddedData.tFrame(round(P/2))+45])
ylim([23 34])
xlabel('Time (s)')
ylabel('X (m)')
title(['Cross-shore Timestack at Y = ' num2str(griddedData.yvec(indCenter))]);

ax3=subplot(3,1,3);
pcolor(griddedData.yvec,griddedData.tFrame,squeeze(griddedData.z(:,find(griddedData.xvec>=28,1,'first'),:))');shading flat
colormap(ax3,'jet');caxis auto;c = colorbar('EastOutside');
ylim([griddedData.tFrame(round(P/2))-10 griddedData.tFrame(round(P/2))+10])
xlabel('Y (m)')
ylabel('Time (s)')
title('Alongshore Timestack at X = 28 m')


