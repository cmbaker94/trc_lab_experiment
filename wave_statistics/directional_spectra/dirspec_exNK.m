%This code will read the gridded LiDAR data and then use that infromation
%to determine the two-dimensional frequency directional spectrum

%% STEP 1: Clear All
clear all
close all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\codes\insitu'))
addpath(genpath('E:\codes\trc_lab_experiment\toolbox')) 

%% STEP 2: Load Bathymetry
A        = load('E:/data/processed/lidar/Riegl/TRC_bathymetry_estimate_line.mat');
dy       = 0.20;
yp       = -14:dy:14;
h        = 1-A.h;
[X,Y]    = meshgrid(A.xp,yp);
X        = X';
Y        = Y';
H        = repmat(h,[1 length(yp)]);

%% STEP 3: Load LiDAR
lpath    = 'E:\data\processed\lidar\Velodyne\2018-09-01-22-13-48\gridded\';
B        = load([lpath,'2018-09-01-22-13-48_Velodyne-HDL-32-Data_gridded.mat']);
z        = B.griddedData.z;
z        = permute(z,[2 1 3]);
[XL,YL]  = meshgrid(B.griddedData.xvec,B.griddedData.yvec);
XL       = XL';
YL       = YL';
hL       = griddata(X,Y,H,XL,YL);

%% STEP 4: Select Points
threshold= 5;              %This is the cutoff percentage
ipos1    = [19];
jpos1    = [35:1:55];
ipos2    = [21];
jpos2    = [23:1:63];
nanpos2  = [5;35;37;38;41];
jpos2(nanpos2) = [];
ipos3    = [23];
nanpos3  = [2;6;7;33];
jpos3    = [23:1:55];
jpos3(nanpos3) = [];

for i =1:1:length(jpos1)
    zt1    = squeeze(z(ipos1,jpos1(i),:));
    T      = 1:1:length(zt1);
    j      = find(isnan(zt1));
    if 100*length(j)/length(zt1)<threshold
        zt2    = zt1;
        T2     = T;
        zt2(j) = [];
        T2(j)  = [];
        zt1    = interp1(T2,zt2,T);
        zt1    = zt1-nanmean(zt1);
    else
       error('Too many NaNs');
    end
    z_1(i,:)   = zt1;
    x_1(i,:)   = XL(ipos1,jpos1(i));
    y_1(i,:)   = YL(ipos1,jpos1(i));
    h_1(i,:)   = hL(ipos1,jpos1(i));  
end
clear i j zt2 T2 zt1

for i =1:1:length(jpos2)
    zt1    = squeeze(z(ipos2,jpos2(i),:));
    T      = 1:1:length(zt1);
    j      = find(isnan(zt1));
    if 100*length(j)/length(zt1)<threshold
        zt2    = zt1;
        T2     = T;
        zt2(j) = [];
        T2(j)  = [];
        zt1    = interp1(T2,zt2,T);
        zt1    = zt1-nanmean(zt1);
    else
       error('Too many NaNs');
    end
    z_2(i,:)   = zt1;
    x_2(i,:)   = XL(ipos2,jpos2(i));
    y_2(i,:)   = YL(ipos2,jpos2(i));
    h_2(i,:)   = hL(ipos2,jpos2(i));  
end
clear i j zt2 T2 zt1

for i =1:1:length(jpos3)
    zt1    = squeeze(z(ipos3,jpos3(i),:));
    T      = 1:1:length(zt1);
    j      = find(isnan(zt1));
    if 100*length(j)/length(zt1)<threshold
        zt2    = zt1;
        T2     = T;
        zt2(j) = [];
        T2(j)  = [];
        zt1    = interp1(T2,zt2,T);
        zt1    = zt1-nanmean(zt1);
    else
       error('Too many NaNs');
    end
    z_3(i,:)   = zt1;
    x_3(i,:)   = XL(ipos3,jpos3(i));
    y_3(i,:)   = YL(ipos3,jpos3(i));
    h_3(i,:)   = hL(ipos3,jpos3(i));  
end
clear i j zt2 T2 zt1 z

% %% Combine Arrays
% X     = [x_1;x_2;x_3];
% Y     = [y_1;y_2;y_3];
% Z     = cat(1,z_1,z_2,z_3);
% j     = find(isnan(Z));
% Z(j)  = 0;
% H     = [h_1;h_2;h_3];
% N     = length(X);
% types = repmat(sensortypeid('n'),N,1);
% bfs   = ones(N,1);
% pos   = [X(:),Y(:),zeros(N,1)];
% Time  = B.griddedData.tFrame-B.griddedData.tFrame(1);
% Z     = Z';
% %Se    = dat2dspec([Time Z],[pos types,bfs],H,512,101); 

%% Try DIWASP
X             = [x_1(1);x_1(end);x_2(1);x_2(end);x_3(7);x_3(end)];
Y             = [y_1(1);y_1(end);y_2(1);y_2(end);y_3(7);y_3(end)];
Z             = cat(1,z_1(1,:),z_1(end,:),z_2(1,:),z_2(end,:),z_3(7,:),z_3(end,:));
SM.dirs       = [0:1:359];
SM.freqs      = [0.001:0.001:3];
SM.S          = zeros(3000,360);
EP.nfft       = 1024;
EP.dres       = 360;
EP.method     = 'IMLM';
ID.datatypes  = {'elev' 'elev' 'elev' 'elev' 'elev' 'elev'};
ID.depth      = 0.1664;
ID.data       = Z';
N             = length(X);
ID.layout     = [X';Y'; zeros(N,1)'];
ID.fs         = 10;
[Smout,EPout] = dirspec(ID,SM,EP);
dtheta        = abs(Smout.dirs(2)-Smout.dirs(1));
df            = abs(Smout.freqs(2)-Smout.freqs(1));
Sf            = nansum(Smout.S,2)*dtheta;
Sd            = nansum(Smout.S,1)*df;
fsmall        = [1:0.01:3];
fsmall4       = 10^-2*fsmall.^-4;

figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
semilogy(Smout.freqs,Sf,'k','linewidth',3);
hold on
semilogy(fsmall,fsmall4,'b','linewidth',2);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
ylim([10^-4.5 10^-1.5])
xlim([0 3]);
h1=gca;
set(h1,'fontsize',15);

subplot(212)
semilogy(Smout.dirs,Sd,'k','linewidth',3);
hold on
grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([0 360]);
h1=gca;
set(h1,'fontsize',15);



