% function calc_directional_spec_iter_lidar(Tinfo)
% %This code will read the gridded LiDAR data and then use that infromation
% %to determine the two-dimensional frequency directional spectrum

% %% STEP 1: Clear All
clear all
close all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
% addpath(genpath('E:\codes\insitu'))
addpath(genpath('E:\codes\trc_lab_experiment\toolbox')) 
% 
% %% Input trial details
% 
Tinfo.spread = 40;
Tinfo.Hs = 0.25;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
Tinfo.filt = 1;
instnum = 8;
iternum = 20;%120;
inHz = 10;

%% SM defaults

SM.dirs       = [0:1:359];
SM.freqs      = [0.01:0.01:3];
SM.S          = zeros(300,360);
EP.nfft       = inHz*(2^5); %  32 sec window
EP.dres       = 360;
EP.method     = 'EMEP';
IDdef.datatypes  = {'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev'};% 'elev' 'elev' 'elev' 'elev' 'elev' 'elev' };
IDdef.fs         = inHz;
% IDdef.depth      = 0.1664;

%% Get trail details

Tinfo = trial_files(Tinfo);

% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);
[camera.time] = cam_time(Tinfo);

dx = 0.5;
dy = 0.25;
% xcam    = [27.5:dx:33];
% ycam    = [-10:dy:4];%10
xlid    = [27.5:dx:33];
ylid    = [-10:dy:4];

%% STEP 2: Load Bathymetry
A        = load('E:/data/processed/lidar/Riegl/TRC_bathymetry_estimate_line.mat');
% dy       = dy;
% yp       = ylid;
h        = Tinfo.tide-A.h;
X = A.xp;
% [X,Y]    = meshgrid(A.xp,yp);
% X        = X';
% Y        = Y';
% H        = repmat(h,[1 length(yp)]);
[~,ihmin] = min(abs(xlid(1)-X));
[~,ihmax] = min(abs(xlid(2)-X));
H = h(ihmin:ihmax);
IDdef.depth = mean(H);

%% Load lidar data

% if Tinfo.filt == 1
%     sname = [Tinfo.savefolder,'lidar_grid_timeseries_regx',num2str(round(xlid(1))),'_',num2str(round(xlid(end))),'_regyneg',num2str(abs(round(ylid(1)))),'_',num2str(round(ylid(end))),'_dx',num2str(dx*100),'_dy',num2str(dy*100),'_filt'];
% else
%     sname = [Tinfo.savefolder,'lidar_grid_timeseries_regx',num2str(round(xlid(1))),'_',num2str(round(xlid(end))),'_regyneg',num2str(abs(round(ylid(1)))),'_',num2str(round(ylid(end))),'_dx',num2str(dx*100),'_dy',num2str(dy*100)];
% end

if Tinfo.filt == 1
    load(['E:\data\processed\lidar\Velodyne\',Tinfo.lidar.tdate,'\gridded\',Tinfo.lidar.tdate,'_Velodyne-HDL-32-Data_gridded_filtered.mat']);
else
    load(['E:\data\processed\lidar\Velodyne\',Tinfo.lidar.tdate,'\gridded\',Tinfo.lidar.tdate,'_Velodyne-HDL-32-Data_gridded.mat']);
end

% lidar record start and end time
starttemp           = datenum(Tinfo.lidar.tdate,'yyyy-mm-dd-HH-MM-SS');
endtemp             = datenum(Tinfo.lidar.tdate,'yyyy-mm-dd-HH-MM-SS')+(length(squeeze(griddedData.z(1,1,:)))*datenum(0,0,0,0,0,1/inHz));
griddedData.time    = starttemp:datenum(0,0,0,0,0,1/inHz):endtemp;

% index of analysis start and end 
[temp,istart]       = nanmin(abs(Tinfo.cam.timevec(1)-griddedData.time));
[temp,iend]         = nanmin(abs(Tinfo.cam.timevec(end)-griddedData.time));

% extract time series of interest
z       = squeeze(griddedData.z(:,:,istart:iend));
z        = permute(z,[2 1 3]);
XL       = griddedData.xvec;
YL       = griddedData.yvec;
XL       = repmat(XL',1,length(YL));
YL       = repmat(YL,size(XL,1),1);
% hL       = griddata(X,Y,H,XL,YL);

numPts = size(z,3)-sum(isnan(z),3);
numPts(numPts==0) = NaN;
numPtsXL = XL;
numPtsYL = YL;

ix = ismember(round(XL(:,1),2),round(xlid,2));
iy = ismember(round(YL(1,:),2),round(ylid,2));
XL = XL(ix,iy);
YL = YL(ix,iy);
z = squeeze(z(ix,iy,:));
clear ix iy griddedData

%% STEP 4: Select Points
threshold = 10;
% code to find good std min
% for it = 1:300
%     ix = randi([1, size(XL,1)],instnum,1);
%     iy = randi([1, size(YL,2)],instnum,1);
%     xd = XL(ix,1);
%     yd = YL(1,iy);
%     testx(it) = nanstd(xd);
%     testy(it) = nanstd(yd);
% end
count = 1;

for it = 1:100000
    isel = randi(length(XL(:)),instnum,1);
    for il = 1:instnum
        [~,ix(il)] = min(abs(XL(:,1)-XL(isel(il)))); 
        [~,iy(il)] = min(abs(YL(1,:)-YL(isel(il)))); 
    end
%     ix = randi([1, size(XL,1)],instnum,1);
%     iy = randi([1, size(YL,2)],instnum,1); 
    xd = XL(ix,1);
    yd = YL(1,iy);
    stdmet = 1;
    linemet = 1;
    lagmet = 1;
    if length(unique(pdist(xd),'stable')) == length(pdist(xd))
    elseif length(unique(pdist(yd),'stable')) == length(yd)
    else
        lagmet = 0;
        display('lag not met')
    end
    if length(unique(round(xd))) < 2
        display('all in one line')
        linemet = 0;
    end
    if std(xd)<1.4 || std(yd)<3.5
        display('not spread enough');
        stdmet = 0;
    end
    if stdmet == 1 && linemet == 1
        for i = 1:instnum
            zt1    = squeeze(z(ix(i),iy(i),:));
            T      = 1:1:length(zt1);
            j      = find(isnan(zt1));
            nanmet(i) = 1;
            if 100*length(j)/length(zt1)<threshold
                zt2    = zt1;
                T2     = T;
                zt2(j) = [];
                T2(j)  = [];
                zt1    = interp1(T2,zt2,T);
                zt1     = fillmissing(zt1,'linear');
                zt1    = zt1-nanmean(zt1);
            else
                display('Too many NaNs');
                nanmet(i) = 0;
            end
            zd(i,:) = zt1;
        end
        if min(nanmet) == 1
            ID              = IDdef;
            
            ID.data       = zd'; % 6000 x 14
            N             = length(xd);
            ID.layout     = [xd';yd; zeros(N,1)']; % 3 x 14
            options = {'PLOTTYPE',0};
            [Smout,EPout] = dirspec(ID,SM,EP,options);
            
            Spec(:,:,count) = Smout.S;
            xSpec(:,count) = xd;
            ySpec(:,count) = yd;
            count = count + 1;
            count
        end
    end
    if count == iternum+1
        break
    end
end

Savg = mean(Spec,3);
dirs = Smout.dirs;
freqs = Smout.freqs;
if Tinfo.filt == 1
    specsave = [Tinfo.savefolder,'dirspec_lid_regx',num2str(round(xlid(1))),'_',num2str(round(xlid(end))),'_regyneg',num2str(abs(round(ylid(1)))),'_',num2str(round(ylid(end))),'_dx',num2str(dx*100),'_dy',num2str(dy*100),'_',num2str(instnum),'inst_',num2str(iternum),'iter_filt'];
end
eval(['save -v7.3 ',specsave,' dirs',' freqs',' Savg',' Spec'])

%% Try DIWASP

xoff = 28;
xon = xlid(end);
xmid = round((xon-xoff)/2,1)+xoff;
ixoff = find(round(XL(:,1),2)==xoff);
ixon = find(round(XL(:,1),2)==xon);
ixmid = find(round(XL(:,1),2)==xmid);
iycent = find(round(YL(1,:),2)==0);
zoff = squeeze(z(ixoff,iycent,:));
zon = squeeze(z(ixon,iycent,:))-nanmean(squeeze(z(ixon,iycent,:)));
zmid = squeeze(z(ixmid,iycent,:))-nanmean(squeeze(z(ixmid,iycent,:)));
zoff = interpNaNseries(zoff,[1:size(z,3)]);
[Soff,foff,Sctemp]   = pwelch(zoff-nanmean(zoff),256,128,[],10,'ConfidenceLevel',0.95); % compute spectra
zon = interpNaNseries(zon,[1:size(z,3)])
[Son,fon,Sctemp]   = pwelch(zon(~isnan(zon))-nanmean(zon),256,128,[],10,'ConfidenceLevel',0.95); % compute spectra
zmid = interpNaNseries(zmid,[1:size(z,3)]);
[Smid,fmid,Sctemp]   = pwelch(zmid-nanmean(zmid),256,128,[],10,'ConfidenceLevel',0.95); % compute spectra

dtheta        = abs(dirs(2)-dirs(1));
df            = abs(freqs(2)-freqs(1));
Sf            = nansum(Savg,2)*dtheta;
Sd            = nansum(Savg,1)*df;
fsmall        = [1:0.01:3];
fsmall4       = 10^-2*fsmall.^-4;

Sname1 = [Tinfo.figfolder,'dirspec_lidar_regx',num2str(round(xlid(1))),'_',num2str(round(xlid(end))),'_regyneg',num2str(abs(round(ylid(1)))),'_',num2str(round(ylid(end))),'_dx',num2str(dx*100),'_dy',num2str(dy*100),'_filt_inst',num2str(instnum),'_iter',num2str(iternum),'_',EP.method];
v = VideoWriter([Sname1,'.avi']);
v.FrameRate=0.5;
open(v)



figure('units','inches','position',[1 1 14 8],'color','w');

ax1 = axes('Position',[0.1 0.55 0.45 0.35]);
ax2 = axes('Position',[0.1 0.1 0.45 0.35]);
ax3 = axes('Position',[0.65 0.1 0.3 0.8]);

axes(ax1)
plot(ax1,freqs,Sf,'k','linewidth',3);
hold on
plot(ax1,fsmall,fsmall4,'Color',[0.5 0.5 0.5],'linewidth',2);
plot(ax1,freqs,nansum(Spec(:,:,1),2)*dtheta,'b','linewidth',1.5);
plot(ax1,foff,Soff,'r','linewidth',1.5);
plot(ax1,foff,Son,'r','linewidth',1.5);
plot(ax1,fmid,Smid,'r','linewidth',1.5);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
ylim([10^-4.5 10^-1.5])
xlim([0 3]);
h1=gca;
set(h1, 'YScale', 'log')
set(h1,'fontsize',15);
title(ax1,[EP.method,', interation #1'],'interpreter','latex')

axes(ax2)
plot(ax2,dirs,Sd,'k','linewidth',3);
hold on
plot(ax2,dirs,nansum(Spec(:,:,1),1)*df,'b','linewidth',1.5);
grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([0 360]);
ylim([10^(-6.5) 10^(-4)])
h1=gca;
set(h1, 'YScale', 'log')
set(h1,'fontsize',15);

axes(ax3)
pcolor(numPtsXL,numPtsYL,numPts/size(z,3))
hold on
plot([xlid(1) xlid(end) xlid(end) xlid(1) xlid(1)],[ylid(1) ylid(1) ylid(end) ylid(end) ylid(1)],...
    'LineStyle','-.','LineWidth',1,'Color','b')
scatter(xSpec(:,1),ySpec(:,1),40,'fill','MarkerFaceColor','b','MarkerEdgeColor','k','LineWidth',0.8)
scatter(xoff,0,80,'x','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',1.5)
scatter(xon,0,80,'x','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',1.5)
scatter(xmid,0,80,'x','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',1.5)
shading interp
colorbar
caxis([.60 1])
grid
xlabel('$x$ (m)','interpreter','latex')
ylabel('$y$ (m)','interpreter','latex')
% xlim([0 360]);
h1  = gca;
% set(h1,'ColorScale','log')
set(h1,'fontsize',15);
writeVideo(v,getframe(gcf))

for i = 2:size(Spec,3)
    axes(ax1)
    semilogy(ax1,freqs,nansum(Spec(:,:,i-1),2)*dtheta,'Color',[0.7 0.7 0.7],'linewidth',1.5);
    plot(ax1,foff,Soff,'r','linewidth',1.5);
    plot(ax1,foff,Son,'r','linewidth',1.5);
    plot(ax1,fmid,Smid,'r','linewidth',1.5);
    semilogy(ax1,freqs,Sf,'k','linewidth',3);
    semilogy(ax1,freqs,nansum(Spec(:,:,i),2)*dtheta,'b','linewidth',1.5);
    ylim([10^-4.5 10^-1.5])
    xlim([0 3]);
    title(ax1,[EP.method,', interation #',num2str(i)],'interpreter','latex')

    axes(ax2)
    semilogy(ax2,dirs,nansum(Spec(:,:,i-1),1)*df,'Color',[0.7 0.7 0.7],'linewidth',1.5);
    semilogy(ax2,dirs,Sd,'k','linewidth',3);
    semilogy(ax2,dirs,nansum(Spec(:,:,i),1)*df,'b','linewidth',1.5);
    xlim([0 360]);
    ylim([10^(-6.5) 10^(-4)])
    
    axes(ax3)
    pcolor(numPtsXL,numPtsYL,numPts/size(z,3))
    hold on
    plot([xlid(1) xlid(end) xlid(end) xlid(1) xlid(1)],[ylid(1) ylid(1) ylid(end) ylid(end) ylid(1)],...
        'LineStyle','-.','LineWidth',1,'Color','b')
    scatter(xoff,0,80,'x','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',1.5)
    scatter(xon,0,80,'x','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',1.5)
    scatter(xmid,0,80,'x','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',1.5)
    scatter(xSpec(:,i),ySpec(:,i),40,'fill','MarkerFaceColor','b','MarkerEdgeColor','k','LineWidth',0.8)
    shading interp
    colorbar
    caxis([0.6 1])
    grid
    xlabel('$x$ (m)','interpreter','latex')
    ylabel('$y$ (m)','interpreter','latex')
    h1  = gca;
    set(h1,'fontsize',15);
    
    writeVideo(v,getframe(gcf))
end
close(v)

%%
dtheta        = abs(dirs(2)-dirs(1));
df            = abs(freqs(2)-freqs(1));
Sf            = nansum(Savg,2)*dtheta;
Sd            = nansum(Savg,1)*df;
fsmall        = [1:0.01:3];
fsmall4       = 10^-2*fsmall.^-4;

figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
for i = 1:size(Spec,3)
    semilogy(Smout.freqs,nansum(Spec(:,:,i),2)*dtheta,'m','linewidth',1);
    hold on
end
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
title(EP.method,'interpreter','latex')

subplot(212)
for i = 1:size(Spec,3)
    semilogy(Smout.dirs,nansum(Spec(:,:,i),1)*df,'m','linewidth',1);
    hold on
end
semilogy(Smout.dirs,Sd,'k','linewidth',3);
hold on
grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([0 360]);
h1=gca;
set(h1,'fontsize',15);

Sname1 = [Tinfo.figfolder,'dirspec_lidar_regx',num2str(round(xlid(1))),'_',num2str(round(xlid(end))),'_regyneg',num2str(abs(round(ylid(1)))),'_',num2str(round(ylid(end))),'_dx',num2str(dx*100),'_dy',num2str(dy*100),'_filt_inst',num2str(instnum),'_iter',num2str(iternum),'_',EP.method];
print(Sname1,'-dpng')


