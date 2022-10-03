% This code will plot rectified images using code from
% D_gridGenExpampleRect.m in the CIRN-Quantitative-Coastal-Imaging-Toolbox


% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\cameras'))
addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\CIRN-Quantitative-Coastal-Imaging-Toolbox\X_CoreFunctions\'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))
%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

xreg = [24 36];
yreg = [-13.5 13.5];
imrange = 10000:10999;
time_range = [datenum(0,0,0,0,0,11364/8) datenum(0,0,0,0,0,12960/8)];
vtime               = 6;  
% PLAN: LOAD TIMEX
    
% Trial info
Tinfo.Hs = 0.30;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
Tinfo.spread = 40;

%% STEP 1: Create paths, files and naming
% general path and names
datapath    = 'E:/';

Tinfo = trial_files(Tinfo);
Tinfo.cam.tstart = '08-30-2018-2129UTC';
Tinfo.cam.tdate = '08-30-2018-2119UTC';
% Stereo Reconstructions
Tinfo.cam = VBAR_GENERAL_TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

%% STEP 2: Create figure folders

datapath = 'E:\';

% figure folder
fssubfolder = datestr(date,'yy-mm-dd');
figfolder   = [datapath,'figures\meas_comp\',Tinfo.cam.trialname,'\vbar\',fssubfolder,'\'];

% make figure folders
eval(['!mkdir ',datapath,'figures\meas_comp\',Tinfo.cam.trialname]);
eval(['!mkdir ',datapath,'figures\meas_comp\',Tinfo.cam.trialname,'\','\vbar\']);
eval(['!mkdir ',figfolder])

%% Camera details
cam.time = time_range(1):datenum(0,0,0,0,0,1/Tinfo.cam.Hz):time_range(end);
cam.istart = round(time_range(1)*24*3600*8);
cam.iend = round(time_range(2)*24*3600*8);
clear *temp

%%
vbarloc = 'E:\data\processed\cameras\TRM-08-30-2018-2119UTC\frames_10000-10999\vbar\';
load([vbarloc,'vbar_fulltank_0.01res_Stack_x2m_t6s.mat'])
meanV = vStack.meanV;
meanV(vStack.QCspan<30 | vStack.prob<0 | vStack.cispan>.2) = NaN;
devVp = vStack.meanV+vStack.stdV;
devVn = vStack.meanV-vStack.stdV;
% values are "good" if:
% vStack.QCspan(j)>30 && vStack.prob(j)>0 && vStack.cispan(j)<.2

meanV(end,:,:) = NaN;

imstack = load([Tinfo.savefolder,'image_stack_x24to36_res2cm.mat']);

%% plot

figure('units','inches','position',[1 1 12 8],'color','w')

xvel = repmat(0,length(vStack.y(:,1)),length(vStack.x(1,:))); % not measuring cross-shore velocity
% permute timex
sf = 4; % scale factor for quivers

% start plotting
figure('units','inches','position',[1 1 9 8],'color','w'); % open figure
count  = vtime/2; % choose middle of time step
icount = vtime; % increase by vtime/2 because the overlap length is half of the window length

savenum = 0:1:length(meanV);

for i = 1:length(meanV(1,1,:))
    
    Iplot = imstack.IR(:,:,(i*vtime*8)+3);
    TXplot = nanmean(imstack.IR(:,:,i*vtime*8:(i+1)*vtime*8),3);
    
    ax1 = axes('Position',[0.095 0.54 0.8 0.4]);
    imagesc(imstack.Y(:,1)',imstack.X(1,:)',flipud(rot90(Iplot)));
    axis image
    hold on
    quiver(vStack.y,vStack.x,squeeze(meanV(:,:,i))*sf,xvel,'AutoScale','off','LineWidth',2,'MaxHeadSize',1,'Color','r')
    % unit quiver
    quiver(9,25,0.5*sf,0,'AutoScale','off','LineWidth',2,'MaxHeadSize',1,'Color','r')
    text(9,26,'0.5 m/s','interpreter','latex','fontsize',20,'color','w')
    text(-2,25,['time: ',num2str(count),' sec'],'interpreter','latex','fontsize',20,'color','w')
    text(-12.6,25,['Vbar'],'interpreter','latex','fontsize',20,'color','w')
    axis equal
    xlim([-13,13]);
    ylim([24,36]);
%     xlabel('alongshore~(m)','interpreter','latex','fontsize',20);
    ylabel('Cross-shore~(m)','interpreter','latex','fontsize',20);
    colormap('gray')
    h1=gca;
    set(h1,'tickdir','out','xminortick','on','yminortick','on');
    set(h1,'ticklength',1*get(h1,'ticklength'));
    set(h1,'ydir','reverse');
    set(h1,'fontsize',18);
    set(h1,'xtick',[-10:5:10],'xticklabel',{''})
    
    ax2 = axes('Position',[0.095 0.095 0.8 0.4]);
    imagesc(imstack.Y(:,1)',imstack.X(1,:)',flipud(rot90(TXplot)));
    axis image; 
    hold on
    axis equal
    text(-12.6,25,['Wave-Averaged'],'interpreter','latex','fontsize',20,'color','w')
    xlim([-13,13]);
    ylim([24,36]);
    xlabel('Alongshore~(m)','interpreter','latex','fontsize',20);
    ylabel('Cross-Shore~(m)','interpreter','latex','fontsize',20);
    colormap('gray')
    h1=gca;
    set(h1,'tickdir','out','xminortick','on','yminortick','on');
    set(h1,'ticklength',1*get(h1,'ticklength'));
    set(h1,'ydir','reverse');
    set(h1,'fontsize',18);
    set(h1,'xtick',[-10:5:10],'xticklabel',{'-10' '-5' '0' '5' '10'})
    
    Sname = [Tinfo.figfolder,'6s/vbartank_',num2str(savenum(i),'%05.f')];
    print(Sname,'-dpng')
%     pause(2)
    clf
    count = count+icount;
end


% count = 0;
% for i = cam.istart:cam.iend
%     count = count+1;
%
%     % find insitu vale
%     tcamtemp = cam.time(count);
%     [temp,iI] = nanmin(abs(cam.time(count)-inst.time));
%     utemp = u(iI,:);
%     ulftemp = uLF(iI,:);
%     vtemp = v(iI,:);
%     vlftemp = vLF(iI,:);
%
%     Ir = IR(:,:,1);
%
%
%     % plot
%     ax1 = axes('Position',[0.07 0.1 0.85 0.85]);
%     imagesc(Y(:,1)',X(1,:)',flipud(rot90(Ir)))
%     hold on
%     scatter(xyz(:,2),xyz(:,1),40,[186, 94, 242]/256,'fill','MarkerEdgeColor','k')
%     quiver(xyz(:,2),xyz(:,1),ulftemp'*sf,vlftemp'*sf,'AutoScale','off','LineWidth',2,'MaxHeadSize',4,'Color',[186, 94, 242]/256);
% %     text(425,950,'(a)','interpreter','latex','fontsize',20,'Color','w');
%     text(11.1,21,'$0.2~\mathrm{m/s}$','interpreter','latex','fontsize',16,'Color','w');
%     text(11.2,22,'$\langle \vec{u} \rangle$','interpreter','latex','fontsize',25,'Color','w');
%     quiver(11.5,20.5,0,-.2*sf,'AutoScale','off','LineWidth',1.2,'MaxHeadSize',0.7,'Color','w');
%     quiver(11.5,20.5,.2*sf,0,'AutoScale','off','LineWidth',1.2,'MaxHeadSize',0.7,'Color','w');
%
% %     imagesc(Ir',X,Y,1)
% %     scatter(sz.UVd(1,:),sz.UVd(2,:),40,'r','fill','MarkerEdgeColor','k')
% %     scatter(is.UVd(1,:),is.UVd(2,:),40,[186, 94, 242]/256,'fill','MarkerEdgeColor','k')
%     axis equal
% %     grid on
%     box on
%     ylim([19 36])
%     xlim([-13.35 13.35])
%     h1=gca;
%     set(h1,'ydir','reverse')
%     set(h1,'tickdir','out','xminortick','on','yminortick','on');
%     set(h1,'ticklength',1*get(h1,'ticklength'));
%     set(h1,'fontsize',15);
%     set(h1,'xtick',[-10:5:10],'xticklabel',{'-10' '-5' '0' '5' '10'});
%     xlabel('Alongshore (m)','interpreter','latex','fontsize',15);
%     ylabel('Cross-Shore (m)','interpreter','latex','fontsize',15);
% %     text(ax1,-13.3,24.8,'(b) Stereo Reconstruction','interpreter','latex','fontsize',15);
% %     h2 = legend('Surf Zone Array','Inner Shelf Array','interpreter','latex','fontsize',15);
% %     set(h2, 'Position', [0.72 0.85 0.19 0.064]);
%     text(ax1,7.9,18.4,['\textbf{Time Step}: ',datestr((imageno(i)/8)/(24*3600),'MM:SS.FFF')],'interpreter','latex','fontsize',15);
%
%     sname = ['innershelfarray_rectimage_vid_',num2str(imageno(i),'%05.f')];
%     print([figfolder,sname],'-dpng')
%     clf
%     clear *temp
% end
