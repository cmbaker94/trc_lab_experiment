function calc_wave_crest_length_image(Tinfo)
% This code will plot rectified images using code from
% D_gridGenExpampleRect.m in the CIRN-Quantitative-Coastal-Imaging-Toolbox


% Set up paths and clear workspace
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\cameras'))
addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\CIRN-Quantitative-Coastal-Imaging-Toolbox/X_CoreFunctions/'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

xloc = [29 30 31];

for ji = 1:length(xloc)
        % Trial info
    xshore = xloc(ji);

%% STEP 1: Create paths, files and naming
 % general path and names
    datapath    = 'E:/';
    
    % 
    Tinfo = trial_files(Tinfo);
%     if Tinfo.spread == 40
%         Tinfo.cam.tdate         = '08-30-2018-2119UTC';  
%     end
    % Stereo Reconstructions
    Tinfo.cam = TRC_camera_info(Tinfo.cam);
    
    % Data and figure storage
    [Tinfo] = wc_comp_store(Tinfo);


if Tinfo.spread == 0
    if Tinfo.Hs == 0.3
        imagepath = ['D:\TRC_Fall_Experiment\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1_JPEG\c2\'];
    else
        imagepath = ['H:\TRC_Fall_Experiment\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1_JPEG\'];
    end
elseif Tinfo.spread == 40
    if Tinfo.Hs == 0.25
        imagepath = ['D:\TRC_Fall_Experiment\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1_JPEG\'];
    else
        imagepath = ['H:\TRC_Fall_Experiment\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1_JPEG\'];
    end
else
    imagepath = ['H:\TRC_Fall_Experiment\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'\',Tinfo.cam.camerasys,'-',Tinfo.cam.tdate,'_Scene1_JPEG\'];
end

%% Establish timeseries range:
[camera.time] = cam_time(Tinfo);

%% Insitu

load([datapath,'data/processed/cameras/c2_intrinsics_extrinsics.mat']);

% extrinsics = [39.35 0.02 11.03 267.7*pi/180 32.82*pi/180 0.13*pi/180];
% extrinsics(1) = extrinsics(1);

%% Create matrix to rectify images
% see code: D_gridGenExampleRect.m

localOrigin = [0, 0]; % [ x y]
localAngle =[0]; % Degrees +CCW from Original World X
localFlagInput=1;

idxdy=0.01;
ixlim=[xshore-idxdy xshore+idxdy];%[19 36];%
iylim=[-14.5 14.5];

iz=Tinfo.tide;

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

%% extract x = 28

%% plot
imageno =  7200:1:11999;

% figure('units','inches','position',[1 1 12 8],'color','w')
for i = 1:length(imageno) 


    %read image
    imagefile = [imagepath,getfield(dir([imagepath,'Movie1_Scene1_c2_',sprintf('%05d',imageno(i)),'_*.jpg']),'name')];
    IM{1} = imread(imagefile);
    
    % World Rectification
    [Ir]= imageRectifier(IM,Intrinsics,Extrinsics,X,Y,Z,teachingMode);
    
    [val,id] = min(abs(X(1,:)-xshore));
    Irtemp = double(rgb2gray(Ir));
    IMtran(i,:) = Irtemp(:,id)';
    IMy(i,:) = Y(:,1);
    dy = Y(2,1)-Y(1,1);
    
    WL = length(IMy(i,:));
    OL = 0;
    
    speccompute = (IMtran(i,:)-nanmean(IMtran(i,:)))/256;
    speccompute             = detrend(speccompute,1);% detrending
    [Stemp,ktemp]   = pwelch(speccompute,WL,OL,[],1/dy,'ConfidenceLevel',0.95); % compute spectra
    [val,id] = min(abs(ktemp-1));
    filter = zeros(size(ktemp));
    filter(1:id-1) = 1;
    filter(id:id+2) = [0.7089 0.2363 0.0225];
    S(i,:) = Stemp.*filter;
    k(i,:) = ktemp;
    
%     % plot
%     if i<20
%         clf
%         plot(IMtran(i,:),'k')
%         hold on
%         plot(Irtemp(:,end-id),'r');
%     end
   
%    figure 
  
%    if i<20
%             clf
%             semilogy(k(i,:),S(i,:),'k','LineWidth',2)
%             xlabel('$L^{-1}$ (m$^{-1}$)','interpreter','latex','fontsize',20);
%             ylabel('$S_{\eta\eta}$ (m$^2$/Hz)','interpreter','latex','fontsize',20);
%             h1=gca;
%             set(h1, 'YScale', 'log')
%             set(h1,'tickdir','in','xminortick','on','yminortick','on');
%             set(h1,'ticklength',1*get(h1,'ticklength'));
%             set(h1,'fontsize',15);
%             %         xlim([0 1/(cam.y(end,1)-cam.y(1,1))])
%             ylim([10^-6 10^-3])
%             pause(0.2)
%         end

end

    %%
    Savg=nanmean(S,1);
%     Scavg = nanmean(Sc,1);
    
    figure
    loglog(k(6,:),Savg,'k','LineWidth',2)
    % semilogy(k(6,:),Savg,'k','LineWidth',2)
    xlabel('$L^{-1}$ (m$^-1$)','interpreter','latex','fontsize',20);
    ylabel('$S_{\eta\eta}$ (m$^2$/Hz)','interpreter','latex','fontsize',20);
    h1=gca;
    % set(h1, 'YScale', 'log')
    set(h1,'tickdir','in','xminortick','on','yminortick','on');
    set(h1,'ticklength',1*get(h1,'ticklength'));
    set(h1,'fontsize',15);
    % xlim([0 1/(cam.y(end,1)-cam.y(1,1))])
    ylim([10^-6 10^-2])
%     title(['x = ',num2str(xshore),'m , avg'],'interpreter','latex','fontsize',20);
    
    
    y = nanmean(IMy,1);
    x = xshore;
    
    
    subname = '';%'onewindow_11movmean_detrended';
%     eval(['!mkdir ',Tcam.datafolder,'AMATH'])
    psname = [Tinfo.cam.datafolder,'image_values_transect_wavecrest_length_x',num2str(xshore),subname,'.mat'];
    eval(['save -v7.3 ',psname,' IMtran',' IMy',' x']);
    psname = [Tinfo.cam.datafolder,'image_transect_wavecrest_length_x',num2str(xshore),subname,'.mat'];
    eval(['save -v7.3 ',psname,' S',' k',' Savg',' y',' x']);
    
    clear IM Ir S k Savg y x IMtran IMy x
    
end
end
