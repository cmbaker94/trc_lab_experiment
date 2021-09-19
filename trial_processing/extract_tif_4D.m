function extract_tif_4D(Tinfo,trialloc,extraction)
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
yROI        = [-13 13];             % alongshore region of interest within the tank
% xloc        = [31;31;31;31;31;31;31];    % x-location of timeseries
% yloc        = [-12;-8;-4;0;4;8;12];      % y-locatino of timeseries
xloctran    = [29; 30; 31];                   % cross-shore location of transect
yloctran    = 0.1;%-0.068;
tide        = 1.07;                 % still water level in tank (m)
range       = 0.05;%0.25;%[0.05 0.05];                 % range in x and y where points are averaged and std is taken
plotFlag    = 0;                    % 1 means make DEM plots, 0 means not to make DEM plots, 2 plots the figure flipped to compare with LiDAR

% extraction: point = 0, y-transect = 1, x-transect = 2, region = 3

extractPT           = 0;                    % 1 means extract the point values, 0 means don't
extractTRANSECTy    = 0;                % 1 means extract alongshore transect, 0 means don't
extractTRANSECTx    = 0;
extractREGION       = 0;

if extraction == 0
    extractPT   = 1;                    % 1 means extract the point values, 0 means don't
elseif extraction == 1
    extractTRANSECTy = 1;                % 1 means extract alongshore transect, 0 means don't
elseif extraction == 2
    extractTRANSECTx = 1;
elseif extraction == 3
    extractREGION   = 1;                % 1 means extract a specific region
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subname = '';

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

% Raw tiff path
tiffpath    = [datapath,trialname,'_Scene',scene,'\',trimname,'dems\'];

%
eval(['!mkdir E:\data\processed\cameras\',trialname]);
eval(['!mkdir E:\data\processed\cameras\',trialname,'\',trimname]);
savefolder = ['E:\data\processed\cameras\',trialname,'\',trimname];

% figure folder
fssubfolder = datestr(date,'yy-mm-dd');
figfolder   = ['E:\figures\cameras\images\',trialname,'\',trimname,'DEM\',fssubfolder,'\'];

% make figure folders
eval(['!mkdir E:\figures\cameras\images\',trialname]);
eval(['!mkdir E:\figures\cameras\images\',trialname,'\',trimname]);
eval(['!mkdir E:\figures\cameras\images\',trialname,'\',trimname,'DEM']);
eval(['!mkdir ',figfolder])

%% STEP 2: general details for processing

imageno     = 1:(numframes+1);
% ftype       = 1; % 1 is for DEM, 2 is dense point cloud
if extractPT == 1
    locxyz = readtable("E:\data\processed\insitu\surfzone_pressure_gages_xyz.csv");
    locxyz = locxyz(:,2:end);
    sz.xyz =  table2array(locxyz);
    inst = convertCharsToStrings(locxyz.Properties.VariableNames);
    xloc = sz.xyz(1,:);
    yloc = sz.xyz(2,:);
elseif extractREGION == 1
    dx      = 0.05;
%     display('dx,dy = 0.1, CHANGE MAY BE NECESSARY, was 0.05')
    dy      = 0.25;
    regx    = xROI(1):dx:xROI(2);
    regy    = yROI(1):dy:yROI(2);
    
%     [X,Y] = meshgrid(xtemp,ytemp,0.5);
    x       = repmat(regx,length(regy),1);
    y       = repmat(regy',1,length(regx));
    z       = NaN(length(regy),length(regx),length(imageno)-1);
    numPts  = NaN(length(regy),length(regx),length(imageno)-1);
    t       = (1:length(imageno)-1)/8;
end

    
if plotFlag == 1
    figure('units','inches','position',[1 1 12 8],'color','w'); % open figure
elseif plotFlag == 2
    figure('units','inches','position',[1 1 8 12],'color','w'); % open figure
end

%% STEP 3: Read DEM and create images

starthere = 1;

for i = starthere:length(imageno)-1
    display(['Image ',num2str(imageno(i))])
    tifffile    = [camerasys,'-',tdate,'_frames_',imagerange,'_DEM_',num2str(imageno(i)),'.tif'];
%     tifffile    = ['test_',num2str(imageno(i)),'.tif'];%[camerasys,'-',tdate,'_frames_',imagerange,'_dem_',num2str(imageno(i)),'.tif'];
    nm          = num2str(imageno(i));
    
    % read in point cloud
    filename    = [tiffpath,tifffile];%[datapath,demfile];
    [xtemp,ytemp,ztemp,resxtemp,resytemp] = read_tiff(filename);
    ytemp = fliplr(ytemp);
    ztemp = flipud(ztemp);
    
    % need to make a 2d for finding nans
    ytemp = repmat(ytemp',1,length(xtemp));
    xtemp = repmat(xtemp,length(ytemp),1);
    
    if extractPT == 1
        for k = 1:length(xloc)
            %         dsx = find(x(:,i)>(xloc(k)-range) & x(:,i)<(xloc(k)+range));
            %         dsy = find(y(:,i)>(yloc(k)-range) & y(:,i)<(yloc(k)+range));
%             dsx = xtemp>(xloc(k)-range) & xtemp<(xloc(k)+range);
%             dsy = find(ytemp>(yloc(k)-range) & ytemp<(yloc(k)+range));
            iloc = xtemp>(xloc(k)-range) & xtemp<(xloc(k)+range) & ytemp>(yloc(k)-range) & ytemp<(yloc(k)+range);
            zc              = ztemp(iloc);
            %         zc              = z(dsy,dsx,i);
            eta_pt_mean(k,i)   = nanmean(zc(:));
            eta_pt_std(k,i)    = nanstd(zc(:));
            clear ds* zc
        end
    end
    
    if extractTRANSECTy == 1
        for k = 1:length(xloctran)
            iloc = xtemp(1,:)>(xloctran(k)-(range/2)) & xtemp(1,:)<(xloctran(k)+(range/2));
%             [temp, ix] = nanmin(abs(xtemp(1,:)-xloctran(k)));
            eval(['x',num2str(xloctran(k)),'(:,i) = nanmean(xtemp(:,iloc),2);'])
            eval(['z',num2str(xloctran(k)),'(:,i) = nanmean(ztemp(:,iloc),2);'])
            eval(['y',num2str(xloctran(k)),'(:,i) = nanmean(ytemp(:,iloc),2);'])
        end
    end
    
   if extractTRANSECTx == 1
        for k = 1:length(yloctran)
            iloc = ytemp(:,1)>(yloctran(k)-(range/2)) & ytemp(:,1)<(yloctran(k)+(range/2));
%             [temp, ix] = nanmin(abs(xtemp(1,:)-xloctran(k)));
            eval(['x',num2str(abs(round(yloctran(k)))),'(i,:) = nanmean(xtemp(iloc,:),1);'])
            eval(['z',num2str(abs(round(yloctran(k)))),'(i,:) = nanmean(ztemp(iloc,:),1);'])
            eval(['y',num2str(abs(round(yloctran(k)))),'(i,:) = nanmean(ytemp(iloc,:),1);'])
        end
    end
    
    if extractREGION == 1
        
        j = ~isnan(ztemp);
        Xtemp = double(xtemp(j));
        Ytemp = double(ytemp(j));
        Ztemp = double(ztemp(j));
%         z(:,:,i) = griddata(Xtemp,Ytemp,Ztemp,x,y);
        
        [z(:,:,i),numPts(:,:,i)]=roundgridfun(Xtemp,Ytemp,Ztemp,x,y,@mean);
        
%         z(:,:,i) = interp2(repmat(xtemp,length(ytemp),1),repmat(ytemp',1,length(xtemp)),ztemp,x,y);


    end
    if plotFlag == 1 && i < 160
        plotFlag = 1;
    else
        plotFlag = 0;
    end
    
    for plot = 1 % plotting section 
        if plotFlag == 1
            % Plot of eta
%             yp   = ytemp;%repmat(ytemp',1,length(xtemp));
%             xp   = xtemp;%repmat(xtemp,length(ytemp),1);
%             zp   = ztemp;
            
            ax1 = axes('Position',[0.13 0.25 0.72 0.5]);
            pcolor(ytemp,xtemp,ztemp-tide-.03)
            shading flat
            colormap(cmocean('deep'));
%             colormap(flipud(winter));
            hc = colorbar('Position', [0.877 0.284 0.03 0.38]);
            caxis([-0.2 0.2])
            axis equal
            grid on
            box on
            shading interp
                    ylim([25.5 36])%[xROI(1) 35])%xROI(2)]) %ylim([20 39]) %xlim([-13 13])
                    xlim([yROI(1) yROI(2)]) %ylim([24 37])
            h1=gca;
            set(h1,'ydir','reverse')
            set(h1,'tickdir','out','xminortick','on','yminortick','on');
            set(h1,'ticklength',1.5*get(h1,'ticklength'));
            set(h1,'fontsize',15);
            %     set(h1,'xtick',[-10:5:10],'xticklabel',{'-10' '-5' '0' '5' '10'});
            text(ax1,13.7,25.9,'$\eta$ (m)','interpreter','latex','fontsize',18);
            xlabel('Alongshore (m)','interpreter','latex','fontsize',18);
            ylabel('Cross-Shore (m)','interpreter','latex','fontsize',18);
%             text(ax1,13.5,24.3,'$\eta~(\mathrm{m})$','interpreter','latex','fontsize',18);
%             ylabel(hc,'$\eta~(\mathrm{m})$','interpreter','latex','fontsize',18);
%             text(ax1,-13.5,21,[trialname(1:22),', ',num2str(imagestart/100,'%03.f'),num2str(imageno(i)-1,'%02.f')],'interpreter','latex','fontsize',18);
            text(ax1,7.3,24.7,['time: ',num2str(round(imageno(i)/8,3),'%.3f'),' sec'],'interpreter','latex','fontsize',18);
%             text(ax1,8,24.3,['frame ',nm,'/',num2str(numframes)],'interpreter','latex','fontsize',18);
            % ['time = ',num2str(round(imageno(i)/7,2)),' sec']
            
            Sname2 = [figfolder,trialname,'_tif_eta_1cmres_',num2str(imagestart/100,'%03.f'),num2str(imageno(i)-1,'%02.f'),subname];
            print(Sname2,'-dpng')
            
            clf
            
            if i < 10
                ax1 = axes('Position',[0.13 0.25 0.72 0.5]);
                pcolor(y,x,z(:,:,i)-tide)
                shading flat
                colormap(cmocean('balance'));
                hc = colorbar('Position', [0.87 0.25 0.03 0.5])
                caxis([-0.4 0.4])
                axis equal
                grid on
                box on
                shading interp
                ylim([xROI(1) xROI(2)]) %ylim([20 39]) %xlim([-13 13])
                xlim([yROI(1) yROI(2)]) %ylim([24 37])
                h1=gca;
                set(h1,'ydir','reverse')
                set(h1,'tickdir','out','xminortick','on','yminortick','on');
                set(h1,'ticklength',1.5*get(h1,'ticklength'));
                set(h1,'fontsize',15);
                %     set(h1,'xtick',[-10:5:10],'xticklabel',{'-10' '-5' '0' '5' '10'});
                xlabel('alongshore (m)','interpreter','latex','fontsize',18);
                ylabel('cross-shore (m)','interpreter','latex','fontsize',18);
                %             text(ax1,13.5,24.3,'$\eta~(\mathrm{m})$','interpreter','latex','fontsize',18);
                ylabel(hc,'$\eta~(\mathrm{m})$','interpreter','latex','fontsize',18);
                %             text(ax1,-13.5,21,[trialname(1:22),', ',num2str(imagestart/100,'%03.f'),num2str(imageno(i)-1,'%02.f')],'interpreter','latex','fontsize',18);
                text(ax1,7.3,24.3,['time: ',num2str(round(imageno(i)/8,3),'%.3f'),' sec'],'interpreter','latex','fontsize',18);
                %             text(ax1,8,24.3,['frame ',nm,'/',num2str(numframes)],'interpreter','latex','fontsize',18);
                % ['time = ',num2str(round(imageno(i)/7,2)),' sec']
                
                Sname2 = [figfolder,trialname,'_tif_eta_',num2str(imagestart/100,'%03.f'),num2str(imageno(i)-1,'%02.f'),'_griddata'];
                print(Sname2,'-dpng')
                
                clf
            end
            clear *temp Sname* nm file filename dr dsr xlen ax
            
        elseif plotFlag == 2 % flipped orientation
            % Plot of eta
            yp   = repmat(ytemp',1,length(xtemp));
            xp   = repmat(xtemp,length(ytemp),1);
            zp   = ztemp;
            
%             ax1 = axes('Position',[0.1 0.1 0.74 0.9]);
            pcolor(xp',yp',flipud((zp-tide)'))
            hold on
            shading interp
%             colormap(cmocean('balance'));
            colormap('thermal');
            colorbar('Position', [0.87 0.218 0.03 0.67])
            caxis([-0.3 0.3])
            axis equal
            grid on
            box on
            shading interp
            xlim([21 35]) %ylim([20 39]) %xlim([-13 13])
            ylim([-13 13]) %ylim([24 37])
            h1=gca;
%             set(h1,'xdir','reverse')
            set(h1,'tickdir','out','xminortick','on','yminortick','on');
            set(h1,'ticklength',1.5*get(h1,'ticklength'));
            set(h1,'fontsize',22);
            %     set(h1,'xtick',[-10:5:10],'xticklabel',{'-10' '-5' '0' '5' '10'});
            ylabel('Alongshore (m)','interpreter','latex','fontsize',18);
            xlabel('Cross-Shore (m)','interpreter','latex','fontsize',18);
            text(35.5,13,'$\eta~(\mathrm{m})$','interpreter','latex','fontsize',18);
%             text(ax1,-,21,[trialname(1:22),', ',num2str(imagestart/100,'%03.f'),num2str(imageno(i)-1,'%02.f')],'interpreter','latex','fontsize',18);
%             text(ax1,8,21,['frame ',nm,'/',num2str(numframes)],'interpreter','latex','fontsize',18);
            % ['time = ',num2str(round(imageno(i)/7,2)),' sec']
            fill([34 34 36 36],[-13 13 13 -13],[0.7 0.7 0.7],'LineWidth',2)
            
            Sname2 = [figfolder,trialname,'_tif_eta_',num2str(imagestart/100,'%03.f'),num2str(imageno(i)-1,'%02.f'),'_test'];
            print(Sname2,'-dpng')
            
            clear *temp Sname* nm file filename dr dsr xlen ax
            
            clf
        end
    end
end


if extractPT == 1
    sensorloc = 1;
    if sensorloc == 1
        psname1 = [savefolder,'dem_pt_sz_array_',num2str(((range*2)^2)*10000),'cm2.mat'];
    else
        if xloc(1)==xloc(2)
            psname1 = [savefolder,'dem_pt_x',num2str(xloc(1)),'m_',num2str(((range*2)^2)*10000),'cm2.mat'];
        elseif yloc(1)==yloc(2)
            psname1 = [savefolder,'dem_pt_y',num2str(yloc(1)),'m_',num2str(((range*2)^2)*10000),'m2.mat'];
        end
    end
eval(['save -v7.3 ',psname1,' eta_pt_mean',' eta_pt_std',' xloc',' yloc',' range',' inst']);
end

if extractTRANSECTy == 1
    for k = 1:length(xloctran)
        %         psname2 = [datafolder,'dem_transect_x',num2str(xloctran(k)),'m_',num2str(((range*2)^2)*10000),'cm2.mat'];
        psname2 = [savefolder,'dem_transect_x',num2str(xloctran(k)),'m_xavg',num2str(range*100),'cm.mat'];
%         xtransect = xloctran(k);
%         eval(['save -v7.3 ',psname2,' eta_mean_x',num2str(xloctran(k)),' eta_std_x',num2str(xloctran(k)),' xtransect',' yint',' dy']);
        eval(['save -v7.3 ',psname2,' x',num2str(xloctran(k)),' y',num2str(xloctran(k)),' z',num2str(xloctran(k))]);
    end
end

if extractTRANSECTx == 1
    for k = 1:length(yloctran)
        %         psname2 = [datafolder,'dem_transect_x',num2str(xloctran(k)),'m_',num2str(((range*2)^2)*10000),'cm2.mat'];
        if round(yloctran(k))>=0
            psname2 = [savefolder,'dem_transect_y',num2str(abs(round(yloctran(k)))),'m_yavg',num2str(range*100),'cm.mat'];
        elseif round(yloctran(k))<0
            psname2 = [savefolder,'dem_transect_yneg',num2str(abs(round(yloctran(k)))),'m_yavg',num2str(range*100),'cm.mat'];
        end
%         ytransect = yloctran(k);
%         eval(['save -v7.3 ',psname2,' eta_mean_x',num2str(xloctran(k)),' eta_std_x',num2str(xloctran(k)),' xtransect',' yint',' dy']);
        eval(['save -v7.3 ',psname2,' x',num2str(abs(round(yloctran(k)))),' y',num2str(abs(round(yloctran(k)))),' z',num2str(abs(round(yloctran(k))))]);
    end
end

if extractREGION == 1
    
        
        wlev = nanmean(z,3);
        Ptsmean = nanmean(numPts,3);
%                 wlev = z-tide;
%         eta = nanmean(wlev,3);
        Hs = 4*nanstd(z,[],3);
        psname2 = [savefolder,'dem_region_x',num2str(regx(1)),'to',num2str(regx(end)),'m_yneg',num2str(abs(regy(1))),'to',num2str(regy(end)),'m_res',num2str((dx)*100),'cm',subname,'.mat'];
        psname3 = [savefolder,'dem_region_x',num2str(regx(1)),'to',num2str(regx(end)),'m_yneg',num2str(abs(regy(1))),'to',num2str(regy(end)),'m_res',num2str((dx)*100),'cm_bulk_stat',subname,'.mat'];
%         xtransect = xloctran(k);
        eval(['save -v7.3 ',psname2,' x',' y',' z',' t',' numPts']);
%         eval(['save -v7.3 ',psname2,' x',' y',' z',' t',' wlev',' eta',' Hs']);
        eval(['save -v7.3 ',psname3,' x',' y',' t',' wlev',' Hs',' tide',' Ptsmean']);
end
close
end

