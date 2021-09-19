function extract_Hs_stats(Tinfo)
% Plot Hs and sea-surface evolutino spectra for remote sensing and in situ
% gages.

% Plot Bulk Statics
% Read in the lower resolution region and compute bulk statistics

% Set up paths and clear workspace
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\codes\insitu'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iarray = 1:2
    
    Tinfo = trial_files(Tinfo);
    
    % general path and names
    Tinfo.cam = TRC_camera_info(Tinfo.cam);
    
    % Data and figure storage
    [Tinfo] = wc_comp_store(Tinfo);
    
    %% STEP 1: Load data
    
    datapath    = 'E:\'; % general path and names
    % Tinfo.comp = ['Hs',num2str(Tinfo.Hs*100),'_Tp',num2str(Tinfo.Tp),'_tide',num2str(Tinfo.tide*100),'_spread',num2str(Tinfo.spread)];
    % procpath = [datapath,'data\processed\conditions\',Tinfo.comp,'\',Tinfo.tstart,'\'];
    
    % bathymetry
    bathymetry = load('E:\data\processed\lidar\Riegl\TRC_bathymetry_estimate_line.mat');
    h = Tinfo.tide-bathymetry.h;
    x = bathymetry.xp;
    
    % Time series 1
    % trange = Sname;
    load([Tinfo.savefolder,'\insitu_Hs_spectra.mat'], 'sz');
    load([Tinfo.savefolder,'\insitu_Hs_spectra.mat'], 'is');
    if Tinfo.filt == 1
        load([Tinfo.savefolder,'\velodyne_Hs_spectra_filt.mat'], 'lidar');
        load([Tinfo.savefolder,'\stereo_Hs_spectra_filt.mat'], 'cam');
    else
        load([Tinfo.savefolder,'\velodyne_Hs_spectra.mat'], 'lidar');
        load([Tinfo.savefolder,'\stereo_Hs_spectra.mat'], 'cam');
    end
    
    %% STEP 2: Define figure folders
    
    % % figure folder
    % fssubfolder = datestr(date,'yy-mm-dd');
    % figfolder   = [datapath,'figures\meas_comp\',savefolder,'\frames_07200-11999\'];
    % eval(['!mkdir ',figfolder])
    % figfolder   = [figfolder,fssubfolder,'\'];
    % eval(['!mkdir ',figfolder])
    
    %% STEP 3: remove sensors\areas with poor data
    
    % remove side walls & offshore
    % LiDAR
    lidar.Hs(lidar.y<-12.8,:)=NaN;
    lidar.Hs(lidar.y>12.8,:)=NaN;
    lidar.Hs(:,lidar.x<24) = NaN;
    lidar.See(lidar.y<-12.8,:,:)=NaN;
    lidar.See(lidar.y>12.8,:)= NaN;
    lidar.See(:,lidar.x<24,:) = NaN;
    
    % Stereo
    cam.Hs(cam.y(:,1)<-12.9,:)=NaN;
    cam.Hs(cam.y(:,1)>13,:)=NaN;
    cam.Hs(:,cam.x(1,:)<28.35) = NaN;
    cam.See(cam.y(:,1)<-12.9,:,:)=NaN;
    cam.See(cam.y(:,1)>13,:,:)=NaN;
    cam.See(:,cam.x(1,:)<28.35,:) = NaN;
    
    % avg and stdev
    lidar.Hs_yavg = nanmean(lidar.Hs,1);
    cam.Hs_yavg = nanmean(cam.Hs,1);
    lidar.See_yavg = squeeze(nanmean(lidar.See,1));
    cam.See_yavg = squeeze(nanmean(cam.See,1));
    
    lidar.Hs_ystd = nanstd(lidar.Hs,1);
    cam.Hs_ystd = nanstd(cam.Hs,1);
    lidar.See_ystd = squeeze(nanstd(lidar.See,1));
    cam.See_ystd = squeeze(nanstd(cam.See,1));
    
    for i = 1:length(sz.Hs)
        if sz.Hs(i)<0.05
            sz.Hs(i) = NaN;
            sz.See(i,:) = NaN;
        end
    end
    for i = 1:length(is.Hs)
        if is.Hs(i)<0.05
            is.Hs(i) = NaN;
            is.See(i,:) = NaN;
        end
    end
    
    %% STEP 4: compare spectra
    
    %     gages = {'wg4';'press2'};
    if iarray == 1
        insitu = sz;
        disp('sz')
        instrum_names = fieldnames(sz.data);
        instrum_names = instrum_names(2:end);
        if Tinfo.filt == 1
            savename = [Tinfo.savefolder,'\Hs_Tp_sz_filt'];
        else
            savename = [Tinfo.savefolder,'\Hs_Tp_sz'];
        end
    elseif iarray == 2
        insitu = is;
        disp('is')
        instrum_names = fieldnames(sz.data);
        instrum_names = instrum_names(2:end);
         if Tinfo.filt == 1
            savename = [Tinfo.savefolder,'\Hs_Tp_is_filt'];
        else
            savename = [Tinfo.savefolder,'\Hs_Tp_is'];
        end
    end
    
    
    grab_tseries = 0;
    for i = 1:length(instrum_names)
        gage = instrum_names(i);
        spectemp = extract_lab_spectra_1gageloc(insitu,cam,lidar,gage,grab_tseries);
        eval(['spec.',gage{1},' = spectemp;'])
    end
    eval(['save -v7.3 ',savename,' spec'])
    
end
end
