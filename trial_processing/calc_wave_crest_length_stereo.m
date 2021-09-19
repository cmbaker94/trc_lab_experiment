function calc_wave_crest_length_stereo(Tinfo)
% Quantify the short-crestedness of waves from stereo reconstructions and
% LiDAR imagery. This code will pick a cross-shore location to compute the
% alongshore spectra of sea-surface elevation at each data timesample and
% average for the data span.

% Set up paths and clear workspace
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\cameras'))

%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%

% % WC comparison
calc_xshore = [29 30 31];
dem = 1; % if 0 means region, if 1 means transect

%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION ABOVE %%%%%%%%%%%%%%%%%%%%%%%%%%

for cx = 1:length(calc_xshore)
    xshore = calc_xshore(cx);
    
    %% STEP 1: Create paths and load files
    
    % general path and names
    datapath    = 'E:/';
    
    % 
    Tinfo = trial_files(Tinfo);

    % Stereo Reconstructions
    Tinfo.cam = TRC_camera_info(Tinfo.cam);
    
    % Data and figure storage
    [Tinfo] = wc_comp_store(Tinfo);

    if dem == 0
        cam         = load([Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_res',num2str((Tinfo.cam.dx)*100),'cm_roundgridfun.mat']);
    elseif dem == 1
        % transect    = load(['/Users/cmbaker9/Documents/Research/Lab_Experiments/data/processed/DEM/TRM-08-30-2018-2119UTC_Scene1/frames_01000-01999/dem_transect_x31m_2.25cm2_redone.mat'])
        transect    = load([Tinfo.cam.datafolder,'dem_transect_x',num2str(xshore),'m_xavg5cm.mat']);
    end
    % Bathymetry
    bathy       = load([datapath,'data/processed/lidar/Riegl/TRC_bathymetry_estimate_line.mat']);
    
    %% STEP 2: Create figure folders
    
    % figure folder
    fssubfolder = datestr(date,'yy-mm-dd');
    figfolder   = [datapath,'figures/meas_comp/',Tinfo.cam.trialname,'/',Tinfo.cam.trimname,fssubfolder,'/'];
    
    % make figure folders
    eval(['!mkdir ',datapath,'figures/meas_comp/',Tinfo.cam.trialname]);
    eval(['!mkdir ',datapath,'figures/meas_comp/',Tinfo.cam.trialname,'/',Tinfo.cam.trimname]);
    eval(['!mkdir ',figfolder])
    
    %% STEP 3: Load Stereo Reconstruction Data
    

    clear *temp
    
    %% STEP 4: Choose cross-shore location
    
    % Let's choose the location based on the bathymetry and the location with
    % maximum points
    watlev = 1.07;
    bathy.depth = bathy.h-watlev;
    
    if dem == 0
        figure('units','inches','position',[1 1 7 14],'Color','w');
        subplot(3,1,[1 2])
        pcolorjw(cam.x,cam.y,nanmean(cam.numPts,3));%,100,'linestyle','none')%.*fliplr(beach2));
        %     axis equal
        shading flat
        xlim([Tinfo.cam.regx(1) Tinfo.cam.regx(end)])
        ylim([Tinfo.cam.regy(1) Tinfo.cam.regy(end)])
        colormap(cmocean('thermal'));
        h2 = colorbar('northoutside');
        ylabel(h2,'Avg. Pts per 0.05x0.05 cm$^2$','interpreter','latex','fontsize',20);
        %     caxis([0 100])
        h1=gca;
        set(h1,'tickdir','in','xminortick','on','yminortick','on');
        set(h1,'ticklength',2*get(h1,'ticklength'));
        set(h1,'ydir','normal');
        set(h1,'fontsize',20);
        xlabel('$x~\mathrm{(m)}$','interpreter','latex','fontsize',20);
        ylabel('$y~\mathrm{(m)}$','interpreter','latex','fontsize',20);
        box on
        
        subplot(3,1,3)
        plot(bathy.xp,bathy.depth,'k','LineWidth',2)
        hold on
        plot(bathy.xp,bathy.h*0,'k--','LineWidth',2);
        fill([bathy.xp; flipud(bathy.xp)],[bathy.depth; -2*ones(size(bathy.depth))],[0.7 0.7 0.7])
        xlim([Tinfo.cam.regx(1) Tinfo.cam.regx(end)])
        ylim([-1.1 0.1])
        xlabel('$x$ (m)','interpreter','latex','fontsize',18);
        ylabel('$h$ (m)','interpreter','latex','fontsize',18);
        grid on
        h1=gca;
        set(h1,'tickdir','in','xminortick','on','yminortick','on');
        set(h1,'ticklength',1*get(h1,'ticklength'));
        set(h1,'fontsize',20);
        
        Sname1 = [figfolder,Tinfo.cam.trialname ,'_',Tinfo.cam.imagerange,'_pick_alongshore_location'];
        print(Sname1,'-dpng')
    end
    
    %% STEP 5: Choose cross-shore location and extract sea-surface elevation
    
    % Ok, so by inpsection of the previous plot it looks like x=30 m may be a
    % good place to start... so here we go!
    if dem == 0
        xshore      = 27;
        dy          = cam.y(2,1)-cam.y(1,1);
        disp(['Here, this code is hard-coded to choose ',num2str(xshore),' m'])
        
        [temp,ix]   = nanmin(abs(xshore-cam.x(1,:)));
        % z           = nanmean(squeeze(cam.z(:,ix-2:ix+2,:)),2);
        % z           = squeeze(z(:,1,:));
    elseif dem == 1
        eval(['dy          = transect.y',num2str(xshore),'(2,1)-transect.y',num2str(xshore),'(1,1);'])
        eval(['z = squeeze(transect.z',num2str(xshore),'(:,:));'])
%         eval(['z = squeeze(transect.z',num2str(xshore),'(:,3,:));'])
        %     eval(['z = nanmean(transect.z',num2str(xshore),',2);'])
        %     z = squeeze(z(:,1,:));
        eval(['x = nanmean(transect.x',num2str(xshore),',1)'';'])
        %     z = transect.eta_mean_x31;
    end
    
    %% STEP 6: Compute spectra
    
    % WL =
    % OL =
    order = 3;
    framelen = 39;
    figure('units','inches','position',[1 1 12 5],'Color','w');
    
    dx = 0;
    
    count = 0
    for i = 6:size(z,2)-5
        count = 1+count;
        ztemp = squeeze(z(:,i-dx:i+dx));
        ztemp = double(nanmean(ztemp,2));
        ztemp(ztemp>1.35)=NaN;
        ztemp(ztemp<0.95)=NaN;
        eval(['ytemp = transect.y',num2str(xshore),'(:,',num2str(i),');'])
        ztemp(ytemp>13)=NaN;
        ztemp(ytemp<-13)=NaN;
        
        zfittemp = movmean(ztemp,11,'omitnan','EndPoints','fill');
%         zfittemp = ztemp;%sgolayfilt(ztemp,order,framelen);
        filt = [nanmean(zfittemp)-(2.5*nanstd(zfittemp)) nanmean(zfittemp)+(4*nanstd(zfittemp))];
        zfittemp(zfittemp<filt(1))=NaN;
        zfittemp(zfittemp>filt(2))=NaN;
        zfit(:,count)=zfittemp;
        % %     zfit(:,i) = smooth(ztemp,2);
        %     plot(ytemp,ztemp,'k','LineWidth',1.5)
        % %     plot(cam.y(:,1),ztemp,'k','LineWidth',1.5)
        %     hold on
        %     plot(ytemp,zfit(:,i),'r','LineWidth',1)
        % %     plot(cam.y(:,1),zfit(:,i),'r','LineWidth',1)
        %     grid on
        %     ylim([1 1.4])
        %     pause(0.2)
        %     clf
    end
    
    %%
    figure('units','inches','position',[1 1 7 7],'Color','w');
    count = 0;
    for i = 1:size(zfit,2)
        count = count+1;
        %     ztemp = squeeze(z(:,i));
        %     runmean = nanmean(squeeze(z(:,i-5:i+5)),2);
        ztemp = zfit(:,i);
        % remove nans and interpolate
        nanz            = isnan(ztemp);             % find location of NaNs
        nant            = [1:numel(ztemp)]';       % make strings for interp
        nanratio        = sum(nanz)/length(nant); % find ratio of NaNs
        ztemp(nanz)     = interp1(nant(~nanz),ztemp(~nanz),nant(nanz)); % interpolate between the NaNs
        zcutnan         = ztemp(~isnan(ztemp)); % if nans still at the beginning or end of list, just remove these points
        eta             = zcutnan-nanmean(zcutnan);
        eta             = detrend(eta,1);% DETRENDING!!!!!!!!
        if i<20
            %         plot(ztemp-nanmean(ztemp(:,i)),'k','LineWidth',1.5)
            %     plot(cam.y(:,1),ztemp,'k','LineWidth',1.5)
            %         hold on
            plot(zfit(:,i)-nanmean(zfit(:,i)),'r','LineWidth',1)
            hold on
            plot(eta,'b','LineWidth',1)
            %     plot(cam.y(:,1),zfit(:,i),'r','LineWidth',1)
            grid on
            %         ylim([1 1.4])
            pause(0.5)
            clf
        end
        
        
        if length(eta)>2350
            WL = length(eta);
            OL = 0;
%             [S(count,:),k(count,:),Sc(count,:,:)]   = pwelch(eta,WL,OL,[],1/dy,'ConfidenceLevel',0.95); % compute spectra
            [Stemp,ktemp,Sc(count,:,:)]   = pwelch(eta,WL,OL,[],1/dy,'ConfidenceLevel',0.95); % compute spectra
            [val,id] = min(abs(ktemp-1));
            filter = zeros(size(ktemp));
            filter(1:id-1) = 1;
            filter(id:id+2) = [0.7089 0.2363 0.0225];
            S(count,:) = Stemp.*filter;
            k(count,:) = ktemp;
        else
            S(count,:) = NaN(size(S(count-1,:)));
            k(count,:) = NaN(size(k(count-1,:)));
            Sc(count,:,:) = NaN(size(Sc(count-1,:,:)));
        end
        
        if i<20
            clf
            semilogy(k(count,:),S(count,:),'k','LineWidth',2)
            xlabel('$L^{-1}$ (m$^{-1}$)','interpreter','latex','fontsize',20);
            ylabel('$S_{\eta\eta}$ (m$^2$/Hz)','interpreter','latex','fontsize',20);
            h1=gca;
            set(h1, 'YScale', 'log')
            set(h1,'tickdir','in','xminortick','on','yminortick','on');
            set(h1,'ticklength',1*get(h1,'ticklength'));
            set(h1,'fontsize',15);
            %         xlim([0 1/(cam.y(end,1)-cam.y(1,1))])
            ylim([10^-6 10^-3])
            title(['x = ',num2str(xshore),'m , t = ',num2str(i/Tinfo.cam.Hz)],'interpreter','latex','fontsize',20);
            pause(0.2)
        end
    end
    
    
    %%
    Savg=nanmean(S,1);
    Scavg = nanmean(Sc,1);
    
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
    title(['x = ',num2str(xshore),'m , avg'],'interpreter','latex','fontsize',20);
    
    %%
    figure('units','inches','position',[1 1 7 10],'Color','w');
    zfitmean = nanmean(nanmean(zfit));
    figure
    pcolor((1:size(zfit(:,1:100),2)),ytemp,zfit(:,1:100)-zfitmean)
    caxis([-0.2 0.2])
    shading flat
    colormap(cmocean('balance'));
    % xlim([0 400])
    
    %%
    figure('units','inches','position',[1 1 7 10],'Color','w');
    zfitmean = nanmean(nanmean(zfit));
    figure
    pcolor((1:size(zfit(:,1:1000),2)),ytemp,zfit(:,1:1000)-zfitmean)
    caxis([-0.2 0.2])
    shading flat
    colormap(cmocean('balance'));
    %%
    
    y = ytemp;
    x = nanmean(x);
    
    
    subname = '';%'onewindow_11movmean_detrended';
    psname = [Tinfo.cam.datafolder,'dem_transect_wavecrest_length_x',num2str(xshore),'_',subname,'.mat'];
    eval(['save -v7.3 ',psname,' S',' k',' Sc',' Savg',' Scavg',' y',' x',' zfit']);

close all
clear zfit
end
end

