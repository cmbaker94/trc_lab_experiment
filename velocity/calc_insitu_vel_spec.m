% Figure of tank geomoetry and in situ gages. 

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

sprd = [0; 10; 20; 30; 40];

%% STEP 2: Define figure folders

% figure folder
datapath = 'E:\';
fssubfolder = datestr(date,'yy-mm-dd');
figfolder   = [datapath,'figures\meas_comp\velocity\',fssubfolder,'\'];
eval(['!mkdir ',figfolder])


%% Loop through two trials

for si = 1:length(sprd)
    
    Tinfo.spread = sprd(si);
    Tinfo.Hs = 0.25;
    Tinfo.Tp = 2;
    Tinfo.tide = 1.07;
    
    Tinfo = trial_files(Tinfo);
    
    % general path and names
    Tinfo.cam = TRC_camera_info(Tinfo.cam);
    
    % Data and figure storage
    [Tinfo] = wc_comp_store(Tinfo);
    
    % STEP 1: Load Data
    
    % load insitu data
    for ia = 1:2
        if ia == 1
            F1          = matfile([datapath,'data/processed/insitu/',Tinfo.sz.tdate,'/',Tinfo.sz.tdate,'-insitu.mat']);
            arrayname = 'sz';
        elseif ia == 2
            F1          = matfile([datapath,'data/processed/insitu/',Tinfo.is.tdate,'/',Tinfo.is.tdate,'-insitu.mat']);
            arrayname = 'is';
        end
    
    % press       = F1.press;
    % wg          = F1.wg;
    Tcam.is_insitu    = F1.Tinfo;
    vel       = F1.vel;
    ptemp   	= F1.press;
    Hz          = ptemp.Hz;
    clear ptemp
    
    time_range = [4 44]*60;
    t_inst = 0:1/Hz:(length(vel.u1)-1)*(1/Hz);
    % find overlapping range for insitu
    [temp,istart] = nanmin(abs(time_range(1)-t_inst));
    [temp,iend] = nanmin(abs(time_range(2)-t_inst));
    inst.time = t_inst(istart:iend)';
    istart = istart - Tinfo.offsets(1);
    iend = iend - Tinfo.offsets(1);
    
    velnames  = fieldnames(vel.xyzd);
    
    mvavgt= 8;% sec
    
    spec2avg = 10;
    % compute bulk statistics
    for i = 1:length(velnames)
        % load pressure data
        instnum = str2double(velnames{i}(4:end));
        %     display(instnum)
        eval(['u = vel.u',num2str(instnum),'(istart:iend);'])
        eval(['v = vel.v',num2str(instnum),'(istart:iend);'])
        WL = 256;%round(length(u)/(spec2avg*Hz));
        OL = WL/2;%round(length(u)/(spec2avg*2*Hz));
        vmag = sqrt(u.^2+v.^2);
        [Svel(:,instnum),f(:,instnum),Svelc]   = pwelch(detrend(vmag),WL*Hz,OL*Hz,[],Hz,'ConfidenceLevel',0.95);
        [Su(:,instnum),f(:,instnum),Suc]   = pwelch(detrend(u),WL*Hz,OL*Hz,[],Hz,'ConfidenceLevel',0.95);
        [Sv(:,instnum),f(:,instnum),Svc(:,:,instnum)]   = pwelch(detrend(v),WL*Hz,OL*Hz,[],Hz,'ConfidenceLevel',0.95);
        %     figure; loglog(f,Svel(:,instnum)); hold on; loglog(f,Su(:,instnum)); loglog(f,Sv(:,instnum)); loglog(f,Su+Sv(:,instnum))
        instname{instnum}             = velnames{i};
        eval(['xyz(instnum,:) = vel.xyzd.',velnames{i},'(1,1:3);'])
        
        fmax = 0.02;%15;
        [val,ifreq] = min(abs(f(:,1)-fmax));
        power(:,instnum) = trapz(f(3:ifreq,instnum),Su(3:ifreq,instnum)+Sv(3:ifreq,instnum));
    end
    % clear vel
    if ia == 1
        irange = 1:8;
    elseif ia == 2
        irange = 4:12;
    end
    Su_avg = nanmean(Su(:,irange),2);
    Sv_avg = nanmean(Sv(:,irange),2);
    Svel_avg = nanmean(Su(:,irange)+Sv(:,irange),2);
    
    eval(['T',num2str(Tinfo.spread),'.',arrayname,'.Svel_avg = Svel_avg;'])
    eval(['T',num2str(Tinfo.spread),'.',arrayname,'.Sv_avg = Sv_avg;'])
    eval(['T',num2str(Tinfo.spread),'.',arrayname,'.Su_avg = Su_avg;'])
    eval(['T',num2str(Tinfo.spread),'.',arrayname,'.Svel = Su+Sv;'])
    eval(['T',num2str(Tinfo.spread),'.',arrayname,'.Sv = Sv;'])
    eval(['T',num2str(Tinfo.spread),'.',arrayname,'.Su = Su;'])
    eval(['T',num2str(Tinfo.spread),'.',arrayname,'.Svc = Svc;'])
    eval(['T',num2str(Tinfo.spread),'.',arrayname,'.f = f(:,1);'])
    eval(['T',num2str(Tinfo.spread),'.',arrayname,'.power = power;'])
    eval(['T',num2str(Tinfo.spread),'.',arrayname,'xyz = xyz;'])
    clear vel u v Su* Sv* power
    end
    clear Tinfo
end

f = f(:,1);
%%

% cm = colormap((cmocean('phase')));
% cm = colormap('jet');
% cmap = [1 round([1:4]*(length(cm(:,1))/4))];
% alt = [5 15 25 30 0];
% cmap = cmap+alt;
% add cm(cmap(1),:) to plots
c1 = [3, 71, 50]/256;
c2 = [236, 0, 252]/256;

figure('units','inches','position',[1 1 18 6],'color','w')
ax1 = axes('Position',[0.17 0.15 0.35 0.79]);
loglog(f,T0.sz.Svel_avg,'Color','k','LineWidth',1.5) 
hold on
loglog(f,T10.sz.Svel_avg,'Color','r','LineWidth',1.5) 
loglog(f,T20.sz.Svel_avg,'Color',[0.2 0.8 0.2],'LineWidth',1.5)
loglog(f,T30.sz.Svel_avg,'Color','b','LineWidth',1.5)
loglog(f,T40.sz.Svel_avg,'Color',c2,'LineWidth',1.5)
hfill = fill([f(3) f(3) f(ifreq) f(ifreq)],[10^-2  1 1 10^-2],[239 242 157]/256,'LineStyle','none','EdgeColor','none')
alpha(hfill,0.4)
conf = [0.84 1.21];
plot([10^(-0.15) 10^(-0.15)],[10^(-0.18)*conf(1) 10^(-0.18)*conf(2)],'k','LineWidth',1)
plot([10^(-0.13) 10^(-0.17)],[10^(-0.18)*conf(1)+0.00002 10^(-0.18)*conf(1)+0.00002],'k','LineWidth',1)
plot([10^(-0.13) 10^(-0.17)],[10^(-0.18)*conf(2)-0.00002 10^(-0.18)*conf(2)-0.00002],'k','LineWidth',1)
scatter(10^(-0.15),10^(-0.18),20,'k','fill')
plot([10^(-2) 10^(-2)],[10^-2 10^0],'Color',[.8 .8 .8],'LineWidth',0.1)
plot([10^(-1) 10^(-1)],[10^-2 10^0],'Color',[.8 .8 .8],'LineWidth',0.1)
plot([10^(-3) 10^0],[10^(-1) 10^(-1)],'Color',[.8 .8 .8],'LineWidth',0.1)
loglog(f,T0.sz.Svel_avg,'Color','k','LineWidth',1.5) 
loglog(f,T10.sz.Svel_avg,'Color','r','LineWidth',1.5) 
loglog(f,T20.sz.Svel_avg,'Color',[0.2 0.8 0.2],'LineWidth',1.5)
loglog(f,T30.sz.Svel_avg,'Color','b','LineWidth',1.5)
loglog(f,T40.sz.Svel_avg,'Color',c2,'LineWidth',1.5)

% grid on
h1=gca;
set(h1, 'XScale', 'log')
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',20);
xlim([0.00305 10^0])
ylim([10^-2 10^0])
xlabel('$f$ (Hz)','interpreter','latex','fontsize',22)
ylabel('$\langle S_{uu}+S_{vv} \rangle$ ((m$^2$/s$^2$)/Hz)','interpreter','latex','fontsize',22)
% set(h1,'yticklabel',[],'xticklabel',[]);
text(ax1,10^(-2.46), 10^(-0.15),'(a) $x = 28.4$ m','interpreter','latex','fontsize',24);
% title(ax3,'$\sigma_{\theta}=40^{\circ}$','interpreter','latex','fontsize',28);
% errorbar(T40.specon.wg.f(67),10^(-1.3),10^-1.85,'Color','k','LineWidth',2.5)%wg.See(3)-wg.Seec(3,1))
h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 10^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$','VLF','interpreter','latex','fontsize',22)
set(h2,'orientation','vertical','Position',[0.005 0.60 0.1 0.3])


szin = [T0.sz.power(:,9:12); T10.sz.power(:,9:12); T20.sz.power(:,9:12); T30.sz.power(:,9:12); T40.sz.power(:,9:12)];
    szout = [T0.sz.power(:,1:8); T10.sz.power(:,1:8); T20.sz.power(:,1:8); T30.sz.power(:,1:8); T40.sz.power(:,1:8)];
    isin = [T0.is.power(:,4:12); T10.is.power(:,4:12); T20.is.power(:,4:12); T30.is.power(:,4:12); T40.is.power(:,4:12)];
    isout = [T0.is.power(:,1:3); T10.is.power(:,1:3); T20.is.power(:,1:3); T30.is.power(:,1:3); T40.is.power(:,1:3)];
    
    szin(szin==0) = NaN;
    szout(szout==0) = NaN;
    isin(isin==0) = NaN;
    isout(isout==0) = NaN;
    
    save(['E:\data\processed\trial_comparison\total_variance.mat'],'szin','szout','isin','isout')
        disp('saved Mat File')
    
    ax2 = axes('Position',[0.575 0.15 0.30 0.79]);
    plot(100-sprd,nanmean(isout,2),'k','LineWidth',3)
    hold on
    plot(100-sprd,nanmean(isin,2),'r','LineWidth',3)
    plot(100-sprd,nanmean(szout,2),'m','LineWidth',3)
    plot(100-sprd,nanmean(szin,2),'b','LineWidth',3)
    
    plot(sprd,nanmean(isout,2),'k','LineWidth',1)
    plot(sprd,nanmean(isin,2),'r','LineWidth',1)
    plot(sprd,nanmean(szout,2),'m','LineWidth',1)
    plot(sprd,nanmean(szin,2),'b','LineWidth',1)
    errorbar(sprd,nanmean(isout,2),nanstd(isout,[],2),'Color','k','LineWidth',1.5,'LineStyle','None')
    errorbar(sprd,nanmean(isin,2),nanstd(isin,[],2),'Color','k','LineWidth',1.5,'LineStyle','None')
    errorbar(sprd,nanmean(szout,2),nanstd(szout,[],2),'Color','k','LineWidth',1.5,'LineWidth',1.5,'LineStyle','None')
    errorbar(sprd,nanmean(szin,2),nanstd(szin,[],2),'Color','k','LineWidth',1.5,'LineStyle','None')
    scatter(sprd,nanmean(isout,2),170,'k','fill','MarkerEdgeColor','k','LineWidth',1)
    scatter(sprd,nanmean(isin,2),170,'r','fill','MarkerEdgeColor','k','LineWidth',1)
    scatter(sprd,nanmean(szout,2),170,'fill','MarkerFaceColor','m','MarkerEdgeColor','k','LineWidth',1)
    scatter(sprd,nanmean(szin,2),170,'b','fill','MarkerEdgeColor','k','LineWidth',1)
%     grid on
    h1=gca;
    % set(h1, 'XScale', 'log')
    set(h1,'tickdir','in','xminortick','on','yminortick','on');
    set(h1,'ticklength',1.5*get(h1,'ticklength'));
    set(h1,'fontsize',20);
    box on
    xlim([-2 42])
    ylim([0 4.2]*10^-3)
    xlabel('$\sigma_{\theta}$ ($^{\circ}$)','interpreter','latex','fontsize',22)
    ylabel('$\sigma_{vel}$ (m$^2$/s$^2$)','interpreter','latex','fontsize',22)
    text(ax2,-0.7, 3.9*10^-3,'(b)','interpreter','latex','fontsize',24)
    % set(h1,'yticklabel',[],'xticklabel',[]);
    % text(ax3,0.05, 10^(-1.3),'(c)','interpreter','latex','fontsize',24);
    % title(ax3,'$\sigma_{\theta}=40^{\circ}$','interpreter','latex','fontsize',28);
    % errorbar(T40.specon.wg.f(67),10^(-1.3),10^-1.85,'Color','k','LineWidth',2.5)%wg.See(3)-wg.Seec(3,1))
    h2 = legend('$x = 24.5$ m','$x = 26.6$ m','$x = 28.4$ m','$x = 30.7$ m','interpreter','latex','fontsize',22)
    set(h2,'orientation','vertical','Position',[0.89 0.68 0.1 0.26])
   sname = ['vel_power'];
    print([figfolder,sname],'-dpng')









%%
inst =4;
figure('units','inches','position',[1 1 12 8],'color','w')
plot(f,T0.sz.Svel(:,inst),'k') 
hold on
plot(f,T10.sz.Svel(:,inst),'r') 
plot(f,T20.sz.Svel(:,inst),'g')
plot(f,T30.sz.Svel(:,inst),'b')
plot(f,T40.sz.Svel(:,inst),'m')
conf = [0.84 1.21];
plot([10^(-0.45) 10^(-0.45)],[10^(-0.15)*conf(1) 10^(-0.15)*conf(2)],'k','LineWidth',1)
plot([10^(-0.43) 10^(-0.47)],[10^(-0.15)*conf(1)+0.00002 10^(-0.15)*conf(1)+0.00002],'k','LineWidth',1)
plot([10^(-0.43) 10^(-0.47)],[10^(-0.15)*conf(2)-0.00002 10^(-0.15)*conf(2)-0.00002],'k','LineWidth',1)
h1=gca;
set(h1, 'XScale', 'log')
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',20);
xlim([min(f) 10^0])
% ylim([10^-5 10^-1])
% set(h1,'yticklabel',[],'xticklabel',[]);
% text(ax3,0.05, 10^(-1.3),'(c)','interpreter','latex','fontsize',24);
% title(ax3,'$\sigma_{\theta}=40^{\circ}$','interpreter','latex','fontsize',28);
% errorbar(T40.specon.wg.f(67),10^(-1.3),10^-1.85,'Color','k','LineWidth',2.5)%wg.See(3)-wg.Seec(3,1))
h2 = legend('0','10','20','30','40','interpreter','latex','fontsize',22)
set(h2,'orientation','vertical')

%%

isxyz = T0.isxyz;
szxyz = T0.szxyz;

x = -3;
if x==0
    inst = [11; 6; 8; 1];
elseif x == 3
    inst = [12; 8; 9];
elseif x == -3
    inst = [10; 4; 7];
elseif x == -6
    inst = [9; 2; 6];
end

    szin = [T0.sz.power(:,inst(1)); T10.sz.power(:,inst(1)); T20.sz.power(:,inst(1)); T30.sz.power(:,inst(1)); T40.sz.power(:,inst(1))];
    szout = [T0.sz.power(:,inst(2)); T10.sz.power(:,inst(2)); T20.sz.power(:,inst(2)); T30.sz.power(:,inst(2)); T40.sz.power(:,inst(2))];
    isin = [T0.is.power(:,inst(3)); T10.is.power(:,inst(3)); T20.is.power(:,inst(3)); T30.is.power(:,inst(3)); T40.is.power(:,inst(3))];
    if length(inst)==4
        isout = [T0.is.power(:,inst(4)); T10.is.power(:,inst(4)); T20.is.power(:,inst(4)); T30.is.power(:,inst(4)); T40.is.power(:,inst(4))];
    else
        isout = NaN(size(szout));
    end
        figure('units','inches','position',[1 1 6 6],'color','w')
%     scatter(0,T0.power(:,inst),'k')
%     hold on
%     scatter(10,T10.power(:,inst),'k')
%     scatter(20,T20.power(:,inst),'k')
%     scatter(30,T30.power(:,inst),'k')
%     scatter(40,T40.power(:,inst),'k')
    scatter(sprd,isout,'k')
    hold on
    scatter(sprd,isin,'r')
    scatter(sprd,szout,'g')
    scatter(sprd,szin,'b')

    h1=gca;
    % set(h1, 'XScale', 'log')
    set(h1,'tickdir','in','xminortick','on','yminortick','on');
    set(h1,'ticklength',1.5*get(h1,'ticklength'));
    set(h1,'fontsize',20);
    box on
    xlim([-5 45])
    % ylim([10^-5 10^-1])
    % set(h1,'yticklabel',[],'xticklabel',[]);
    % text(ax3,0.05, 10^(-1.3),'(c)','interpreter','latex','fontsize',24);
    % title(ax3,'$\sigma_{\theta}=40^{\circ}$','interpreter','latex','fontsize',28);
    % errorbar(T40.specon.wg.f(67),10^(-1.3),10^-1.85,'Color','k','LineWidth',2.5)%wg.See(3)-wg.Seec(3,1))
    h2 = legend('isout','isin','szout','szin','interpreter','latex','fontsize',22)
    set(h2,'orientation','vertical','Location','southeast')
    
%%
% sprd = [0; 10; 20; 30; 40];
% T10.sz.power = NaN(size(T20.sz.power));
% T10.is.power = NaN(size(T20.is.power));

% x = -3;
    

%%

figure('units','inches','position',[1 1 12 8],'color','w')
plot(f,T0.sz.Svel_avg,'k') 
hold on
plot(f,T10.sz.Svel_avg,'r') 
plot(f,T20.sz.Svel_avg,'g')
plot(f,T30.sz.Svel_avg,'b')
plot(f,T40.sz.Svel_avg,'m')
h1=gca;
set(h1, 'XScale', 'log')
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',20);
xlim([min(f) 10^0])
% ylim([10^-5 10^-1])
% set(h1,'yticklabel',[],'xticklabel',[]);
% text(ax3,0.05, 10^(-1.3),'(c)','interpreter','latex','fontsize',24);
% title(ax3,'$\sigma_{\theta}=40^{\circ}$','interpreter','latex','fontsize',28);
% errorbar(T40.specon.wg.f(67),10^(-1.3),10^-1.85,'Color','k','LineWidth',2.5)%wg.See(3)-wg.Seec(3,1))
h2 = legend('0','10','20','30','40','interpreter','latex','fontsize',22)
set(h2,'orientation','vertical')

%%



%%
figure; loglog(f,T30.sz.Sv(:,4)); hold on; loglog(f,T30.sz.Svc(:,1,4));loglog(f,T30.sz.Svc(:,2,4))
figure; plot(f,T30.sz.Sv(:,4)); hold on; plot(f,T30.sz.Svc(:,1,4));plot(f,T30.sz.Svc(:,2,4)); h1 = gca; set(h1, 'XScale', 'log')

%%


figure; loglog(f,T0.Sv,'k'); hold on; loglog(f,T20.Sv,'r'); loglog(f,T40.Sv,'b')
figure; loglog(f,T0.Su,'k'); hold on; loglog(f,T20.Su,'r'); loglog(f,T40.Su,'b')