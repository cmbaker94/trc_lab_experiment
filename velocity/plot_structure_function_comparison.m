%% Read WAMFlow output
%  Makes WAMs for i2Rgus Camera Dec2020 Set Up

%% Housekeeping
close all
clear all
addpath('/Users/cmbaker9/Documents/MTOOLS')
addpath(genpath('/Users/cmbaker9/Documents/Research/swashzone_dynamics/codes/cameras'))
genpath('/Users/cmbaker9/Documents/Research/Lab_Experiments/codes/CIRN-Quantitative-Coastal-Imaging-Toolbox')
% addpath(genpath('/Users/cmbaker9/Documents/Research/Lab_Experiments/codes/CIRN-Quantitative-Coastal-Imaging-Toolbox/X_CoreFunctions'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))
cd E:\code\trc_lab_experiment\trial_processing

%% User Input

 Tinfo.Hs = 0.3;
    Tinfo.Tp = 2;
    Tinfo.tide = 1.07;
sprd = [0 20 30 40];
    
    idxdy=0.05;
ixlim=[25 35];
iylim=[-14 14];
    awin = 2;
    ovr  = 1;

    for isprd = 1:length(sprd)
        Tinfo.spread = sprd(isprd);
Tinfo = trial_files(Tinfo);

% Stereo Reconstructions
Tinfo.cam = TRC_camera_info(Tinfo.cam);
camera.time        = cam_time(Tinfo);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);


%% Load File Names + Times

subname = '';
geoname = ['x',num2str(ixlim(1)),'to',num2str(round(ixlim(2))),'_y',num2str(round(iylim(2))),'_res',num2str(idxdy*100),'cm',subname];
odir = [Tinfo.savefolder(1:74),'orthos\',geoname];

cnames = getfnames(odir,['c2_']);
L = cnames(:,4:end-5);

if length(L)>4802
    display('Temporarily picking subset of times')
    L = L(7201:12000,:);
end

t = str2num(L)/Tinfo.cam.Hz;
dt=mode(diff(t));
ts=(t-t(1));
Tbin= ((awin/2)+(ovr/2)):ovr: (ts(end)-awin/2-ovr/2); % Centered On

oname=['c2_',L(1,:),'_',L(end,:),'_win',num2str(awin),'_ovr',num2str(ovr)];
fname = fullfile(Tinfo.savefolder(1:74),'WAMS',geoname,oname,[oname,'R.',num2str(awin),'s.',num2str(ovr),'dt.optFlow.']);

%% Load data
load([fname,'inversecascade.mat']);
sf{isprd} = SF;
clear SF
    end

    eval(['!mkdir ',Tinfo.figfolder,'\',geoname,'\',oname])
%%

figure('units','inches','position',[1 1 8 5],'color','w')
scatter(sf{1}.cam.dy,mean(sf{1}.cam.Sdy,1),'b','LineWidth',1)
hold on
scatter(sf{2}.cam.dy,mean(sf{2}.cam.Sdy,1),'r','LineWidth',1)
scatter(sf{3}.cam.dy,mean(sf{3}.cam.Sdy,1),'g','LineWidth',1)
scatter(sf{4}.cam.dy,mean(sf{4}.cam.Sdy,1),'m','LineWidth',1)
% scatter(sf{5}.cam.dy,mean(sf{5}.cam.Sdy(25:51,:),1),'k','LineWidth',1)
% scatter(cam.dy,cam.Sdy(25,:))
% scatter(cam.dy,cam.Sdy(end,:))
scfit = [0.2:0.01:10];
% plot(scfit,0.00055*(scfit.^(2/3)),'Color',[0.7 0.7 0.7],'LineWidth',2,'LineStyle','-')
scale1 = 0.0013;
plot(scfit,scale1*(scfit.^(2/3)),'Color',[0.3 0.3 0.3],'LineWidth',2,'LineStyle','-')
scfit = [0:0.01:0.6];
% plot(scfit,0.00055*(scfit.^(2/3)),'Color',[0.7 0.7 0.7],'LineWidth',2,'LineStyle','-')
plot(scfit,0.008*(scfit.^(2)),'Color',[0.7 0.7 0.7],'LineWidth',2,'LineStyle','-')
scatter(sf{1}.ins.out.dy,sf{1}.ins.out.Sdy,'b','^','LineWidth',2)
scatter(sf{1}.ins.in.dy,sf{1}.ins.in.Sdy,'b','^','LineWidth',2)
scatter(sf{2}.ins.out.dy,sf{2}.ins.out.Sdy,'r','^','LineWidth',2)
scatter(sf{2}.ins.in.dy,sf{2}.ins.in.Sdy,'r','^','LineWidth',2)
scatter(sf{3}.ins.out.dy,sf{3}.ins.out.Sdy,'g','^','LineWidth',2)
scatter(sf{3}.ins.in.dy,sf{3}.ins.in.Sdy,'g','^','LineWidth',2)
scatter(sf{4}.ins.out.dy,sf{4}.ins.out.Sdy,'m','^','LineWidth',2)
scatter(sf{4}.ins.in.dy,sf{4}.ins.in.Sdy,'m','^','LineWidth',2)
set(gca,'xscale','log')
set(gca,'yscale','log')
% ylim([-0.3 0.3])
xlim([0.05 10])
box on
% grid on
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',1*get(h1,'ticklength'));
set(h1,'fontsize',20);
% title(['Velocity at p',num2str(ia),' $x$ = ',num2str(round(xy(1,ia),1)),' m, $y$ = ',num2str(round(xy(2,ia),2)),' m (WAMFlow*res and filtered in situ)'],'interpreter','latex','fontsize',15);
ylabel('Structure Function','interpreter','latex','fontsize',20);
xlabel('Lag (m)','interpreter','latex','fontsize',20);
h2 = legend('$\sigma_{\theta}=0^{\circ}$','$\sigma_{\theta}=20^{\circ}$','$\sigma_{\theta}=30^{\circ}$','$\sigma_{\theta}=40^{\circ}$','$S\sim(\Delta y)^{2/3}$','$S\sim(\Delta y)^{2}$','interpreter','latex','fontsize',15)
    set(h2, 'Location','southeast');%'Position', [0.88 0.777 0.10 0.14])
% title('$u$ (m/s)','interpreter','latex','fontsize',20);
sname =['WAMvsinsitu_res',geoname,'_awin',num2str(awin),'_ovr',(num2str(ovr)),'_structurefn'];
    print([Tinfo.figfolder,geoname,'\',oname,'\',sname],'-dpng')
    
    %%
    
    figure('units','inches','position',[1 1 8 5],'color','w')
scatter(sf{1}.cam.dy,mean(sf{1}.cam.Sdy(20:end,:),1),'b','LineWidth',1)
hold on
scatter(sf{2}.cam.dy,mean(sf{2}.cam.Sdy(20:end,:),1),'r','LineWidth',1)
scatter(sf{3}.cam.dy,mean(sf{3}.cam.Sdy(20:end,:),1),'g','LineWidth',1)
scatter(sf{4}.cam.dy,mean(sf{4}.cam.Sdy(20:end,:),1),'m','LineWidth',1)
% scatter(sf{5}.cam.dy,mean(sf{5}.cam.Sdy(25:51,:),1),'k','LineWidth',1)
% scatter(cam.dy,cam.Sdy(25,:))
% scatter(cam.dy,cam.Sdy(end,:))
scfit = [0.2:0.01:10];
% plot(scfit,0.00055*(scfit.^(2/3)),'Color',[0.7 0.7 0.7],'LineWidth',2,'LineStyle','-')
scale1 = 0.00141;
plot(scfit,scale1*(scfit.^(2/3)),'Color',[0.3 0.3 0.3],'LineWidth',2,'LineStyle','-')
scfit = [0:0.01:0.6];
% plot(scfit,0.00055*(scfit.^(2/3)),'Color',[0.7 0.7 0.7],'LineWidth',2,'LineStyle','-')
plot(scfit,0.008*(scfit.^(2)),'Color',[0.7 0.7 0.7],'LineWidth',2,'LineStyle','-')
scatter(sf{1}.ins.out.dy,sf{1}.ins.out.Sdy,'b','^','LineWidth',2)
scatter(sf{1}.ins.in.dy,sf{1}.ins.in.Sdy,'b','^','LineWidth',2)
scatter(sf{2}.ins.out.dy,sf{2}.ins.out.Sdy,'r','^','LineWidth',2)
scatter(sf{2}.ins.in.dy,sf{2}.ins.in.Sdy,'r','^','LineWidth',2)
scatter(sf{3}.ins.out.dy,sf{3}.ins.out.Sdy,'g','^','LineWidth',2)
scatter(sf{3}.ins.in.dy,sf{3}.ins.in.Sdy,'g','^','LineWidth',2)
scatter(sf{4}.ins.out.dy,sf{4}.ins.out.Sdy,'m','^','LineWidth',2)
scatter(sf{4}.ins.in.dy,sf{4}.ins.in.Sdy,'m','^','LineWidth',2)
set(gca,'xscale','log')
set(gca,'yscale','log')
% ylim([-0.3 0.3])
xlim([0.05 10])
box on
% grid on
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',1*get(h1,'ticklength'));
set(h1,'fontsize',20);
% title(['Velocity at p',num2str(ia),' $x$ = ',num2str(round(xy(1,ia),1)),' m, $y$ = ',num2str(round(xy(2,ia),2)),' m (WAMFlow*res and filtered in situ)'],'interpreter','latex','fontsize',15);
ylabel('Structure Function','interpreter','latex','fontsize',20);
xlabel('Lag (m)','interpreter','latex','fontsize',20);
h2 = legend('$\sigma_{\theta}=0^{\circ}$','$\sigma_{\theta}=20^{\circ}$','$\sigma_{\theta}=30^{\circ}$','$\sigma_{\theta}=40^{\circ}$','$S\sim(\Delta y)^{2/3}$','$S\sim(\Delta y)^{2}$','interpreter','latex','fontsize',15)
    set(h2, 'Location','southeast');%'Position', [0.88 0.777 0.10 0.14])
% title('$u$ (m/s)','interpreter','latex','fontsize',20);
sname =['WAMvsinsitu_res',geoname,'_awin',num2str(awin),'_ovr',(num2str(ovr)),'_structurefn'];
    print([Tinfo.figfolder,geoname,'\',oname,'\',sname],'-dpng')
   

    %%
    
figure('units','inches','position',[1 1 8 5],'color','w')
scatter(sf{1}.cam.dy,mean(sf{1}.cam.Sdy,1),10,'MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',2)
hold on
scatter(sf{1}.cam.dy,sf{1}.cam.Sdy(1,:),10,'MarkerEdgeColor',[66, 221, 245]/256,'LineWidth',2)
scatter(sf{1}.cam.dy,sf{1}.cam.Sdy(25,:),10,'MarkerEdgeColor',[66, 117, 219]/256,'LineWidth',2)
scatter(sf{1}.cam.dy,sf{1}.cam.Sdy(51,:),10,'MarkerEdgeColor',[28, 16, 201]/256,'LineWidth',2)
% scatter(cam.dy,cam.Sdy(25,:))
% scatter(cam.dy,cam.Sdy(end,:))
scatter(sf{1}.ins.out.dy,sf{1}.ins.out.Sdy,'r','LineWidth',2)
scatter(sf{1}.ins.in.dy,sf{1}.ins.in.Sdy,'m','LineWidth',2)
scfit = [0.3:0.01:10];
% plot(scfit,0.00055*(scfit.^(2/3)),'Color',[0.7 0.7 0.7],'LineWidth',2)
plot(scfit,0.004*(scfit.^(2/3)),'k','LineWidth',2,'LineStyle','-')
scfit = [0:0.01:0.8];
% plot(scfit,0.00055*(scfit.^(2/3)),'Color',[0.7 0.7 0.7],'LineWidth',2)
plot(scfit,0.01*(scfit.^(2)),'Color',[0.7 0.7 0.7],'LineWidth',2,'LineStyle','-')
set(gca,'xscale','log')
set(gca,'yscale','log')
ylim([10^(-5) 10^(-1.5)])
xlim([0.05 10])
box on
% grid on
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',1*get(h1,'ticklength'));
set(h1,'fontsize',20);
% title(['Velocity at p',num2str(ia),' $x$ = ',num2str(round(xy(1,ia),1)),' m, $y$ = ',num2str(round(xy(2,ia),2)),' m (WAMFlow*res and filtered in situ)'],'interpreter','latex','fontsize',15);
ylabel('Structure Function','interpreter','latex','fontsize',20);
xlabel('Lag (m)','interpreter','latex','fontsize',20);
h2 = legend('wamflow: sz avg',['wamflow: ','$x$=',num2str(x(ixrange(1))),' m'],...
    ['wamflow: ','$x$=',num2str(x(ixrange(25))),' m'],['wamflow: ','$x$=',num2str(x(ixrange(end))),' m'],...
    'in situ: outer sz','in situ: inner sz','$S\sim(\Delta y)^{2/3}$','$S\sim(\Delta y)^{2}$','interpreter','latex','fontsize',12)
    set(h2, 'Location','southeast');%'Position', [0.88 0.777 0.10 0.14])
% title('$u$ (m/s)','interpreter','latex','fontsize',20);
sname =['WAMvsinsitu_res',geoname,'_awin',num2str(awin),'_ovr',(num2str(ovr)),'_structurefn_difreg'];
    print([Tinfo.figfolder,geoname,'\',oname,'\',sname],'-dpng')