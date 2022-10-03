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
projpath = '/Users/cmbaker9/Documents/Research/swashzone_dynamics/';

%% User Input

 Tinfo.Hs = 0.3;
    Tinfo.Tp = 2;
    Tinfo.tide = 1.07;
    Tinfo.spread = 40;
    scale = 0; % if scale output is 0 then it's already been scaled for resolution, if 1 then needs to be still (aka all previously processed videos)
    
    idxdy=0.05;
ixlim=[25 35];
iylim=[-14 14];
    awin = 4;
    ovr  = 2;

Tinfo = trial_files(Tinfo);

% Stereo Reconstructions
Tinfo.cam = TRC_camera_info(Tinfo.cam);
% camera.time        = cam_time(Tinfo);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

%% Load File Names + Times

subname = '';
geoname = ['x',num2str(ixlim(1)),'to',num2str(round(ixlim(2))),'_y',num2str(round(iylim(2))),'_res',num2str(idxdy*100),'cm',subname];
odir = [Tinfo.savefolder(1:74),'orthos\',geoname];

cnames = getfnames(odir,['c2_']);
L = cnames(:,4:end-5);

display('Temporarily picking subset')
L = L(7201:12000,:);

t = str2num(L)/Tinfo.cam.Hz;
dt=mode(diff(t));
ts=(t-t(1));
Tbin= (awin/2):ovr: (ts(end)-awin/2-ovr); % Centered On
% Tbin= ((awin/2)+(ovr/2)):ovr: (ts(end)-awin/2-ovr/2); % Centered On

oname=['c2_',L(1,:),'_',L(end,:),'_win',num2str(awin),'_ovr',num2str(ovr)];
fname = fullfile(Tinfo.savefolder(1:74),'WAMS',geoname,oname,[oname,'R.',num2str(awin),'s.',num2str(ovr),'dt.optFlow.']);

%% Read data
tic
for i = 1:length(Tbin) 
    A = readtable([fname,'u',num2str(i-1),'.txt']);
    vel_u(:,:,i) = A{:,:};
%     utemp = movmedian(A{:,:},3,1);
%     utemp = movmedian(utemp,3,2);
%     vel_u(:,:,i) = utemp;
    A = readtable([fname,'v',num2str(i-1),'.txt']);
    vel_v(:,:,i) = A{:,:};
%     vtemp = movmedian(A{:,:},3,1);
%     vtemp = movmedian(vtemp,3,2);
%     vel_v(:,:,i) = vtemp;
%     A = readtable([fname,'stdFrames',num2str(i-1),'.txt']);
%     stdFrames(:,:,i) = A{:,:};
    display(i)
end
toc

if scale == 1
vel_u = vel_u.*(idxdy/ovr);
vel_v = vel_v.*(idxdy/ovr);
end


starttemp       = datenum(Tinfo.cam.tstart(1:end-3),'mm-dd-yyyy-HHMM')+datenum(0,0,0,0,0,str2double(L(1,:))/Tinfo.cam.Hz);
endtemp         = datenum(Tinfo.cam.tstart(1:end-3),'mm-dd-yyyy-HHMM')+datenum(0,0,0,0,0,(str2double(L(end,:)))/Tinfo.cam.Hz);
camera.time        = starttemp:datenum(0,0,0,0,0,1/Tinfo.cam.Hz):endtemp;
cam_t = camera.time(Tbin*Tinfo.cam.Hz); 
x = ixlim(1):idxdy:ixlim(2);
y = iylim(1):idxdy:iylim(2);
time = cam_t;

tic
save([fname,'compiled.mat'],'vel_u','vel_v','x','y','time')
disp('saved Mat File')
clear time
toc

%% Load WAM

load(fullfile(Tinfo.savefolder(1:74),'WAMS',geoname,oname,[oname 'R.mat']))

for i = 1:size(Iwam,3)-1
    Itest = movmean(Iwam(:,:,i),1,5);
    Itest = movmean(Itest,2,5);
    Ipixm(:,:,i) = Itest;
    if i>1
        Idiff(:,:,i-1) = Ipixm(:,:,i)-Ipixm(:,:,i-1); 
    end
end

%% Load insitu data

% [~,~,vel] = load_insitu(Tinfo);
[~,~,vel] = load_insitu_fullrun(Tinfo,camera.time);

%% Define wamflow properties and extract location with in situ gages

tsel = cam_t(end)+datenum(0,0,0,0,0,awin/2);
[~,itime] = nanmin(abs(tsel-vel.time));
ins_t  = vel.time;

% vel_uf=vel_u;
% vel_uf(abs(vel_uf)<0.01)=NaN;
% vel_vf=vel_v;
% vel_vf(abs(vel_vf)<0.01)=NaN;

for i = 1:12
    ip = find(contains(vel.loc,['vel',num2str(i,'%02.f')]));
    if ~isempty(ip)
        xyfind = vel.xyz(1:2,ip);
        [~,ix] = nanmin(abs(x-xyfind(1)));
        [~,iy] = nanmin(abs(y-xyfind(2)));
        xy(:,i) = xyfind;
        camIdiff(:,i) = squeeze(Idiff(iy,ix,:));
        cam_u(:,i) = squeeze(median(median(vel_u(iy-2:iy+2,ix-12:ix+2,:),1,'omitnan'),2,'omitnan'));
        cam_v(:,i) = squeeze(median(median(vel_v(iy-2:iy+2,ix-2:ix+2,:),1,'omitnan'),2,'omitnan'));
        ins_u(:,i) = vel.u(:,ip);
        ins_v(:,i) = vel.v(:,ip);
%         ins_u(:,i) = vel.u(1:itime,ip);
%         ins_v(:,i) = vel.v(1:itime,ip);
    else
        cam_u(:,i) = NaN(1,size(vel_u,3));
        cam_v(:,i) = NaN(1,size(vel_v,3));
        ins_u(:,i) = NaN(1,size(vel.u,1));
        ins_v(:,i) = NaN(1,size(vel.v,1));
%         ins_u(:,i) = NaN(1,itime);
%         ins_v(:,i) = NaN(1,itime);
    end
end

%% Time-average in situ data
% 
% for j = 1:length(Tbin)
%     irange = round(1+((100*(Tbin(j)-awin/2):(100*(Tbin(j)+awin/2)))));
% %     irange = 1+((100*(Tbin(j)):(100*(Tbin(j)))));
%     insmm_u(j,:) = nanmean(ins_u(irange,:),1);
%     insmm_v(j,:) = nanmean(ins_v(irange,:),1);
% end

% for j = 1:12
%     insmm2
% end

% ins_vel =  sqrt((ins_u).^2+(ins_v).^2);
% insmm_vel = sqrt((insmm_u).^2+(insmm_v).^2);
% cam_vel = sqrt((cam_u).^2+(cam_v).^2);

%%


for i = 1:12
	[r2(1,i),bias(1,i),rmse(1,i),pfit(1,i),slope(1,i),r(1,i)] = r2_bias_stats(insmm_u(:,i),cam_u(:,i));
    [r2(2,i),bias(2,i),rmse(2,i),pfit(2,i),slope(2,i),r(2,i)] = r2_bias_stats(insmm_v(:,i),cam_v(:,i));
    [r2(3,i),bias(3,i),rmse(3,i),pfit(3,i),slope(3,i),r(3,i)] = r2_bias_stats(insmm_vel(:,i),cam_vel(:,i));
end

%%
tpass = 32;
mm = 10/ovr;
irange = round(100*Tbin(1):100*ovr:100*Tbin(end));
inspl_t = ins_t(irange);
insplm_t = inspl_t(mm/2:mm:end);
camplm_t = cam_t(mm/2:mm:end);
% it = find(cam_t==ins_t);
for i = 1:12
    insrl_u(:,i) = smoothdata(ins_u(:,i),'rloess',12*20*ovr);
    insrl_v(:,i) = smoothdata(ins_v(:,i),'rloess',12*20*ovr);
    
    
    tempu = pl64ta(insrl_u(:,i),100*tpass);
    tempv = pl64ta(insrl_v(:,i),100*tpass);
    inspl_u(:,i) = tempu(irange);
    inspl_v(:,i) = tempv(irange);
    tempu = movmedian(tempu(irange),mm);
    tempv = movmedian(tempv(irange),mm);
    insplm_u(:,i) = tempu(mm/2:mm:end);
    insplm_v(:,i) = tempv(mm/2:mm:end);
    insplm2_u(:,i) = inspl_u(mm/2:mm:end,i);
    insplm2_v(:,i) = inspl_v(mm/2:mm:end,i);
    
    
    camrl_u(:,i) = smoothdata(cam_u(:,i),'rloess',12*ovr);
    camrl_v(:,i) = smoothdata(cam_v(:,i),'rloess',12*ovr);
    camrlpl_u(:,i) = pl64ta(camrl_u(:,i),tpass/ovr);
    camrlpl_v(:,i) = pl64ta(camrl_v(:,i),tpass/ovr);
    tempu = movmedian(camrlpl_u(:,i),mm);
    tempv = movmedian(camrlpl_v(:,i),mm);
    camrlplm_u(:,i) = tempu(mm/2:mm:end);
    camrlplm_v(:,i) = tempv(mm/2:mm:end);
    camrlplm2_u(:,i) = camrlpl_u(mm/2:mm:end,i);
    camrlplm2_v(:,i) = camrlpl_v(mm/2:mm:end,i);
    
    campl_u(:,i) = pl64ta(cam_u(:,i),camf);
    campl_v(:,i) = pl64ta(cam_v(:,i),camf);
    tempu = movmedian(campl_u(:,i),mm);
    tempv = movmedian(campl_v(:,i),mm);
    camplm_u(:,i) = tempu(mm/2:mm:end);
    camplm_v(:,i) = tempv(mm/2:mm:end);
    camplm2_u(:,i) = campl_u(mm/2:mm:end,i);
    camplm2_v(:,i) = campl_v(mm/2:mm:end,i);
end

campl_vel = sqrt((campl_u).^2+(campl_v).^2);
inspl_vel = sqrt((inspl_u).^2+(inspl_v).^2);

camplm_vel = sqrt((camplm_u).^2+(camplm_v).^2);
insplm_vel = sqrt((insplm_u).^2+(insplm_v).^2);

figure; plot(inspl_t,inspl_v(:,10)); hold on; plot(cam_t,camrlpl_v(:,10))

%%
save([fname,'insitu.mat'],'ins_u','ins_v','insmm_u','insmm_v','cam_u','cam_v','cam_t','ins_t','xy')
disp('saved Mat File')

%% Compute spectra
WL = 64;
OL = WL/2;

for i = 1:12
    if sum(~isnan(ins_u(:,i))) > 0
        Hz = 100;
        [ins_Svel(:,i),ins_f(:,i)]              = pwelch(detrend(ins_vel(:,i)),WL*Hz,OL*Hz,[],Hz,'ConfidenceLevel',0.95);
        [ins_Su(:,i),ins_f(:,i)]                = pwelch(detrend(ins_u(:,i)),WL*Hz,OL*Hz,[],Hz,'ConfidenceLevel',0.95);
        [ins_Sv(:,i),ins_f(:,i),ins_Svc(:,:,i)] = pwelch(detrend(ins_v(:,i)),WL*Hz,OL*Hz,[],Hz,'ConfidenceLevel',0.95);
        Hz = 1/ovr;
        [cam_Svel(:,i),cam_f(:,i)]              = pwelch(detrend(fillmissing(cam_vel(:,i),'linear')),WL*Hz,OL*Hz,[],Hz,'ConfidenceLevel',0.95);
        [cam_Su(:,i),cam_f(:,i)]                = pwelch(detrend(fillmissing(cam_u(:,i),'linear')),WL*Hz,OL*Hz,[],Hz,'ConfidenceLevel',0.95);
        [cam_Sv(:,i),cam_f(:,i),cam_Svc(:,:,i)] = pwelch(detrend(fillmissing(cam_v(:,i),'linear')),WL*Hz,OL*Hz,[],Hz,'ConfidenceLevel',0.95);
    end
end


%% Try lowess filter
% I think this'll only work for xy domain

% test2 = fit(cam_t',test,'lowess');
% test2 = fit(cam_t',test,'smoothingspline');
% test3 = pl64ta(test,5);
% plot(test2,cam_t,test)
% hold on
% scatter(cam_t,test3)

%% Make PLOTS!

eval(['!mkdir ',Tinfo.figfolder,geoname,'\',oname])

figure('units','inches','position',[1 1 12 8],'color','w')
for ia = 1:12
subplot(2,1,1)
plot(ins_t,ins_u(:,ia),'k')
hold on
plot(cam_t,cam_u(:,ia),'r')
plot(cam_t,zeros(size(cam_t)),'LineStyle','-.','LineWidth',0.5,'Color',[0.5 0.5 0.5])
datetick('x','MM:SS')
xlim([ins_t(1) ins_t(end)])
ylim([-2 2])
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
title(['Velocity at p',num2str(ia),' $x$ = ',num2str(round(xy(1,ia),1)),' m, $y$ = ',num2str(round(xy(2,ia),2)),' m'],'interpreter','latex','fontsize',15);
xlabel('$t$ (MM:SS)','interpreter','latex','fontsize',15);
ylabel('$u$ (m/s)','interpreter','latex','fontsize',15);

subplot(2,1,2)
plot(ins_t,ins_v(:,ia),'k')
hold on
plot(cam_t,cam_v(:,ia),'r')
plot(cam_t,zeros(size(cam_t)),'LineStyle','-.','LineWidth',0.5,'Color',[0.5 0.5 0.5])
datetick('x','MM:SS')
xlim([ins_t(1) ins_t(end)])
ylim([-2 2])
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
% text(ax1,25,450,'(a) Image','interpreter','latex','fontsize',15,'Color','w');
xlabel('$t$ (MM:SS)','interpreter','latex','fontsize',15);
ylabel('$v$ (m/s)','interpreter','latex','fontsize',15);

sname =['WAMvsinsitu_p',num2str(ia),'_res',geoname,'_awin',num2str(awin),'_ovr',(num2str(ovr)),'_raw'];
    print([Tinfo.figfolder,geoname,'\',oname,'\',sname],'-dpng')
    clf
end
close

%%
fnum = 5;
   figure('units','inches','position',[1 1 12 8],'color','w')
for ia = 1:12
subplot(3,1,1)
% plot(ins_t,movmean(ins_u(:,ia),6*100),'k')
plot(cam_t,insmm_u(:,ia),'k','LineWidth',1.5)
hold on
plot(cam_t,cam_u(:,ia),'r')
plot(cam_t,medfilt1(cam_u(:,ia),fnum),'b','LineWidth',1.5)
plot(cam_t,zeros(size(cam_t)),'LineStyle','-.','LineWidth',0.5,'Color',[0.5 0.5 0.5])
datetick('x','MM:SS')
xlim([ins_t(1) ins_t(end)])
ylim([-0.5 0.5])
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
title(['Velocity at p',num2str(ia),' $x$ = ',num2str(round(xy(1,ia),1)),' m, $y$ = ',num2str(round(xy(2,ia),2)),' m (WAMFlow*res and filtered in situ)'],'interpreter','latex','fontsize',15);
xlabel('$t$ (MM:SS)','interpreter','latex','fontsize',15);
ylabel('$u$ (m/s)','interpreter','latex','fontsize',15);

subplot(3,1,2)
% plot(ins_t,movmean(ins_v(:,ia),6*100),'k')
plot(cam_t,insmm_v(:,ia),'k','LineWidth',1.5)
hold on
plot(cam_t,cam_v(:,ia),'r')
plot(cam_t,medfilt1(cam_v(:,ia),fnum),'b','LineWidth',1.5)
plot(cam_t,zeros(size(cam_t)),'LineStyle','-.','LineWidth',0.5,'Color',[0.5 0.5 0.5])
datetick('x','MM:SS')
xlim([ins_t(1) ins_t(end)])
ylim([-0.5 0.5])
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
% text(ax1,25,450,'(a) Image','interpreter','latex','fontsize',15,'Color','w');
xlabel('$t$ (MM:SS)','interpreter','latex','fontsize',15);
ylabel('$v$ (m/s)','interpreter','latex','fontsize',15);

subplot(3,1,3)
% plot(ins_t,movmean(ins_v(:,ia),6*100),'k')
plot(cam_t,insmm_vel(:,ia),'k','LineWidth',1.5)
hold on
plot(cam_t,cam_vel(:,ia),'r')
plot(cam_t,calc_mag(medfilt1(cam_u(:,ia),fnum),medfilt1(cam_v(:,ia),fnum)),'b','LineWidth',1.5)
plot(cam_t,zeros(size(cam_t)),'LineStyle','-.','LineWidth',0.5,'Color',[0.5 0.5 0.5])
datetick('x','MM:SS')
xlim([ins_t(1) ins_t(end)])
ylim([0 0.6])
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
% text(ax1,25,450,'(a) Image','interpreter','latex','fontsize',15,'Color','w');
xlabel('$t$ (MM:SS)','interpreter','latex','fontsize',15);
ylabel('$|\vec{u}|$ (m/s)','interpreter','latex','fontsize',15);

sname =['WAMvsinsitu_p',num2str(ia),'_res',geoname,'_awin',num2str(awin),'_ovr',(num2str(ovr)),'_resampled'];
    print([Tinfo.figfolder,geoname,'\',oname,'\',sname],'-dpng')
clf
end
close

%%

%%
fnum = 5;
   figure('units','inches','position',[1 1 12 8],'color','w')
for ia = 1:12
subplot(3,1,1)
% plot(ins_t,movmean(ins_u(:,ia),6*100),'k')
plot(inspl_t,inspl_u(:,ia)-nanmean(inspl_u(:,ia)),'k','LineWidth',1.5)
hold on
plot(cam_t,campl_u(:,ia)-nanmean(campl_u(:,ia)),'r')
plot(cam_t,zeros(size(cam_t)),'LineStyle','-.','LineWidth',0.5,'Color',[0.5 0.5 0.5])
datetick('x','MM:SS')
xlim([ins_t(1) ins_t(end)])
ylim([-0.5 0.5])
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
title(['Velocity at p',num2str(ia),' $x$ = ',num2str(round(xy(1,ia),1)),' m, $y$ = ',num2str(round(xy(2,ia),2)),' m (WAMFlow*res and filtered in situ)'],'interpreter','latex','fontsize',15);
xlabel('$t$ (MM:SS)','interpreter','latex','fontsize',15);
ylabel('$u$ (m/s)','interpreter','latex','fontsize',15);

subplot(3,1,2)
% plot(ins_t,movmean(ins_v(:,ia),6*100),'k')
plot(inspl_t,inspl_v(:,ia)-nanmean(inspl_v(:,ia)),'k','LineWidth',1.5)
hold on
plot(cam_t,campl_v(:,ia)-nanmean(campl_v(:,ia)),'r')
plot(cam_t,zeros(size(cam_t)),'LineStyle','-.','LineWidth',0.5,'Color',[0.5 0.5 0.5])
datetick('x','MM:SS')
xlim([ins_t(1) ins_t(end)])
ylim([-0.5 0.5])
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
% text(ax1,25,450,'(a) Image','interpreter','latex','fontsize',15,'Color','w');
xlabel('$t$ (MM:SS)','interpreter','latex','fontsize',15);
ylabel('$v$ (m/s)','interpreter','latex','fontsize',15);

subplot(3,1,3)
% plot(ins_t,movmean(ins_v(:,ia),6*100),'k')
plot(inspl_t,inspl_vel(:,ia)-nanmean(inspl_vel(:,ia)),'k','LineWidth',1.5)
hold on
plot(cam_t,campl_vel(:,ia)-nanmean(campl_vel(:,ia)),'r')
plot(cam_t,zeros(size(cam_t)),'LineStyle','-.','LineWidth',0.5,'Color',[0.5 0.5 0.5])
datetick('x','MM:SS')
xlim([ins_t(1) ins_t(end)])
ylim([0 0.6])
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
% text(ax1,25,450,'(a) Image','interpreter','latex','fontsize',15,'Color','w');
xlabel('$t$ (MM:SS)','interpreter','latex','fontsize',15);
ylabel('$|\vec{u}|$ (m/s)','interpreter','latex','fontsize',15);

sname =['WAMvsinsitu_p',num2str(ia),'_res',geoname,'_awin',num2str(awin),'_ovr',(num2str(ovr)),'_pl64'];
    print([Tinfo.figfolder,geoname,'\',oname,'\',sname],'-dpng')
clf
end
close

%%

fnum = 5;
   figure('units','inches','position',[1 1 12 8],'color','w')
for ia = 1:12
subplot(3,1,1)
% plot(ins_t,movmean(ins_u(:,ia),6*100),10,'k')
scatter(cam_t,insmm_u(:,ia),10,'k','LineWidth',1.5)
hold on
scatter(cam_t,cam_u(:,ia),10,'r')
% scatter(cam_t,medfilt1(cam_u(:,ia),fnum),'b','LineWidth',1.5)
plot(cam_t,zeros(size(cam_t)),'LineStyle','-.','LineWidth',0.5,'Color',[0.5 0.5 0.5])
datetick('x','MM:SS')
xlim([ins_t(1) ins_t(end)])
ylim([-0.5 0.5])
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
title(['Velocity at p',num2str(ia),' $x$ = ',num2str(round(xy(1,ia),1)),' m, $y$ = ',num2str(round(xy(2,ia),2)),' m (WAMFlow*res and filtered in situ)'],'interpreter','latex','fontsize',15);
xlabel('$t$ (MM:SS)','interpreter','latex','fontsize',15);
ylabel('$u$ (m/s)','interpreter','latex','fontsize',15);

subplot(3,1,2)
% plot(ins_t,movmean(ins_v(:,ia),6*100),'k')
scatter(cam_t,insmm_v(:,ia),10,'k','LineWidth',1.5)
hold on
scatter(cam_t,cam_v(:,ia),10,'r')
% scatter(cam_t,medfilt1(cam_v(:,ia),fnum),10,'b','LineWidth',1.5)
plot(cam_t,zeros(size(cam_t)),'LineStyle','-.','LineWidth',0.5,'Color',[0.5 0.5 0.5])
datetick('x','MM:SS')
xlim([ins_t(1) ins_t(end)])
ylim([-0.5 0.5])
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
% text(ax1,25,450,'(a) Image','interpreter','latex','fontsize',15,'Color','w');
xlabel('$t$ (MM:SS)','interpreter','latex','fontsize',15);
ylabel('$v$ (m/s)','interpreter','latex','fontsize',15);

subplot(3,1,3)
% plot(ins_t,movmean(ins_v(:,ia),6*100),'k')
scatter(cam_t,insmm_vel(:,ia),10,'k','LineWidth',1.5)
hold on
scatter(cam_t,cam_vel(:,ia),10,'r')
% scatter(cam_t,calc_mag(medfilt1(cam_u(:,ia),fnum),medfilt1(cam_v(:,ia),fnum)),10,'b','LineWidth',1.5)
plot(cam_t,zeros(size(cam_t)),'LineStyle','-.','LineWidth',0.5,'Color',[0.5 0.5 0.5])
datetick('x','MM:SS')
xlim([ins_t(1) ins_t(end)])
ylim([0 0.6])
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
% text(ax1,25,450,'(a) Image','interpreter','latex','fontsize',15,'Color','w');
xlabel('$t$ (MM:SS)','interpreter','latex','fontsize',15);
ylabel('$|\vec{u}|$ (m/s)','interpreter','latex','fontsize',15);

sname =['WAMvsinsitu_p',num2str(ia),'_res',geoname,'_awin',num2str(awin),'_ovr',(num2str(ovr)),'_resampled_scatter'];
    print([Tinfo.figfolder,geoname,'\',oname,'\',sname],'-dpng')
clf
end
close

%%
figure('units','inches','position',[1 1 12 5],'color','w')
for ia = 1:12
subplot(1,3,1)
scatter(insmm_u(:,ia)-nanmean(insmm_u(:,ia)),cam_u(:,ia)-nanmean(cam_u(:,ia)),'k')
hold on
scatter(insmm_u(:,ia)-nanmean(insmm_u(:,ia)),camrl_u(:,ia)-nanmean(camrl_u(:,ia)),'r')
plot([-1 1],[-1 1],'LineStyle','-.','LineWidth',1.5,'Color',[0.5 0.5 0.5])
axis equal
ylim([-0.3 0.3])
xlim([-0.3 0.3])
box on
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
% title(['Velocity at p',num2str(ia),' $x$ = ',num2str(round(xy(1,ia),1)),' m, $y$ = ',num2str(round(xy(2,ia),2)),' m (WAMFlow*res and filtered in situ)'],'interpreter','latex','fontsize',15);
xlabel('in situ','interpreter','latex','fontsize',15);
ylabel('wamflow','interpreter','latex','fontsize',15);
title('$u$ (m/s)','interpreter','latex','fontsize',15);

subplot(1,3,2)
scatter(insmm_v(:,ia)-nanmean(insmm_v(:,ia)),cam_v(:,ia)-nanmean(cam_v(:,ia)),'k')
hold on
scatter(insmm_v(:,ia)-nanmean(insmm_v(:,ia)),camrl_v(:,ia)-nanmean(camrl_v(:,ia)),'r')
plot([-1 1],[-1 1],'LineStyle','-.','LineWidth',1.5,'Color',[0.5 0.5 0.5])
axis equal
ylim([-0.3 0.3])
xlim([-0.3 0.3])
box on
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
% text(ax1,25,450,'(a) Image','interpreter','latex','fontsize',15,'Color','w');
xlabel('in situ','interpreter','latex','fontsize',15);
% ylabel('wamflow','interpreter','latex','fontsize',15);
title('$v$ (m/s)','interpreter','latex','fontsize',15);

subplot(1,3,3)
scatter(insmm_vel(:,ia),cam_vel(:,ia),'k')
hold on
plot([-1 1],[-1 1],'LineStyle','-.','LineWidth',1.5,'Color',[0.5 0.5 0.5])
axis equal
ylim([0 0.6])
xlim([0 0.6])
box on
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
xlabel('in situ','interpreter','latex','fontsize',15);
% ylabel('wamflow','interpreter','latex','fontsize',15);
title('$|\vec{u}|$ (m/s)','interpreter','latex','fontsize',15);


sname =['WAMvsinsitu_p',num2str(ia),'_res',geoname,'_awin',num2str(awin),'_ovr',(num2str(ovr)),'_scatter'];
    print([Tinfo.figfolder,geoname,'\',oname,'\',sname],'-dpng')
    clf
end
close

%% trying new pl64 method
figure('units','inches','position',[1 1 12 5],'color','w')
for ia = 1:12
subplot(1,3,1)
scatter(inspl_u(:,ia)-nanmean(inspl_u(:,ia)),campl_u(:,ia)-nanmean(campl_u(:,ia)),'k')
hold on
plot([-1 1],[-1 1],'LineStyle','-.','LineWidth',1.5,'Color',[0.5 0.5 0.5])
axis equal
ylim([-0.3 0.3])
xlim([-0.3 0.3])
box on
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
% title(['Velocity at p',num2str(ia),' $x$ = ',num2str(round(xy(1,ia),1)),' m, $y$ = ',num2str(round(xy(2,ia),2)),' m (WAMFlow*res and filtered in situ)'],'interpreter','latex','fontsize',15);
xlabel('in situ','interpreter','latex','fontsize',15);
ylabel('wamflow','interpreter','latex','fontsize',15);
title('$u$ (m/s)','interpreter','latex','fontsize',15);

subplot(1,3,2)
scatter(inspl_v(:,ia)-nanmean(inspl_v(:,ia)),campl_v(:,ia)-nanmean(campl_v(:,ia)),'k')
hold on
plot([-1 1],[-1 1],'LineStyle','-.','LineWidth',1.5,'Color',[0.5 0.5 0.5])
axis equal
ylim([-0.3 0.3])
xlim([-0.3 0.3])
box on
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
% text(ax1,25,450,'(a) Image','interpreter','latex','fontsize',15,'Color','w');
xlabel('in situ','interpreter','latex','fontsize',15);
% ylabel('wamflow','interpreter','latex','fontsize',15);
title('$v$ (m/s)','interpreter','latex','fontsize',15);

subplot(1,3,3)
scatter(inspl_vel(:,ia)-nanmean(inspl_vel(:,ia)),cam_vel(:,ia)-nanmean(cam_vel(:,ia)),'k')
hold on
plot([-1 1],[-1 1],'LineStyle','-.','LineWidth',1.5,'Color',[0.5 0.5 0.5])
axis equal
ylim([0 0.6])
xlim([0 0.6])
box on
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
xlabel('in situ','interpreter','latex','fontsize',15);
% ylabel('wamflow','interpreter','latex','fontsize',15);
title('$|\vec{u}|$ (m/s)','interpreter','latex','fontsize',15);

sname =['WAMvsinsitu_p',num2str(ia),'_res',geoname,'_awin',num2str(awin),'_ovr',(num2str(ovr)),'_scatter'];
    print([Tinfo.figfolder,geoname,'\',oname,'\',sname],'-dpng')
    clf
end
close

%% filtered but subsampled

figure('units','inches','position',[1 1 12 5],'color','w')
for ia = 1:12
subplot(1,3,1)
scatter(insplm_u(:,ia)-nanmean(insplm_u(:,ia)),camplm_u(:,ia)-nanmean(camplm_u(:,ia)),'r')
hold on
scatter(insplm_u(:,ia)-nanmean(insplm_u(:,ia)),camrlplm_u(:,ia)-nanmean(camrlplm_u(:,ia)),'k')
hold on
plot([-1 1],[-1 1],'LineStyle','-.','LineWidth',1.5,'Color',[0.5 0.5 0.5])
axis equal
ylim([-0.3 0.3])
xlim([-0.3 0.3])
box on
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
% title(['Velocity at p',num2str(ia),' $x$ = ',num2str(round(xy(1,ia),1)),' m, $y$ = ',num2str(round(xy(2,ia),2)),' m (WAMFlow*res and filtered in situ)'],'interpreter','latex','fontsize',15);
xlabel('in situ','interpreter','latex','fontsize',15);
ylabel('wamflow','interpreter','latex','fontsize',15);
title('$u$ (m/s)','interpreter','latex','fontsize',15);

subplot(1,3,2)
scatter(insplm_v(:,ia)-nanmean(insplm_v(:,ia)),camplm_v(:,ia)-nanmean(camplm_v(:,ia)),'r')
hold on
scatter(insplm_v(:,ia)-nanmean(insplm_v(:,ia)),camrlplm_v(:,ia)-nanmean(camrlplm_v(:,ia)),'k')
hold on
plot([-1 1],[-1 1],'LineStyle','-.','LineWidth',1.5,'Color',[0.5 0.5 0.5])
axis equal
ylim([-0.3 0.3])
xlim([-0.3 0.3])
box on
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
% text(ax1,25,450,'(a) Image','interpreter','latex','fontsize',15,'Color','w');
xlabel('in situ','interpreter','latex','fontsize',15);
% ylabel('wamflow','interpreter','latex','fontsize',15);
title('$v$ (m/s)','interpreter','latex','fontsize',15);

subplot(1,3,3)
scatter(insplm_vel(:,ia)-nanmean(insplm_vel(:,ia)),camplm_vel(:,ia)-nanmean(camplm_vel(:,ia)),'k')
hold on
plot([-1 1],[-1 1],'LineStyle','-.','LineWidth',1.5,'Color',[0.5 0.5 0.5])
axis equal
ylim([0 0.6])
xlim([0 0.6])
box on
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
xlabel('in situ','interpreter','latex','fontsize',15);
% ylabel('wamflow','interpreter','latex','fontsize',15);
title('$|\vec{u}|$ (m/s)','interpreter','latex','fontsize',15);

% pause(2)
sname =['WAMvsinsitu_p',num2str(ia),'_res',geoname,'_awin',num2str(awin),'_ovr',(num2str(ovr)),'_scatter'];
    print([Tinfo.figfolder,geoname,'\',oname,'\',sname],'-dpng')
    clf
end
close

%%
figure('units','inches','position',[1 1 12 8],'color','w')
for ia = 1:12
plot(ins_f,ins_Svel(:,ia),'k') 
hold on
plot(cam_f,cam_Svel(:,ia),'k','LineStyle','-.') 
plot(ins_f,ins_Su(:,ia),'r') 
plot(cam_f,cam_Su(:,ia),'r','LineStyle','-.') 
plot(ins_f,ins_Sv(:,ia),'b') 
plot(cam_f,cam_Sv(:,ia),'b','LineStyle','-.') 
% conf = [0.84 1.21];
% plot([10^(-0.45) 10^(-0.45)],[10^(-0.15)*conf(1) 10^(-0.15)*conf(2)],'k','LineWidth',1)
% plot([10^(-0.43) 10^(-0.47)],[10^(-0.15)*conf(1)+0.00002 10^(-0.15)*conf(1)+0.00002],'k','LineWidth',1)
% plot([10^(-0.43) 10^(-0.47)],[10^(-0.15)*conf(2)-0.00002 10^(-0.15)*conf(2)-0.00002],'k','LineWidth',1)
h1=gca;
set(h1, 'XScale', 'log', 'YScale', 'log')
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',20);
% xlim([min(f) 10^0])
 ylim([10^-3 10^-0.5])
% set(h1,'yticklabel',[],'xticklabel',[]);
% text(ax3,0.05, 10^(-1.3),'(c)','interpreter','latex','fontsize',24);
% title(ax3,'$\sigma_{\theta}=40^{\circ}$','interpreter','latex','fontsize',28);
% errorbar(T40.specon.wg.f(67),10^(-1.3),10^-1.85,'Color','k','LineWidth',2.5)%wg.See(3)-wg.Seec(3,1))
% h2 = legend('0','10','20','30','40','interpreter','latex','fontsize',22)
% set(h2,'orientation','vertical') 
pause(1)
    
sname =['WAMvsinsitu_p',num2str(ia),'_res',geoname,'_awin',num2str(awin),'_ovr',(num2str(ovr)),'_spectra'];
    print([Tinfo.figfolder,geoname,'\',oname,'\',sname],'-dpng')
    clf
end
close

%%

figure('units','inches','position',[1 1 8 8],'color','w')
for ia = 1:9
scatter(nanmean(insmm_vel(:,ia)),nanmean(cam_vel(:,ia)),'k')
hold on
scatter(nanmean(insmm_u(:,ia)),nanmean(cam_u(:,ia)),'r')
scatter(nanmean(insmm_v(:,ia)),nanmean(cam_v(:,ia)),'b')
end
for ia = 10:12
scatter(nanmean(insmm_vel(:,ia)),nanmean(cam_vel(:,ia)),'k','fill')
hold on
scatter(nanmean(insmm_u(:,ia)),nanmean(cam_u(:,ia)),'r','fill')
scatter(nanmean(insmm_v(:,ia)),nanmean(cam_v(:,ia)),'b','fill')
end
plot([-1 1],[-1 1],'LineStyle','-.','LineWidth',1.5,'Color',[0.5 0.5 0.5])
axis equal
ylim([-0.3 0.3])
xlim([-0.3 0.3])
box on
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','off','yminortick','off');
set(h1,'ticklength',0.5*get(h1,'ticklength'));
set(h1,'fontsize',15);
xlabel('in situ','interpreter','latex','fontsize',15);
ylabel('wamflow','interpreter','latex','fontsize',15);
title('Time-Averaged Velocity (m/s)','interpreter','latex','fontsize',15);
h2 = legend('$|\vec{u}|$','$u$','$v$','interpreter','latex','fontsize',22)
set(h2,'orientation','vertical') 

sname =['WAMvsinsitu_res',geoname,'_awin',num2str(awin),'_ovr',(num2str(ovr)),'_scatter_mean'];
    print([Tinfo.figfolder,geoname,'\',oname,'\',sname],'-dpng')

    
