% function calc_directional_spec_iter_press(Tinfo,hyd)
%This code will read the gridded camera data and then use that infromation
%to determine the two-dimensional frequency directional spectrum

%% STEP 1: Clear All
close all
clc
clear all
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
% addpath(genpath('E:\codes\insitu'))
addpath(genpath('E:\codes\trc_lab_experiment\toolbox')) 

datapath = 'E:\';
fssubfolder = datestr(date,'yy-mm-dd');
figfolder   = [datapath,'figures\methods_manuscript\results\dirspec_cond_comp\',fssubfolder,'\'];
eval(['!mkdir ',figfolder])
% 
% %% Input trial details
% 

Tinfo.Hs = 0.25;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
Tinfo.filt = 1;
if Tinfo.Hs == 0.25
    spread = [0 10 20 30 40];
elseif  Tinfo.Hs  == 0.3
    spread = [0 20 30 40];
end

inHz = 100;
geom =1;
frange = [0.25 1.2];% freq to calc Sd over
complags = 0; %base selection off computed lags
array =  1; %surfzone array  wg: 1, inner shelf array wg: 2

for iT = 1:length(spread)
%% Get trail details

Tinfo.spread = spread(iT);

Tinfo = trial_files(Tinfo);
% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);
% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);
[camera.time] = cam_time(Tinfo);

%% Load data

[~,wg] = load_insitu(Tinfo,array);
IDdef.depth = Tinfo.tide;
z = wg.z;
XL       = wg.xyz(1,:);
YL       = wg.xyz(2,:);
ZL       = mean(wg.z,1);%press.xyz(3,:);

%% STEP 4: Select Points
if complags == 1
    isel = [];
    lagx=[];
    lagy=[];
    lag = [];
    insnum = 5;
    C = nchoosek([1:size(wg.z,2)],insnum);
    for i = 1:size(C,1)
        xt = XL(C(i,:));
        yt = YL(C(i,:));
        count = 1;
        for j = 1:insnum-1
            for k = j+1:insnum
                lx(count) = xt(j)-xt(k);
                ly(count) = yt(j)-yt(k);
                lg(count)=sqrt(lx(count)^2+ly(count)^2);
                count = count+1;
            end
        end
        [v,w]=unique(round(lg,1));
        if length(w)==length(lg)
            if length(unique(round(xt)))>1
                isel = [isel; C(i,:)];
                lag = [lag; lg];
                lagx = [lagx; lx];
                lagy = [lagy; lx];
            end
        end
    end
    clear lx ly lg
    
    
    [~,ilag] = max(std(lag,[],2));
    ix = isel(ilag,:); % pick max lag
elseif complags == 0
    % pick sensors:
    if geom == 1
        sensorpick = [4  12  13  10  6 3 1]; %  BEST SO FAR
    elseif geom == 2
        sensorpick = [2 4  5 7 11  12  14]; %  BEST SO FAR 2?
    elseif geom == 3 
        sensorpick = [4 12 11 15 5 2 9]; % try this one next
    elseif geom == 4
        sensorpick = [4  12  13  10  6 3 1 15]; %  BEST SO FAR
    elseif geom == 5
        sensorpick = [2 4  5 7 11  12  14  9]; %  BEST SO FAR
    end
    
    for ip = 1:length(sensorpick)
        sensorname = ['wg',num2str(sensorpick(ip),'%02.f')];
        %     ix(ip) = ismember(sensorname,wg.loc);
        ix(ip)  = find(strcmp(sensorname,wg.loc));
    end
    % ix = [2 4 7 13 10 11];
end
%% SM defaults

SM.dirs       = [0:1:359];
SM.freqs      = [0.01:0.01:3];%[0.01:0.01:3];
SM.S          = zeros(300,360);
EP.nfft       = inHz*(2^6); %  32 sec window, previously 2^5
EP.dres       = 360;
EP.method     = 'EMLM';%'EMLM' 'IMLM' 'EMEP' 'BDM'

IDdef.datatypes  = {'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev'};% 'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev'}

IDdef.fs         = inHz;
% IDdef.depth      = 0.1664;


%%

xd = XL(ix);
yd = YL(ix);
zd = ZL(ix);

ID              = IDdef;
tpass = 1/10;
insHz = 100;
for ip = 1:length(ix)
    ztemp = detrend(pl64ta(z(:,ix(ip)),insHz*tpass));
    ID.data(:,ip) = ztemp;
%     ID.data         = z(:,ix(ip)); % 6000 x 14
end
N             = size(xd,2);
ID.layout     = [xd;yd; zeros(1,N)]; % 3 x 14
options = {'PLOTTYPE',0};
[Smout,EPout] = dirspec(ID,SM,EP,options);

dtheta = diff(Smout.dirs(1:2));
df      = diff(Smout.freqs(1:2));
Spec = Smout.S;

S(:,:,iT) =Spec;
dirs = Smout.dirs;
freqs = Smout.freqs;
dtheta        = abs(dirs(2)-dirs(1));
df            = abs(freqs(2)-freqs(1));
% dfrange 

[~,ilow] = nanmin(abs(freqs-frange(1)));
[~,ihigh] = nanmin(abs(freqs-frange(2)));
Sf(:,iT)      = nansum(Spec,2)*dtheta;
% Sd(:,iT)      = nansum(Spec,1)*df;
for f = 1:size(S,2)
    Sd(f,iT) = trapz(freqs(ilow:ihigh),Spec(ilow:ihigh,f));
end

%% comput with pwelch
sensorpickautospec = [4]; % try this one next
sensorname = ['wg',num2str(sensorpickautospec,'%02.f')];
ixauto  = find(strcmp(sensorname,wg.loc));
ztemp = detrend(pl64ta(z(:,ix(ip)),insHz*tpass));
[Sauto(:,iT),fauto]   = pwelch(ztemp,insHz*2^6,insHz*2^5,[],insHz,'ConfidenceLevel',0.95); % compute spectra

clear ID SM options Smout EPout
end

% specsave = [Tinfo.savefolder,'dirspec_wg'];
% eval(['save -v7.3 ',specsave,' dirs',' freqs',' Savg',' Spec'])


%% find param
% centerspec = 0;
% [~,ilow] = nanmin(abs(freqs-0.25));
% [~,ihigh] = nanmin(abs(freqs-1.2));
for iT = 1:length(spread)
%     for f = 1:size(S,2)
%         Sd(f,iT) = trapz(freqs(ilow:ihigh),S(ilow:ihigh,f,iT));
%     end 
    Tinfo.spread = spread(iT);
     [dirsprd_fit(iT),Sfit(iT),dirang(iT)] = fit_cos2s(Tinfo,dirs,abs(Sd(:,iT))');
%      [dirsprd(iT),Sfit(iT),dirang(iT)] = fit_cos2s_angsprd(Tinfo,dirs,abs(Sd(:,iT))');
        %  this one has a different set of iterations
    rrange = [-pi/2 pi/2]; 
    rad = deg2rad(dirs-180);
    [val,idmin] = min(abs(rad-rrange(1)));
    [val,idmax] = min(abs(rad-rrange(2)));
    rad = rad(idmin:idmax);
    Srad = Sd(idmin:idmax,iT)';
    dirsprd_obs(iT) = calc_spread(rad,Srad);
end

if complags ==  1
    psname = ['E:\data\processed\conditions\','dirspec_cond_comp\dirspec_wg_Hs',num2str(100*Tinfo.Hs),'m_instno',num2str(length(ix)),'_',EP.method,'_innershelf.mat'];
elseif complags == 0
    psname = ['E:\data\processed\conditions\','dirspec_cond_comp\dirspec_wg_Hs',num2str(100*Tinfo.Hs),'m_geom',num2str(geom),'_',EP.method,'.mat'];
end

eval(['save -v7.3 ',psname,' sensorpick',' spread',' S',' freqs',' dirs',' Sf',' Sd',' dirsprd_fit',' dirsprd_obs',' Sfit',' dirang',' Sauto',' fauto']);

%% Try DIWASP

S = abs(S);
fsmall        = [1:0.01:3];
fsmall4       = 10^-2*fsmall.^-4;

figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
semilogy(freqs,Sf(:,1),'k','linewidth',2);
hold on
if Tinfo.Hs  == 0.25
    semilogy(freqs,Sf(:,2),'r','linewidth',2);
    semilogy(freqs,Sf(:,3),'b','linewidth',2);
    semilogy(freqs,Sf(:,4),'g','linewidth',2);
    semilogy(freqs,Sf(:,5),'m','linewidth',2);
    semilogy(fauto,Sauto(:,1),'k','linewidth',1,'LineStyle','-.');
    semilogy(fauto,Sauto(:,2),'r','linewidth',1,'LineStyle','-.');
    semilogy(fauto,Sauto(:,3),'b','linewidth',1,'LineStyle','-.');
    semilogy(fauto,Sauto(:,4),'g','linewidth',1,'LineStyle','-.');
    semilogy(fauto,Sauto(:,5),'m','linewidth',1,'LineStyle','-.');
elseif  Tinfo.Hs  == 0.3
    semilogy(freqs,Sf(:,2),'b','linewidth',2);
    semilogy(freqs,Sf(:,3),'g','linewidth',2);
    semilogy(freqs,Sf(:,4),'m','linewidth',2);
end
semilogy(fsmall,fsmall4,'b','linewidth',2);
if Tinfo.Hs == 0.25
    h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 10^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
elseif Tinfo.Hs == 0.3
    h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
end
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
ylim([10^-4.5 10^-1])
xlim([0 3]);
h1=gca;
set(h1,'fontsize',15);
% title(EP.method,'interpreter','latex')

subplot(212)
semilogy(dirs,Sd(:,1),'k','linewidth',2);
hold on
if Tinfo.Hs == 0.25
    semilogy(dirs,Sd(:,2),'r','linewidth',2);
    semilogy(dirs,Sd(:,3),'b','linewidth',2);
    semilogy(dirs,Sd(:,4),'g','linewidth',2);
    semilogy(dirs,Sd(:,5),'m','linewidth',2);
elseif Tinfo.Hs == 0.3
    semilogy(dirs,Sd(:,2),'b','linewidth',2);
    semilogy(dirs,Sd(:,3),'g','linewidth',2);
    semilogy(dirs,Sd(:,4),'m','linewidth',2);
end
grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([0 360]);
h1=gca;
set(h1,'fontsize',15);

% figfolder   = [Tinfo.figfolder,'dirspec/'];
% eval(['!mkdir ',figfolder])
Sname1 = [figfolder,'dirspec_wg_Hs',num2str(100*Tinfo.Hs),'m_instno',num2str(length(ix)),'_',EP.method];
print(Sname1,'-dpng')

%%

S = abs(S);
fsmall        = [1:0.01:3];
fsmall4       = 10^-2*fsmall.^-4;

figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
plot(freqs,Sf(:,1),'k','linewidth',1);
hold on
if Tinfo.Hs  == 0.25
    plot(freqs,Sf(:,2),'r','linewidth',1);
    plot(freqs,Sf(:,3),'b','linewidth',1);
    plot(freqs,Sf(:,4),'g','linewidth',1);
    plot(freqs,Sf(:,5),'m','linewidth',1);
    plot(fauto,Sauto(:,1),'k','linewidth',1,'LineStyle','-.');
    plot(fauto,Sauto(:,1),'r','linewidth',1,'LineStyle','-.');
    plot(fauto,Sauto(:,1),'b','linewidth',1,'LineStyle','-.');
    plot(fauto,Sauto(:,1),'g','linewidth',1,'LineStyle','-.');
    plot(fauto,Sauto(:,1),'m','linewidth',1,'LineStyle','-.');
elseif  Tinfo.Hs  == 0.3
    plot(freqs,Sf(:,2),'b','linewidth',1);
    plot(freqs,Sf(:,3),'g','linewidth',1);
    plot(freqs,Sf(:,4),'m','linewidth',1);
end
plot(fsmall,fsmall4,'b','linewidth',2);
if Tinfo.Hs == 0.25
    h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 10^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
elseif Tinfo.Hs == 0.3
    h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
end
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
ylim([10^-4.5 10^-1])
xlim([0 3]);
h1=gca;
set(h1,'fontsize',15);
% title(EP.method,'interpreter','latex')

subplot(212)
plot(dirs,Sd(:,1),'k','linewidth',1);
hold on
if Tinfo.Hs == 0.25
    plot(dirs,Sd(:,2),'r','linewidth',2);
    plot(dirs,Sd(:,3),'b','linewidth',2);
    plot(dirs,Sd(:,4),'g','linewidth',2);
    plot(dirs,Sd(:,5),'m','linewidth',2);
elseif Tinfo.Hs == 0.3
    plot(dirs,Sd(:,2),'b','linewidth',2);
    plot(dirs,Sd(:,3),'g','linewidth',2);
    plot(dirs,Sd(:,4),'m','linewidth',2);
end
grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([0 360]);
h1=gca;
set(h1,'fontsize',15);

% figfolder   = [Tinfo.figfolder,'dirspec/'];
% eval(['!mkdir ',figfolder])
Sname1 = [figfolder,'dirspec_wg_Hs',num2str(100*Tinfo.Hs),'m_instno',num2str(length(ix)),'_',EP.method];

if complags ==  1
    sname = [figfolder,'Sf_Sd_wg_Hs',num2str(100*Tinfo.Hs),'m_instno',num2str(length(ix)),'_',EP.method];
elseif complags == 0
    sname = [figfolder,'Sf_Sd_wg_Hs',num2str(100*Tinfo.Hs),'m_geom',num2str(geom),'_',EP.method];
end
print(sname,'-dpng')

%%
% Trun = 1;
figure('units','inches','position',[1 1 17 5],'color','w');

conf = [0.9631 1.03156];

ax0 = axes('Position',[0.03 -2 0.28 0.96]);
hcb = pcolor(S(:,:,1))
caxis([10^(-8) 10^(-3.5)])
set(gca,'ColorScale','log')
colormap(jet)
cb = colorbar('Position',[0.035 0.095 0.022 0.65])
cb.Ruler.MinorTick = 'on';
h1=gca;
set(h1,'fontsize',18)


ax1 = axes('Position',[0.09 0.08 0.24 0.75]);
[h,c] = polarPcolor(freqs(1:100),[dirs],S((1:100),:,1),'ncolor',10,'colormap','jet','Ncircles',4,'Nspokes',9,'colBar',0);%,'RtickLabel',['0' '0.5' '1.0' '1.5' '2.0'])%,'typeRose','default')
caxis([10^(-5) 10^(-3.5)])

ax2 = axes('Position',[0.39 0.08 0.24 0.75]);
if Tinfo.Hs == 0.25
    [h,c] = polarPcolor(freqs(1:100),[dirs],S((1:100),:,4),'ncolor',10,'colormap','jet','Ncircles',4,'Nspokes',9,'colBar',0);%,'RtickLabel',['0' '0.5' '1.0' '1.5' '2.0'])%,'typeRose','default')
elseif Tinfo.Hs ==0.3
    [h,c] = polarPcolor(freqs(1:100),[dirs],S((1:100),:,3),'ncolor',10,'colormap','jet','Ncircles',4,'Nspokes',9,'colBar',0);%,'RtickLabel',['0' '0.5' '1.0' '1.5' '2.0'])%,'typeRose','default')
end
caxis([10^(-5) 10^(-3.5)])

ax3 = axes('Position',[0.72 0.15 0.26 0.83]);
% ADD THE CONFIDENCE INTERVAL
semilogy(dirs-180,Sd(:,1),'k','linewidth',1);
hold on
if Tinfo.Hs == 0.25
    semilogy(dirs-180,Sd(:,2),'r','linewidth',2);
    semilogy(dirs-180,Sd(:,3),'b','linewidth',2);
    semilogy(dirs-180,Sd(:,4),'g','linewidth',2);
    semilogy(dirs-180,Sd(:,5),'m','linewidth',2);
elseif Tinfo.Hs == 0.3
    semilogy(dirs-180,Sd(:,2),'b','linewidth',2);
    semilogy(dirs-180,Sd(:,3),'g','linewidth',2);
    semilogy(dirs-180,Sd(:,4),'m','linewidth',2);
end
xlabel('$\theta~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{\eta\eta}(\theta)~(\mathrm{m^{2}/^{\circ}})$','interpreter','latex')
xlim([-180 180]);
ylim([10^(-6) 10^(-2)])
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'YScale','log')
set(h1,'fontsize',15);
if Tinfo.Hs == 0.25
    h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 10^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
elseif Tinfo.Hs == 0.3
    h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
end
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
text(-340-700,10^(-2.3),'(a)$\sigma_{\theta} = 0^{\circ}$','interpreter','latex','color','k','fontsize',22)
text(-130-530,10^(-2.3),'(b) $\sigma_{\theta} = 30^{\circ}$','interpreter','latex','color','k','fontsize',22)
text(95-270,10^(-2.3),'(c)','interpreter','latex','color','k','fontsize',22)
text(-410-750,10^(-2.7),'$S_{\eta\eta}(f,\theta)$','interpreter','latex','color','k','fontsize',22)


if complags ==  1
    sname = [figfolder,'dirspec_wg_Hs',num2str(100*Tinfo.Hs),'m_instno',num2str(length(ix)),'_',EP.method];
elseif complags == 0
    sname = [figfolder,'dirspec_wg_Hs',num2str(100*Tinfo.Hs),'m_geom',num2str(geom),'_',EP.method];
end
print(sname,'-dpng')