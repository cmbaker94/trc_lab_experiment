% Plot directional spectra properties and compare for multiple trials.

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\codes\insitu'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

spread = [0];
inst = 14;
iter = 200;
% rad2deg  = 180/pi;
frange = [0.3 0.8];%[0.25 1.25];
dirrange = [0 359];

for isprd = 1:length(spread)

% Trial info
Tinfo.Hs = 0.3;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
Tinfo.spread = spread(isprd);

Tinfo = trial_files(Tinfo);

% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

%% STEP 2: Define figure folders

% figure folder
fssubfolder = datestr(date,'yy-mm-dd');
Tinfo.figfolder   = [Tinfo.datapath,'figures\wave_statistics\directional\',fssubfolder,'\'];

% make figure folders
eval(['mkdir ',Tinfo.datapath,'figures\wave_statistics\directional\;'])
eval(['mkdir ',Tinfo.figfolder])


% STEP 1: Load data

geninfo = 'regx28-33_regy-10-4_dx50_dy50';
Sftheta = csvread([Tinfo.savefolder,'spec_cam_',geninfo,'_grid_rand_',num2str(inst),'inst_',num2str(iter),'iter.csv']);
freq = csvread([Tinfo.savefolder,'freq_cam_',geninfo,'_grid_rand_',num2str(inst),'inst_',num2str(iter),'iter.csv'])';
dir = csvread([Tinfo.savefolder,'dirs_cam_',geninfo,'_grid_rand_',num2str(inst),'inst_',num2str(iter),'iter.csv'])';
cam = calc_dir_prop(Sftheta,freq,dir,frange,dirrange);
[cam.sprd,cam.S] = fit_cos2s(cam.dir,cam.Sd)
cam.sprd
clear Sftheta freq dir

%% lidar

geninfo = 'regx28-33_regy-10-4_dx50_dy50';
Sftheta = csvread([Tinfo.savefolder,'spec_lid_',geninfo,'_grid_rand_',num2str(inst),'inst_',num2str(iter),'iter.csv']);
freq = csvread([Tinfo.savefolder,'freq_lid_',geninfo,'_grid_rand_',num2str(inst),'inst_',num2str(iter),'iter.csv'])';
dir = csvread([Tinfo.savefolder,'dirs_lid_',geninfo,'_grid_rand_',num2str(inst),'inst_',num2str(iter),'iter.csv'])';
lid = calc_dir_prop(Sftheta,freq,dir,frange,dirrange);
clear Sftheta freq dir
[lid.sprd,lid.S] = fit_cos2s(lid.dir,lid.Sd);
lid.sprd

%%

Sftheta = csvread([Tinfo.savefolder,'spec_press_press_sz_5inst_200iter.csv']);
freq = csvread([Tinfo.savefolder,'freq_press_press_sz_5inst_200iter.csv'])';
dir = csvread([Tinfo.savefolder,'dirs_press_press_sz_5inst_200iter.csv'])';
pg = calc_dir_prop(Sftheta,freq,dir,frange,dirrange);
clear Sftheta
[pg.sprd,pg.S] = fit_cos2s(pg.dir,pg.Sd);
pg.sprd

%%

puvfile     = ['E:\data\processed\insitu\',Tinfo.sz.tdate,'\PUV_p06.csv'];
puv06 = read_puv_csv(puvfile);
puvfile     = ['E:\data\processed\insitu\',Tinfo.sz.tdate,'\PUV_p11.csv'];
eval(['puv11 = read_puv_csv(puvfile);'])
% freqrange = [0.2 1.25];
% spec = T0puv06.th;
% [emean] = calc_energy_weighted_mean(spec,freq,freqrange)

end

%%
S = 85;%[1:65]%4, 7, 16, 65];
TH = 0;
spreading = gamma(S+1)/2/sqrt(pi)./gamma(S+1/2).*cos(TH/2).^(2*S);
spreading2 = ((2.^(2*S-1))/pi).*((gamma(S+1)).^2)./gamma(2*S+1).*cos(TH/2).^(2*S);


theta = [-(pi/2):0.1:pi/2]%[-4:0.01:4];%deg2rad([-40:1:40]);
thetad = rad2deg(theta)+180;
for i = 1:length(S)
%     sprd(i,:) = gamma(S(i)+1)/2/sqrt(pi)./gamma(S(i)+1/2).*cos(theta/2).^(2*S(i));
    sdist(i,:) = ((2.^(2*S(i)-1))/pi).*((gamma(S(i)+1)).^2)./gamma(2*S(i)+1).*cos(theta/2).^(2*S(i));
end

plot(theta,sdist(1,:))
% eqn from Squire and Montiel, 2016
pick = 1;
left = trapz(theta,sdist(pick,:).*cos(theta))^2;
right = trapz(theta,sdist(pick,:).*sin(theta))^2;
sprd = sqrt(2*(1-sqrt(left+right)));
test = rad2deg(sprd)


%% find floor

% floor =sort (cam.Sd,2, 'ascend');
% cam.floor = mean(floor(1:30));
% 
% floor =sort (lid.Sd,2, 'ascend');
% lid.floor = mean(floor(1:30));
% 
% floor =sort (pg.Sd,2, 'ascend');
% pg.floor = mean(floor(1:30));

%%


% rrange = [-pi/2 pi/2]; 
% rad = deg2rad(dir-180);
% [val,idmin] = min(abs(rad-rrange(1)));
% [val,idmax] = min(abs(rad-rrange(2)));
% rad = rad(idmin:idmax);
% Srad = Sd(idmin:idmax);
% fascend =sort (Srad,2, 'ascend');
% floor = mean(fascend(1:30));
% Srad = Srad-floor;
% 
% intS = trapz(rad,Srad);
% Srnorm = Srad/intS;
% 
% 
% S = 1:85;
% TH = 0;
% pltflg = 0;
% for i = 1:length(S)
%     sdist = ((2.^(2*S(i)-1))/pi).*((gamma(S(i)+1)).^2)./gamma(2*S(i)+1).*cos(rad/2).^(2*S(i));
% %     serror(i) = fit1dregerror(rad,sdist,rad,Srnorm,pltflg);
%     y0 = sdist;
%     y1 = Srnorm;
% %     abs_dy(i,:) = abs(y0-y1);   % absolute error
% %     relerr(i,:) = abs(y0-y1)./y0 ;  % relative error
% %     pererr(i,:) = abs(y0-y1)./y0*100 ;   % percentage error
%     mean_err(i) = mean(abs(y0-y1)) ;    % mean absolute error
%     MSE(i) = mean((y0-y1).^2) ;        % Mean square error
%     RMSE(i) = sqrt(mean((y0-y1).^2)) ; % Root mean square error
% end
% [val,iS] = nanmin(RMSE);
% Sfit = S(iS);
% Sdist = ((2.^(2*Sfit-1))/pi).*((gamma(Sfit+1)).^2)./gamma(2*Sfit+1).*cos(rad/2).^(2*Sfit);
% 
% 
% figure
% plot(rad,Srnorm)
% hold on
% plot(rad,Sdist)
% 
% pick = 4;
% left = trapz(rad,Sdist.*cos(rad))^2;
% right = trapz(rad,Sdist.*sin(rad))^2;
% sprd = sqrt(2*(1-sqrt(left+right)));
% dirsprd = rad2deg(sprd)



%%
fsmall        = [1:0.01:3];
fsmall4       = 10^-2.5*fsmall.^-4;
figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
semilogy(cam.freq,cam.Sf,'k','linewidth',1.5);
hold on
semilogy(lid.freq,lid.Sf,'r','linewidth',1.5);
semilogy(pg.freq(1:151),pg.Sf(1:151),'b','linewidth',1.5);
semilogy(puv11.freq,puv11.SSE,'g','linewidth',1.5);
% semilogy(fsmall,fsmall4,'b','linewidth',2);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
ylim([10^-4.5 10^-1.5])
xlim([0 3]);
h1=gca;
set(h1,'fontsize',15);
h2 = legend('Stereo','Lidar','Press','puv')
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');


subplot(212)
semilogy(cam.dir,cam.Sd,'k','linewidth',1.5);
hold on
semilogy(lid.dir,lid.Sd,'r','linewidth',1.5);
semilogy(pg.dir,pg.Sd,'b','linewidth',1.5);
plot([min(cam.dir) max(cam.dir)],[cam.floor cam.floor],'k','linestyle','-.')
plot([min(lid.dir) max(lid.dir)],[lid.floor lid.floor],'r','linestyle','-.')
plot([min(pg.dir) max(pg.dir)],[pg.floor pg.floor],'b','linestyle','-.')
grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([0 360]);
h1=gca;
set(h1,'fontsize',15);
sname = [Tinfo.figfolder,'Sf_Sd'];
print(sname,'-dpng')

%%
% freqrange = [0.2 1.2];
% dirrange = [90 270];
% 
% % def calc_energy_weighted_mean(spec,freq,freqrange):
% %     idmin = min(range(len(freq)), key=lambda i: abs(freq[i]-freqrange[0]))
% %     idmax = min(range(len(freq)), key=lambda i: abs(freq[i]-freqrange[1]))
% %     specint = np.trapz(spec[idmin:idmax],x=freq[idmin:idmax])
% %     emean = specint/(freq[idmax]-freq[idmin])
% %     return emean
%     
% % [val, ifmin] = min(abs(freq-freqrange(1)));
% % [val, ifmax] = min(abs(freq-freqrange(2)));
% % [val, idmin] = min(abs(dir-dirrange(1)));
% % [val, idmax] = min(abs(dir-dirrange(2)));
% % specf = nanmean(Sftheta(:,idmin:idmax),2);
% % spec = (specf(ifmin:ifmax).*freq(ifmin:ifmax))/(freq(ifmax)-freq(ifmin));
% spec = sig_2;
% [val, ifmin] = min(abs(freq-freqrange(1)));
% [val, ifmax] = min(abs(freq-freqrange(2)));
% % specf = nanmean(sig_2(ifmin:ifmax),2);
% test = trapz(spec(ifmin:ifmax),freq(ifmin:ifmax)')/(freq(ifmax)-freq(ifmin));
% 
% % specint = trapz(

%%
fsmall        = [1:0.01:3];
fsmall4       = 10^-2.5*fsmall.^-4;
figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
semilogy(T0cam.freq,T0cam.Sf,'k','linewidth',1.5);
hold on
semilogy(T0lid.freq,T0lid.Sf,'k','linewidth',1.5,'linestyle','-.');
semilogy(T0pg.freq(1:151),T0pg.Sf(1:151),'k','linewidth',1.5,'linestyle',':');
semilogy(T0cam.freq,T0cam.Sf,'k','linewidth',1.5);
hold on
semilogy(T20cam.freq,T20cam.Sf,'r','linewidth',1.5);
semilogy(T30cam.freq,T30cam.Sf,'Color',[0.1 0.7 0.15],'linewidth',1.5);
semilogy(T40cam.freq,T40cam.Sf,'b','linewidth',1.5);
semilogy(T0lid.freq,T0lid.Sf,'k','linewidth',1.5,'linestyle','-.');
semilogy(T20lid.freq,T20lid.Sf,'r','linewidth',1.5,'linestyle','-.');
semilogy(T30lid.freq,T30lid.Sf,'Color',[0.1 0.7 0.15],'linewidth',1.5,'linestyle','-.');
semilogy(T40lid.freq,T40lid.Sf,'b','linewidth',1.5,'linestyle','-.');
semilogy(T0pg.freq(1:151),T0pg.Sf(1:151),'k','linewidth',1.5,'linestyle',':');
semilogy(T20pg.freq(1:151),T20pg.Sf(1:151),'r','linewidth',1.5,'linestyle',':');
semilogy(T30pg.freq(1:151),T30pg.Sf(1:151),'Color',[0.1 0.7 0.15],'linewidth',1.5,'linestyle',':');
semilogy(T40pg.freq(1:151),T40pg.Sf(1:151),'b','linewidth',1.5,'linestyle',':');
% semilogy(fsmall,fsmall4,'b','linewidth',2);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
ylim([10^-4.5 10^-1.5])
xlim([0 3]);
h1=gca;
set(h1,'fontsize',15);
h2 = legend('Stereo','Lidar','Press','$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');


subplot(212)
semilogy(T0cam.dir,T0cam.Sd,'k','linewidth',1.5);
hold on
semilogy(T20cam.dir,T20cam.Sd,'r','linewidth',1.5);
semilogy(T30cam.dir,T30cam.Sd,'Color',[0.1 0.7 0.15],'linewidth',1.5);
semilogy(T40cam.dir,T40cam.Sd,'b','linewidth',1.5);
semilogy(T0lid.dir,T0lid.Sd,'k','linewidth',1.5,'linestyle','-.');
semilogy(T20lid.dir,T20lid.Sd,'r','linewidth',1.5,'linestyle','-.');
semilogy(T30lid.dir,T30lid.Sd,'Color',[0.1 0.7 0.15],'linewidth',1.5,'linestyle','-.');
semilogy(T40lid.dir,T40lid.Sd,'b','linewidth',1.5,'linestyle','-.');
semilogy(T0pg.dir,T0pg.Sd,'k','linewidth',1.5,'linestyle',':');
semilogy(T20pg.dir,T20pg.Sd,'r','linewidth',1.5,'linestyle',':');
semilogy(T30pg.dir,T30pg.Sd,'Color',[0.1 0.7 0.15],'linewidth',1.5,'linestyle',':');
semilogy(T40pg.dir,T40pg.Sd,'b','linewidth',1.5,'linestyle',':');
grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([0 360]);
h1=gca;
set(h1,'fontsize',15);
sname = [Tinfo.figfolder,'Sf_Sd'];
print(sname,'-dpng')

%%

figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
plot(T0cam.freq,T0cam.Sf,'k','linewidth',1.5);
hold on
plot(T20cam.freq,T20cam.Sf,'r','linewidth',1.5);
plot(T30cam.freq,T30cam.Sf,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40cam.freq,T40cam.Sf,'b','linewidth',1.5);
% plot(T50cam.freq,T50cam.Sf,'m','linewidth',1.5);
plot(T0lid.freq,T0lid.Sf,'k','linewidth',1.5,'linestyle','-.');
plot(T20lid.freq,T20lid.Sf,'r','linewidth',1.5,'linestyle','-.');
plot(T30lid.freq,T30lid.Sf,'Color',[0.1 0.7 0.15],'linewidth',1.5,'linestyle','-.');
plot(T40lid.freq,T40lid.Sf,'b','linewidth',1.5,'linestyle','-.');
plot(T0pg.freq,T0pg.Sf,'k','linewidth',1.5,'linestyle',':');
plot(T20pg.freq,T20pg.Sf,'r','linewidth',1.5,'linestyle',':');
plot(T30pg.freq,T30pg.Sf,'Color',[0.1 0.7 0.15],'linewidth',1.5,'linestyle',':');
plot(T40pg.freq,T40pg.Sf,'b','linewidth',1.5,'linestyle',':');
% semilogy(fsmall,fsmall4,'b','linewidth',2);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
% ylim([0 10^-4])
xlim([0 3]);
h1=gca;
set(h1,'fontsize',15);
h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');


subplot(212)
plot(T0cam.dir,T0cam.Sd,'k','linewidth',1.5);
hold on
plot(T20cam.dir,T20cam.Sd,'r','linewidth',1.5);
plot(T30cam.dir,T30cam.Sd,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40cam.dir,T40cam.Sd,'b','linewidth',1.5);
% plot(T50cam.dir,T50cam.Sd,'m','linewidth',1.5);
plot(T0lid.dir,T0lid.Sd,'k','linewidth',1.5,'linestyle','-.');
plot(T20lid.dir,T20lid.Sd,'r','linewidth',1.5,'linestyle','-.');
plot(T30lid.dir,T30lid.Sd,'Color',[0.1 0.7 0.15],'linewidth',1.5,'linestyle','-.');
plot(T40lid.dir,T40lid.Sd,'b','linewidth',1.5,'linestyle','-.');
plot(T0pg.dir,T0pg.Sd,'k','linewidth',1.5,'linestyle',':');
plot(T20pg.dir,T20pg.Sd,'r','linewidth',1.5,'linestyle',':');
plot(T30pg.dir,T30pg.Sd,'Color',[0.1 0.7 0.15],'linewidth',1.5,'linestyle',':');
plot(T40pg.dir,T40pg.Sd,'b','linewidth',1.5,'linestyle',':');
grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([0 360]);
ylim([0 0.000015])
h1=gca;
set(h1,'fontsize',15);
sname = [Tinfo.figfolder,'Sf_Sd_lin'];
print(sname,'-dpng')


%%
figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
%plot(freq,th_1,'k','linewidth',1.5);
plot(T0cam.freq,T0cam.th_2,'k','linewidth',1.5);
hold on
plot(T0lid.freq,T0lid.th_2,'k','linewidth',1.5,'linestyle','-.');
plot(T0pg.freq,T0pg.th_2,'k','linewidth',1.5,'linestyle',':');
plot(T0cam.freq,T0cam.th_2,'k','linewidth',1.5);
plot(T20cam.freq,T20cam.th_2,'r','linewidth',1.5);
plot(T30cam.freq,T30cam.th_2,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40cam.freq,T40cam.th_2,'b','linewidth',1.5);
% plot(T50cam.freq,T50cam.th,'m','linewidth',1.5);
plot(T0lid.freq,T0lid.th_2,'k','linewidth',1.5,'linestyle','-.');
plot(T20lid.freq,T20lid.th_2,'r','linewidth',1.5,'linestyle','-.');
plot(T30lid.freq,T30lid.th_2,'Color',[0.1 0.7 0.15],'linewidth',1.5,'linestyle','-.');
plot(T40lid.freq,T40lid.th_2,'b','linewidth',1.5,'linestyle','-.');
plot(T0pg.freq,T0pg.th_2,'k','linewidth',1.5,'linestyle',':');
plot(T20pg.freq,T20pg.th_2,'r','linewidth',1.5,'linestyle',':');
plot(T30pg.freq,T30pg.th_2,'Color',[0.1 0.7 0.15],'linewidth',1.5,'linestyle',':');
plot(T40pg.freq,T40pg.th_2,'b','linewidth',1.5,'linestyle',':');
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$\theta_2~(^{\circ})$','interpreter','latex')
h2 = legend('Stero','Camera','Lidar','$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
%ylim([10^-4.5 10^-1.5])
xlim([0 1.5]);
ylim([-90 90]);
h1=gca;
set(h1,'fontsize',15);

subplot(212)
%plot(freq,sig_1,'k','linewidth',1.5);
hold on
plot(T0cam.freq,T0cam.sig_2,'k','linewidth',1.5);
plot(T20cam.freq,T20cam.sig_2,'r','linewidth',1.5);
plot(T30cam.freq,T30cam.sig_2,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40cam.freq,T40cam.sig_2,'b','linewidth',1.5);
% plot(T50cam.freq,T50cam.sig,'m','linewidth',1.5);
plot(T0lid.freq,T0lid.sig_2,'k','linewidth',1.5,'linestyle','-.');
plot(T20lid.freq,T20lid.sig_2,'r','linewidth',1.5,'linestyle','-.');
plot(T30lid.freq,T30lid.sig_2,'Color',[0.1 0.7 0.15],'linewidth',1.5,'linestyle','-.');
plot(T40lid.freq,T40lid.sig_2,'b','linewidth',1.5,'linestyle','-.');
plot(T0pg.freq,T0pg.sig_2,'k','linewidth',1.5,'linestyle',':');
plot(T20pg.freq,T20pg.sig_2,'r','linewidth',1.5,'linestyle',':');
plot(T30pg.freq,T30pg.sig_2,'Color',[0.1 0.7 0.15],'linewidth',1.5,'linestyle',':');
plot(T40pg.freq,T40pg.sig_2,'b','linewidth',1.5,'linestyle',':');
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$\sigma_{\theta,2}~(^{\circ})$','interpreter','latex')
xlim([0 1.5]);
ylim([0 45]);
h1=gca;
set(h1,'fontsize',15);
% h2 = legend('Stero','Camera','Lidar','$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
% set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
sname = [Tinfo.figfolder,'theta_sig'];
print(sname,'-dpng')

%%

figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
%plot(freq,th_1,'k','linewidth',1.5);
hold on
plot(T0cam.freq,T0cam.th_1,'k','linewidth',1.5);
plot(T20cam.freq,T20cam.th_1,'r','linewidth',1.5);
plot(T30cam.freq,T30cam.th_1,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40cam.freq,T40cam.th_1,'b','linewidth',1.5);
% plot(T50cam.freq,T50cam.th,'m','linewidth',1.5);
plot(T0lid.freq,T0lid.th_1,'k','linewidth',1.5,'linestyle','-.');
plot(T20lid.freq,T20lid.th_1,'r','linewidth',1.5,'linestyle','-.');
plot(T30lid.freq,T30lid.th_1,'Color',[0.1 0.7 0.15],'linewidth',1.5,'linestyle','-.');
plot(T40lid.freq,T40lid.th_1,'b','linewidth',1.5,'linestyle','-.');
plot(T0pg.freq,T0pg.th_1,'k','linewidth',1.5,'linestyle',':');
plot(T20pg.freq,T20pg.th_1,'r','linewidth',1.5,'linestyle',':');
plot(T30pg.freq,T30pg.th_1,'Color',[0.1 0.7 0.15],'linewidth',1.5,'linestyle',':');
plot(T40pg.freq,T40pg.th_1,'b','linewidth',1.5,'linestyle',':');
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$\theta_1~(^{\circ})$','interpreter','latex')
%ylim([10^-4.5 10^-1.5])
xlim([0 1.5]);
ylim([-90 90]);
h1=gca;
set(h1,'fontsize',15);

subplot(212)
%plot(freq,sig_1,'k','linewidth',1.5);
hold on
plot(T0cam.freq,T0cam.sig_1,'k','linewidth',1.5);
plot(T20cam.freq,T20cam.sig_1,'r','linewidth',1.5);
plot(T30cam.freq,T30cam.sig_1,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40cam.freq,T40cam.sig_1,'b','linewidth',1.5);
% plot(T50cam.freq,T50cam.sig,'m','linewidth',1.5);
plot(T0lid.freq,T0lid.sig_1,'k','linewidth',1.5,'linestyle','-.');
plot(T20lid.freq,T20lid.sig_1,'r','linewidth',1.5,'linestyle','-.');
plot(T30lid.freq,T30lid.sig_1,'Color',[0.1 0.7 0.15],'linewidth',1.5,'linestyle','-.');
plot(T40lid.freq,T40lid.sig_1,'b','linewidth',1.5,'linestyle','-.');
plot(T0pg.freq,T0pg.sig_1,'k','linewidth',1.5,'linestyle',':');
plot(T20pg.freq,T20pg.sig_1,'r','linewidth',1.5,'linestyle',':');
plot(T30pg.freq,T30pg.sig_1,'Color',[0.1 0.7 0.15],'linewidth',1.5,'linestyle',':');
plot(T40pg.freq,T40pg.sig_1,'b','linewidth',1.5,'linestyle',':');
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$\sigma_{\theta,1}~(^{\circ})$','interpreter','latex')
xlim([0 1.5]);
ylim([0 45]);
h1=gca;
set(h1,'fontsize',15);
h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
sname = [Tinfo.figfolder,'theta_sig_1'];
print(sname,'-dpng')

%%

figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
plot(T0cam.freq,T0cam.Sf,'k','linewidth',1.5);
hold on
plot(T20cam.freq,T20cam.Sf,'r','linewidth',1.5);
plot(T30cam.freq,T30cam.Sf,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40cam.freq,T40cam.Sf,'b','linewidth',1.5);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
% ylim([0 10^-4])
xlim([0 3]);
h1=gca;
set(h1,'fontsize',15);
h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');


subplot(212)
plot(T0cam.dir,T0cam.Sd,'k','linewidth',1.5);
hold on
plot(T20cam.dir,T20cam.Sd,'r','linewidth',1.5);
plot(T30cam.dir,T30cam.Sd,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40cam.dir,T40cam.Sd,'b','linewidth',1.5);
grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([0 360]);
ylim([0 0.000015])
h1=gca;
set(h1,'fontsize',15);
sname = [Tinfo.figfolder,'Sf_Sd_cam'];
print(sname,'-dpng')

figure
plot(T0cam.dir,T0cam.Sd,'k','linewidth',1.5);
hold on
plot(T20cam.dir,T20cam.Sd,'r','linewidth',1.5);
plot(T30cam.dir,T30cam.Sd,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40cam.dir,T40cam.Sd,'b','linewidth',1.5);
plot(thetad,sprd(1,:),'k','LineWidth',2)
plot(thetad,sprd(2,:),'r','LineWidth',2)
plot(thetad,sprd(3,:),'b','LineWidth',2)
plot(thetad,sprd(4,:),'g','LineWidth',2)
grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([0 360]);
ylim([0 0.000015])
h1=gca;
set(h1,'fontsize',15);
sname = [Tinfo.figfolder,'Sf_Sd_cam'];
print(sname,'-dpng')

%%

figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
plot(T0cam.freq,T0cam.Sf,'k','linewidth',1.5);
hold on
plot(T20cam.freq,T20cam.Sf,'r','linewidth',1.5);
plot(T30cam.freq,T30cam.Sf,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40cam.freq,T40cam.Sf,'b','linewidth',1.5);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
% ylim([0 10^-4])
xlim([0 3]);
h1=gca;
set(h1,'fontsize',15);
h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');


subplot(212)
plot(T0cam.dir,T0cam.Sd,'k','linewidth',1.5);
hold on
plot(T20cam.dir,T20cam.Sd,'r','linewidth',1.5);
plot(T30cam.dir,T30cam.Sd,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40cam.dir,T40cam.Sd,'b','linewidth',1.5);
grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([0 360]);
ylim([0 0.000015])
h1=gca;
set(h1,'fontsize',15);
sname = [Tinfo.figfolder,'Sf_Sd_cam'];
print(sname,'-dpng')

%%

figure('units','inches','position',[1 1 6 6],'color','w');
plot(T0cam.dir,T0cam.Sd,'k','linewidth',2);
hold on
plot(T20cam.dir,T20cam.Sd,'r','linewidth',2,'LineStyle','-');
% plot(T30lid.dir,T30lid.Sd,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40cam.dir,T40cam.Sd,'b','linewidth',2,'LineStyle','-');
% grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([90 270]);
ylim([0 0.000016])
h1=gca;
set(h1,'fontsize',15);
h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 40^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');

sname = [Tinfo.figfolder,'Sf_Sd_cam_sd'];
print(sname,'-dpng')

%%

%%

figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
plot(T0pg.freq,T0pg.Sf,'k','linewidth',1.5);
hold on
plot(T20pg.freq,T20pg.Sf,'r','linewidth',1.5);
plot(T30pg.freq,T30pg.Sf,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40pg.freq,T40pg.Sf,'b','linewidth',1.5);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
% ylim([0 10^-4])
xlim([0 3]);
h1=gca;
set(h1,'fontsize',15);
h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');


subplot(212)
plot(T0pg.dir,T0pg.Sd,'k','linewidth',1.5);
hold on
plot(T20pg.dir,T20pg.Sd,'r','linewidth',1.5);
plot(T30pg.dir,T30pg.Sd,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40pg.dir,T40pg.Sd,'b','linewidth',1.5);
grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([0 360]);
ylim([0 0.000015])
h1=gca;
set(h1,'fontsize',15);
sname = [Tinfo.figfolder,'Sf_Sd_pg'];
print(sname,'-dpng')

%%

%%
figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
plot(T0cam.freq,T0cam.th_2,'k','linewidth',1.5);
hold on
plot(T20cam.freq,T20cam.th_2,'r','linewidth',1.5);
plot(T30cam.freq,T30cam.th_2,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40cam.freq,T40cam.th_2,'b','linewidth',1.5);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$\theta_2~(^{\circ})$','interpreter','latex')
h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
%ylim([10^-4.5 10^-1.5])
xlim([0 1.5]);
ylim([-90 90]);
h1=gca;
set(h1,'fontsize',15);

subplot(212)
%plot(freq,sig_1,'k','linewidth',1.5);
hold on
plot(T0cam.freq,T0cam.sig_2,'k','linewidth',1.5);
plot(T20cam.freq,T20cam.sig_2,'r','linewidth',1.5);
plot(T30cam.freq,T30cam.sig_2,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40cam.freq,T40cam.sig_2,'b','linewidth',1.5);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$\sigma_{\theta,2}~(^{\circ})$','interpreter','latex')
xlim([0 1.5]);
ylim([0 45]);
h1=gca;
set(h1,'fontsize',15);
% h2 = legend('Stero','Camera','Lidar','$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
% set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
sname = [Tinfo.figfolder,'theta_sig_cam'];
print(sname,'-dpng')

%%

figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
plot(T0lid.freq,T0lid.th_2,'k','linewidth',1.5);
hold on
plot(T20lid.freq,T20lid.th_2,'r','linewidth',1.5);
plot(T30lid.freq,T30lid.th_2,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40lid.freq,T40lid.th_2,'b','linewidth',1.5);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$\theta_2~(^{\circ})$','interpreter','latex')
h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
%ylim([10^-4.5 10^-1.5])
xlim([0 1.5]);
ylim([-90 90]);
h1=gca;
set(h1,'fontsize',15);

subplot(212)
%plot(freq,sig_1,'k','linewidth',1.5);
hold on
plot(T0lid.freq,T0lid.sig_2,'k','linewidth',1.5);
plot(T20lid.freq,T20lid.sig_2,'r','linewidth',1.5);
plot(T30lid.freq,T30lid.sig_2,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40lid.freq,T40lid.sig_2,'b','linewidth',1.5);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$\sigma_{\theta,2}~(^{\circ})$','interpreter','latex')
xlim([0 1.5]);
ylim([0 45]);
h1=gca;
set(h1,'fontsize',15);
% h2 = legend('Stero','lidera','Lidar','$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
% set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
sname = [Tinfo.figfolder,'theta_sig_lid'];
print(sname,'-dpng')

%%

%%

figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
plot(T0pg.freq,T0pg.th_2,'k','linewidth',1.5);
hold on
plot(T20pg.freq,T20pg.th_2,'r','linewidth',1.5);
plot(T30pg.freq,T30pg.th_2,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40pg.freq,T40pg.th_2,'b','linewidth',1.5);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$\theta_2~(^{\circ})$','interpreter','latex')
h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
%ylim([10^-4.5 10^-1.5])
xlim([0 1.5]);
ylim([-90 90]);
h1=gca;
set(h1,'fontsize',15);

subplot(212)
%plot(freq,sig_1,'k','linewidth',1.5);
hold on
plot(T0pg.freq,T0pg.sig_2,'k','linewidth',1.5);
plot(T20pg.freq,T20pg.sig_2,'r','linewidth',1.5);
plot(T30pg.freq,T30pg.sig_2,'Color',[0.1 0.7 0.15],'linewidth',1.5);
plot(T40pg.freq,T40pg.sig_2,'b','linewidth',1.5);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$\sigma_{\theta,2}~(^{\circ})$','interpreter','latex')
xlim([0 1.5]);
ylim([0 45]);
h1=gca;
set(h1,'fontsize',15);
% h2 = legend('Stero','pgera','pgar','$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
% set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
sname = [Tinfo.figfolder,'theta_sig_pg'];
print(sname,'-dpng')
% %%
% %%
% fsmall        = [1:0.01:3];
% fsmall4       = 10^-2.5*fsmall.^-4;
% figure('units','inches','position',[1 1 10 8],'color','w');
% subplot(211)
% semilogy(T0pg.freq,T0pg.Sf,'k','linewidth',1.5);
% hold on
% semilogy(T20pg.freq,T20pg.Sf,'r','linewidth',1.5);
% semilogy(T30pg.freq,T30pg.Sf,'b','linewidth',1.5);
% semilogy(T40pg.freq,T40pg.Sf,'Color',[0.1 0.7 0.15],'linewidth',1.5);
% % semilogy(fsmall,fsmall4,'b','linewidth',2);
% grid
% xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
% ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
% ylim([10^-4.5 10^-1.5])
% xlim([0 3]);
% h1=gca;
% set(h1,'fontsize',15);
% h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
% set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
% 
% 
% subplot(212)
% semilogy(T0pg.dir,T0pg.Sd,'k','linewidth',1.5);
% hold on
% semilogy(T20pg.dir,T20pg.Sd,'r','linewidth',1.5);
% semilogy(T30pg.dir,T30pg.Sd,'b','linewidth',1.5);
% semilogy(T40pg.dir,T40pg.Sd,'Color',[0.1 0.7 0.15],'linewidth',1.5);
% grid
% xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
% ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
% xlim([0 360]);
% h1=gca;
% set(h1,'fontsize',15);
% sname = [Tinfo.figfolder,'Sf_Sd'];
% print(sname,'-dpng')
% 
% %%
% figure('units','inches','position',[1 1 10 8],'color','w');
% subplot(211)
% %plot(freq,th_1,'k','linewidth',1.5);
% hold on
% plot(freq,T0pg.th,'k','linewidth',1.5);
% plot(freq,T20pg.th,'r','linewidth',1.5);
% plot(freq,T30pg.th,'b','linewidth',1.5);
% plot(freq,T40pg.th,'Color',[0.1 0.7 0.15],'linewidth',1.5);
% grid
% xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
% ylabel('$\theta~(^{\circ})$','interpreter','latex')
% %ylim([10^-4.5 10^-1.5])
% xlim([0 3]);
% ylim([-90 90]);
% h1=gca;
% set(h1,'fontsize',15);
% 
% subplot(212)
% %plot(freq,sig_1,'k','linewidth',1.5);
% hold on
% plot(freq,T0pg.sig,'k','linewidth',1.5);
% plot(freq,T20pg.sig,'r','linewidth',1.5);
% plot(freq,T30pg.sig,'b','linewidth',1.5);
% plot(freq,T40pg.sig,'Color',[0.1 0.7 0.15],'linewidth',1.5);
% grid
% xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
% ylabel('$\sigma_{\theta}~(^{\circ})$','interpreter','latex')
% xlim([0 3]);
% ylim([0 45]);
% h1=gca;
% set(h1,'fontsize',15);
% h2 = legend('$\sigma_{\theta} = 0^{\circ}$','$\sigma_{\theta} = 20^{\circ}$','$\sigma_{\theta} = 30^{\circ}$','$\sigma_{\theta} = 40^{\circ}$')
% set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
% sname = [Tinfo.figfolder,'theta_sig'];
% print(sname,'-dpng')
% 
% %%
% [val, ifmax] = min(abs(freq-1.5));
% 
% fsmall        = [1:0.01:3];
% fsmall4       = 10^-2.5*fsmall.^-4;
% figure('units','inches','position',[1 1 10 8],'color','w');
% subplot(211)
% semilogy(T0cam.freq,T0cam.Sf,'k','linewidth',2);
% hold on
% % semilogy(freq,T20cam.Sf,'r','linewidth',1.5);
% % semilogy(freq,T30cam.Sf,'Color',[0.1 0.7 0.15],'linewidth',1.5);
% semilogy(T40cam.freq,T40cam.Sf,'r','linewidth',2);
% semilogy(T0pg.freq,T0pg.Sf,'k','linewidth',2,'linestyle',':');
% semilogy(T40pg.freq,T40pg.Sf,'r','linewidth',2,'linestyle',':');
% % semilogy(fsmall,fsmall4,'b','linewidth',2);
% grid
% xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
% ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
% ylim([10^-4.5 10^-1.5])
% xlim([0 3]);
% h1=gca;
% set(h1,'fontsize',15);
% h2 = legend('stereo, $\sigma_{\theta} = 0^{\circ}$','stereo, $\sigma_{\theta} = 40^{\circ}$','press, $\sigma_{\theta} = 0^{\circ}$','press, $\sigma_{\theta} = 40^{\circ}$')
% set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
% 
% 
% subplot(212)
% semilogy(T0cam.dir,T0cam.Sd,'k','linewidth',2);
% hold on
% semilogy(T40cam.dir,T40cam.Sd,'r','linewidth',2);
% semilogy(T0pg.dir,T0pg.Sd,'k','linewidth',2,'linestyle',':');
% semilogy(T40pg.dir,T40pg.Sd,'r','linewidth',2,'linestyle',':');
% grid
% xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
% ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
% xlim([0 360]);
% h1=gca;
% set(h1,'fontsize',15);
% sname = [Tinfo.figfolder,'Sf_Sd'];
% print(sname,'-dpng')