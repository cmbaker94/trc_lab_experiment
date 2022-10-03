% This code will extract and save time series from the camera dems and
% lidar mat file to compute the directional spectra as an 'array'.

% Set up paths and clear workspace
clear all
close all
clc
addpath(genpath('/Users/cmbaker9/Documents/MTOOLS'))
addpath(genpath('/Users/cmbaker9/Documents/Research/Lab_Experiments/codes/insitu'))
addpath(genpath('/Users/cmbaker9/Documents/Research/Lab_Experiments/codes/trc_lab_experiment/toolbox'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trial info
Tinfo.Hs = 0.25;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
Tinfo.spread = 0;

%%%%%%%%%%%%%%%%%% END OF USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%

%% STEP 0: Prepare file structures, folders, etc.

Tinfo = trial_files(Tinfo);

% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

%% STEP 1: Load insitu data

% [press,wg] = load_insitu(Tinfo);
% 
% sensors = {'11','6'};
% for i = 1:length(sensors)
%     sensor = sensors(i);
%     [x,y,eta] = grab_press_sensor(press,sensor);
% end

F1          = matfile([Tinfo.datapath,'data/processed/insitu/',Tinfo.sz.tdate,'/',Tinfo.sz.tdate,'-insitu.mat']);
press          = F1.press;
wg       = F1.wg;
tstart      = F1.Tinfo;
Hz          = press.Hz;

rho = 1000;
g = 9.81;
plotcompspec=1;

pressno = '6';
eval(['tempress = press.press',pressno,';'])
eval(['xy = press.xyzd.press',pressno,'(1:2);'])
time = 0:0.01:(length(tempress)-1)/100;

figure('units','inches','position',[1 1 15 5],'Color','w');
plot(time,tempress,'k','LineWidth',0.5)
title(['\textbf{Pressure Gage}: pg',pressno,', $x$ = ',num2str(round(xy(1))),' m, $y$ = ',num2str(round(xy(2))),' m'],'interpreter','latex','fontsize',20);
xlabel('$t$ (s)','interpreter','latex','fontsize',20);
ylabel('Pressure (N/m$^2$)','interpreter','latex','fontsize',20);
h1=gca;
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1*get(h1,'ticklength'));
set(h1,'fontsize',15);
xlim([0 600])
grid on

Sname1 = [Tinfo.figfolder,Tinfo.comp,'_timeseries_comparison_press',pressno];
print(Sname1,'-dpng')

clf
temz = (tempress/(rho*g))-nanmean((tempress(1:20)/(rho*g)));
plot(time,temz,'k','LineWidth',0.5)
title(['\textbf{Pressure Gage}: pg',pressno,', $x$ = ',num2str(round(xy(1))),' m, $y$ = ',num2str(round(xy(2))),' m'],'interpreter','latex','fontsize',20);
xlabel('$t$ (s)','interpreter','latex','fontsize',20);
ylabel('$z$ (m)','interpreter','latex','fontsize',20);
h1=gca;
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1*get(h1,'ticklength'));
set(h1,'fontsize',15);
xlim([0 600])
ylim([-0.3001 0.3001])
grid on

Sname1 = [Tinfo.figfolder,Tinfo.comp,'_timeseries_comparison_hyd',pressno];
print(Sname1,'-dpng')