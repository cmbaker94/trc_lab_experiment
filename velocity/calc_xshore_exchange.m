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

path = [datapath,'data/processed/insitu/'];
f1 = dir([path,'09-06*']);
for i = 1:length(f1)
    t1{i} = f1(i).name;
end
f2 = dir([path,'09-07*']);
for i = 1:length(f2)
    t2{i} = f2(i).name;
end
f3 = dir([path,'09-08*']);
for i = 1:length(f3)
    t3{i} = f3(i).name;
end
f5 = dir([path,'09-10*']);
for i = 1:length(f5)
    t5{i} = f5(i).name;
end

trials = [t1 t2 t3 t5];



%% Loop through two trials

for si = 1:length(trials)
    
%     Tinfo.spread = sprd(si);
%     Tinfo.Hs = 0.3;
%     Tinfo.Tp = 2;
%     Tinfo.tide = 1.07;
    
%     Tinfo = trial_files(Tinfo);
    
%     % general path and names
%     Tinfo.cam = TRC_camera_info(Tinfo.cam);
%     
%     % Data and figure storage
%     [Tinfo] = wc_comp_store(Tinfo);
    
    % STEP 1: Load Data
    
    % load insitu data
            F1          = matfile([datapath,'data/processed/insitu/',char(trials(1,si)),'/',char(trials(1,si)),'-insitu.mat']);
            arrayname = 'is';
    
    % press       = F1.press;
    % wg          = F1.wg;
    INinfo    = F1.Tinfo;
    cond = INinfo.cond;
    condstore{si} = cond;
    temp = strfind(cond,'H=');
    hs = str2double(cond(temp+2:temp+5));
    
    temp = strfind(cond,'T=');
    tp = cond(temp+2:temp+4);
    if isnan(str2double(tp(3)))
        tp = str2double(INinfo.cond(temp+2));
    else
        tp = str2double(tp);
    end
    
    temp = strfind(cond,'h=');
    h = str2double(cond(temp+2:temp+5));
    
    temp = strfind(cond,'s=');
    s = cond(temp+2:temp+3);
    
    if isempty(s)
        s = NaN;
    else
        if isnan(str2double(s(2)))
            s = str2double(INinfo.cond(temp+2));
        else
            s = str2double(s);
        end
    end
    condextract(si,:) = [hs tp h s];
%     pause(1)
    if isnan(s)
        dir = 0;
    elseif s == 4
        dir = 40;
    elseif s == 16
        dir = 20;
    elseif s == 65
        dir = 10;
    elseif s == 7
        dir = 30;
    end
    
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

    
    velnames  = fieldnames(vel.xyzd);
    
    mvavgt=50;%12 sec
    
    % compute bulk statistics
    for i = 1:length(velnames)
        instnum = str2double(velnames{i}(4:end));
        eval(['u(:,',num2str(instnum),') = vel.u',num2str(instnum),'(istart:iend);'])
%         ulf(:,instnum) = movmean(u(:,instnum),100*mvavgt)-nanmean(u(:,instnum));
%         uLF2(:,instnum) = pl64ta(ulf(:,instnum),mvavgt*30);
        uLF(:,instnum) = pl64ta(u(:,instnum),mvavgt*100)-pl64ta(u(:,instnum),(1/0.0031)*100);
%         uLF(:,instnum) = pl64ta(u(:,instnum),mvavgt*100);
        instname{instnum}             = velnames{i};
        eval(['xyz(instnum,:) = vel.xyzd.',velnames{i},'(1,1:3);'])
    end
    
%     figure;
%     plot(u(1:30000,3),'k')
%     hold on
%      plot(ulf(1:30000,3),'r')
%       plot(uLF2(1:30000,3),'g')
%        plot(uLF(1:30000,3),'b')
    
    uLFoff = uLF-nanmean(uLF);
%     uLFoff = uLF;
    uLFoff(uLFoff>0) = NaN;
    yloc = xyz(4:12,2);
       
    Uex = nanmean(uLFoff(4:12,:),'all');
    Uoff = [];
    uLFt = uLF-nanmean(uLF);
    for ti = 1:length(uLFt(:,1))
        uLF1 = uLFt(ti,4:12);
        uLF2 = uLF1(uLF1<0);
        y2 = yloc(uLF1<0);
        if length(y2)==1
            Uoff(ti) = uLF2;
        elseif length(y2)==0
            Uoff(ti) = NaN;
        else
            Uoff(ti) = trapz(y2,uLF2)/(max(y2)-min(y2));
        end
    end
    Uex2 = nanmean(Uoff);
%     display(Uex)
    
%     eval(['T',num2str(Tinfo.spread),'.uLF = uLF;'])
%     eval(['T',num2str(Tinfo.spread),'.Uex = Uex;'])
    uoff(:,:,si) = uLF;
    Uexchange(si) = Uex;
    Uexchange2(si) = Uex2;
    watlev(si) = h;
    Hs(si) = hs;
    Tp(si) = tp;
    sprd(si) = dir;
    clear vel u v Su* Sv* power
    clear Tinfo
end

%%
% f = f(:,1);
Uexchange = Uexchange2;

figure
subplot(2,2,1)
scatter(sprd,Uexchange)
hold on
Usprd(1) = median(Uexchange(sprd==0));
Usprd(2) = median(Uexchange(sprd==10));
Usprd(3) = median(Uexchange(sprd==20));
Usprd(4) = median(Uexchange(sprd==30));
Usprd(5) = median(Uexchange(sprd==40));
scatter([0 10 20 30 40],Usprd,'fill')

subplot(2,2,2)
scatter(Hs,Uexchange)
hold on
UHs(1) = median(Uexchange(Hs==0.2));
UHs(2) = median(Uexchange(Hs==0.25));
UHs(3) = median(Uexchange(Hs==0.3));
UHs(4) = median(Uexchange(Hs==0.35));
scatter([0.20 0.25 0.3 0.35],UHs,'fill')

subplot(2,2,3)
scatter(Tp,Uexchange)
hold on
UTp(1) = median(Uexchange(Tp==1.5));
UTp(2) = median(Uexchange(Tp==2));
UTp(3) = median(Uexchange(Tp==2.5));
UTp(4) = median(Uexchange(Tp==3));
scatter([1.5 2 2.5 3],UTp,'fill')

subplot(2,2,4)
scatter(watlev,Uexchange)
hold on
Uh(1) = median(Uexchange(watlev==1.0));
Uh(2) = median(Uexchange(watlev==1.07));
scatter([1.0 1.07],Uh,'fill')

%%
Uexchange = -Uexchange2;

figure('units','inches','position',[1 1 12 8],'color','w')
subplot(2,2,1)
scatter(Hs,Uexchange,60,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',1.5)
hold on
UHs(1) = nanmean(Uexchange(Hs==0.2));
UHs(2) = nanmean(Uexchange(Hs==0.25));
UHs(3) = nanmean(Uexchange(Hs==0.3));
UHs(4) = nanmean(Uexchange(Hs==0.35));
% scatter([0.20 0.25 0.3 0.35],UHs,'fill')
Hsplt = [0.20 0.25 0.3 0.35];
c = polyfit(Hsplt,UHs,1);
U_est = polyval(c,Hsplt);
plot(Hsplt,U_est,'Color',[1 0 0 0.5],'LineWidth',1.5)
scatter(Hsplt,UHs,50,'fill','r','MarkerEdgeColor','k')
h1=gca;
box on
grid on
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',18);
xlim([0.19 0.36])
ylim([0 0.03])
xlabel('$H_s~($m)','interpreter','latex','fontsize',20)
ylabel('$U_{ex}~($m/s)','interpreter','latex','fontsize',20)
text(0.195, 0.027,'(a)','interpreter','latex','fontsize',20);

subplot(2,2,2)
scatter(Tp,Uexchange,60,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',1.5)
hold on
UTp(1) = nanmean(Uexchange(Tp==1.5));
UTp(2) = nanmean(Uexchange(Tp==2));
UTp(3) = nanmean(Uexchange(Tp==2.5));
UTp(4) = nanmean(Uexchange(Tp==3));
Tpplt = [1.5 2 2.5 3];
c = polyfit(Tpplt,UTp,1);
U_est = polyval(c,Tpplt);
plot(Tpplt,U_est,'Color',[1 0 0 0.5],'LineWidth',1.5)
scatter(Tpplt,UTp,50,'fill','r','MarkerEdgeColor','k')
h1=gca;
box on
grid on
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',18);
xlim([1.4 3.1])
ylim([0 0.03])
xlabel('$T_p~($s)','interpreter','latex','fontsize',20)
ylabel('$U_{ex}~($m/s)','interpreter','latex','fontsize',20)
text(1.45, 0.027,'(b)','interpreter','latex','fontsize',20);

subplot(2,2,3)
scatter(sprd,Uexchange,60,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',1.5)
hold on
Usprd(1) = nanmean(Uexchange(sprd==0));
Usprd(2) = nanmean(Uexchange(sprd==10));
Usprd(3) = nanmean(Uexchange(sprd==20));
Usprd(4) = nanmean(Uexchange(sprd==30));
Usprd(5) = nanmean(Uexchange(sprd==40));
sprdplt = [0 10 20 30 40];
c = polyfit(sprdplt,Usprd,1);
U_est = polyval(c,sprdplt);
plot(sprdplt,U_est,'Color',[1 0 0 0.5],'LineWidth',1.5)
scatter(sprdplt,Usprd,50,'fill','r','MarkerEdgeColor','k')
h1=gca;
box on
grid on
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',18);
xlim([-2 42])
ylim([0 0.03])
xlabel('$\sigma_{\theta}~(^{\circ})$','interpreter','latex','fontsize',20)
ylabel('$U_{ex}~($m/s)','interpreter','latex','fontsize',20)
% set(h1,'yticklabel',[],'xticklabel',[]);
text(-1, 0.027,'(c)','interpreter','latex','fontsize',20);
% title(ax3,'$\sigma_{\theta}=40^{\circ}$','interpreter','latex','fontsize',28);
% errorbar(T40.specon.wg.f(67),10^(-1.3),10^-1.85,'Color','k','LineWidth',2.5)%wg.See(3)-wg.Seec(3,1))

subplot(2,2,4)
scatter(watlev,Uexchange,60,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',1.5)
hold on
Uh(1) = nanmean(Uexchange(watlev==1.0));
Uh(2) = nanmean(Uexchange(watlev==1.07));
% scatter([1.0 1.07],Uh,'fill')
hplt = [1.0 1.07];
c = polyfit(hplt,Uh,1);
U_est = polyval(c,hplt);
plot(hplt,U_est,'Color',[1 0 0 0.5],'LineWidth',1.5)
scatter(hplt,Uh,50,'fill','r','MarkerEdgeColor','k')
h1=gca;
box on
grid on
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',18);
xlim([0.995 1.075])
ylim([0 0.03])
xlabel('$h~($m)','interpreter','latex','fontsize',20)
ylabel('$U_{ex}~($m/s)','interpreter','latex','fontsize',20)
text(0.997, 0.027,'(d)','interpreter','latex','fontsize',20);
Sname1 = [figfolder,'exchange'];
print(Sname1,'-dpng')

%%
Uex25 = Uexchange(Hs == 0.25);
Hs25 = Hs(Hs == 0.25);
sprd25 = sprd(Hs == 0.25);
watlev25 = watlev(Hs == 0.25);
Tp25 = Tp(Hs == 0.25);


figure
subplot(2,2,1)
scatter(sprd25,Uex25,50,'MarkerEdgeColor',[0.6 0.6 0.6],'LineWidth',1.5)
hold on
Usprd(1) = nanmean(Uex25(sprd25==0));
Usprd(2) = nanmean(Uex25(sprd25==10));
Usprd(3) = nanmean(Uex25(sprd25==20));
Usprd(4) = nanmean(Uex25(sprd25==30));
Usprd(5) = nanmean(Uex25(sprd25==40));
scatter([0 10 20 30 40],Usprd,40,'fill','r')
h1=gca;
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',20);
xlim([-5 45])
% ylim([10^-5 10^-1])
% set(h1,'yticklabel',[],'xticklabel',[]);
% text(ax3,0.05, 10^(-1.3),'(c)','interpreter','latex','fontsize',24);
% title(ax3,'$\sigma_{\theta}=40^{\circ}$','interpreter','latex','fontsize',28);
% errorbar(T40.specon.wg.f(67),10^(-1.3),10^-1.85,'Color','k','LineWidth',2.5)%wg.See(3)-wg.Seec(3,1))
h2 = legend('0','10','20','30','40','interpreter','latex','fontsize',22)
set(h2,'orientation','vertical')

subplot(2,2,2)
scatter(Hs25,Uex25)
hold on
UHs(1) = nanmean(Uex25(Hs25==0.2));
UHs(2) = nanmean(Uex25(Hs25==0.25));
UHs(3) = nanmean(Uex25(Hs25==0.3));
UHs(4) = nanmean(Uex25(Hs25==0.35));
scatter([0.20 0.25 0.3 0.35],UHs,'fill')

subplot(2,2,3)
scatter(Tp25,Uex25)
hold on
UTp(1) = nanmean(Uex25(Tp25==1.5));
UTp(2) = nanmean(Uex25(Tp25==2));
UTp(3) = nanmean(Uex25(Tp25==2.5));
UTp(4) = nanmean(Uex25(Tp25==3));
scatter([1.5 2 2.5 3],UTp,'fill')

subplot(2,2,4)
scatter(watlev25,Uex25)
hold on
Uh(1) = nanmean(Uex25(watlev25==1.0));
Uh(2) = nanmean(Uex25(watlev25==1.07));
scatter([1.0 1.07],Uh,'fill')

