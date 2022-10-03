% Test directional spectra method with WAFO toolbox

% create synthetic spectra
close all
clear all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\wafo'))
addpath(genpath('E:\codes\trc_lab_experiment\toolbox')) 


%% Define variables

Hmo     = 0.25;
Tp      = 2;
s       = 16;% maximum: 170

% Define simulation domain
dx  = 0.05;
dt  = 1/8;
L   = 30;%14
T   = 60*10;

% calc method for  remote sensing in lab 
selrand  =  0;
if selrand ==  1
    instnum = 16;
    iternum = 50;
elseif selrand == 0
end

inHz = 1/dt;
testnum = 13;

%% Create spectra

% translate to wavelengths in shallow depth
Lsz = Tp*sqrt(0.5*9.81);
% use this as offshore wavelength
% Tp = sqrt((2*pi*Lsz)/9.81);

sdata   = [Hmo Tp];% 1.5];
wc      = 33/Tp;%(1/Tp)*2*pi;
w       = linspace(0,wc,257);
S       = jonswap(w,sdata);

data    = [s, nan, wc];
D       = spreading(linspace(-pi,pi,101),'cos2s',0,data,S.w,0);
W       = mkdspec(S,D);
S.note  = ['Demospec: ',S.note,', truncated at 2*wp'];

%% Perform Simulation

Nx  = round(L/dx);
Ny  = Nx+1;
Nt  = round(T/dt);
x   = linspace(0,L-dx,Nx); % m
y   = linspace(0,L,Ny); % m
t   = linspace(0,T-dt,Nt); % s
dy  = dx;

% Call WAFO codes to simulate surface trace
% if fftdim = 2, 2D ifft at each time.  if fftdim = 1, 1D ifft at each space
fftdim =1;
% [Y,~] = seasim(W,Nx,Ny,Nt,dx,dy,dt,fftdim,1);

dfwafo = W.w(2)-W.w(1);
figure
plot(W.theta*(180/pi),nansum(W.S,2)*dfwafo)

W.dir = W.theta*(180/pi);
W.f = W.w/(2*pi);

testnum = 9
specsave = ['E:\data\processed\dirspec_wafo\','simulated_Hmo',num2str(Hmo*100),'_Tp',num2str(floor(Tp)),'_s',num2str(s),'_wafo_output_test',num2str(testnum)];
% eval(['save -v7.3 ',specsave,' Y',' W',' D',' S',' Hmo',' Tp',' s'])
load([specsave,'.mat'])
testnum =13

%% SM defaults

SM.dirs       = [0:1:359];
SM.freqs      = [0.01:0.01:3];
SM.S          = zeros(300,360);
EP.nfft       = inHz*(2^5); %  32 sec window
EP.dres       = 360;
EP.method     = 'EMEP';
IDdef.datatypes  = {'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev' 'elev'};
IDdef.fs         = inHz;
% IDdef.depth      = 0.1664;

%% Get trail details

if selrand ==  1
    dx = 0.05;
    dy = 0.05;
    xL    = [0:dx:33-27.49];%[27.5:dx:33];
    yL    = [0:dx:4--10];
    xloc = size(Y.Z,2)/2+[-length(xL)/2 length(xL)/2-1];
    yloc = round(size(Y.Z,1)/2)+[-round(length(yL))/2 round(length(yL))/2-1];
    z = Y.Z(yloc(1):yloc(2),xloc(1):xloc(2),:);
elseif selrand == 0
    xL = Y.x';
    yL = Y.y';
    z = Y.Z(:,1:length(xL),:);
end
XL = repmat(xL,length(yL),1)';
YL = repmat(yL',1,length(xL))';
IDdef.depth = 100;


% 
z = permute(z,[2 1 3]);

if selrand == 0
    xd = [-1 0 0.3 0.7 1.5 5 10 15 25 0 0 0 0 0 0 0]'+2;
    yd = [8 8 8 8 8 8 8 8 8 5.5 6.7 8 10 13 15 25];  
    for il = 1:length(xd)
        [~,ix] = min(abs(XL(:,1)-xd(il)));
        [~,iy] = min(abs(YL(1,:)-yd(il)));
        zd(il,:)    = squeeze(z(ix,iy,:));
    end
end

figure('units','inches','position',[1 1 6 8],'color','w');
v = VideoWriter(['E:\data\processed\dirspec_wafo\','simulated_Hmo',num2str(Hmo*100),'_Tp',num2str(floor(Tp)),'_s',num2str(s),'_test',num2str(testnum),'_seasurfaceelev.avi']);
v.FrameRate=8; %?????
open(v)
for i = 1%:60*8%length(Y.t)
    pcolor(XL,YL,z(:,:,i))
    shading interp
    colorbar
    colormap(cmocean('balance'));
    hold on
    if selrand == 0
        scatter(xd,yd,'m','MarkerFaceColor','m','MarkerEdgeColor','k')
    end
    caxis([-0.25 0.25])
    xlabel('$x$ (m)','interpreter','latex','fontsize',15)
    ylabel('$y$ (m)','interpreter','latex','fontsize',15)
    text(.2+max(xL),.5+max(yL),'$\eta$ (m)','interpreter','latex','fontsize',15)
    text(0,0.5+max(yL),['Time: ',num2str(round(i/8)),' s'],'interpreter','latex','fontsize',15)
    axis equal
    ylim([min(yL) max(yL)])
    xlim([min(xL) max(xL)]);
    h1=gca;
    set(h1,'fontsize',15);
    writeVideo(v,getframe(gcf))
%     pause(0.2)
    Sname1 = ['E:\data\processed\dirspec_wafo\','simulated_Hmo',num2str(Hmo*100),'_Tp',num2str(floor(Tp)),'_s',num2str(s),'_test',num2str(testnum),'_instruments_',num2str(i)];
    print(Sname1,'-dpng')

    clf
end
close(v)

%% STEP 4: Select Points
threshold = 10;
count = 1;

fsmall        = [1:0.01:3];
fsmall4       = 10^-2*fsmall.^-4;
figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
semilogy(W.f,nansum(W.S,1)*diff(W.dir(1:2))/9,'r','linewidth',3)
hold on
semilogy(fsmall,fsmall4,'b','linewidth',2);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
ylim([10^-4.5 10^-1.2])
xlim([0 3]);
h1=gca;
set(h1,'fontsize',15);

subplot(212)
semilogy(W.dir,nansum(W.S,2)*diff(W.f(1:2))/8,'r','linewidth',3)
hold on
grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([-180 180]);
ylim([10^-8 10^-3])
h1=gca;
set(h1,'fontsize',15);


v = VideoWriter(['E:\data\processed\dirspec_wafo\','simulated_Hmo',num2str(Hmo*100),'_Tp',num2str(floor(Tp)),'_s',num2str(s),'_test',num2str(testnum),'.avi']);
v.FrameRate=8; %?????
open(v)

if selrand == 1
    itlen  = 10000;
elseif selrand == 0
    itlen = 1;
    nanmet =  1;
end

for it = 1:itlen
    if selrand == 1
        isel = randi(length(XL(:)),instnum,1);%length(XL(:)),instnum);
        for il = 1:instnum
            [~,ix(il)] = min(abs(XL(:,1)-XL(isel(il))));
            [~,iy(il)] = min(abs(YL(1,:)-YL(isel(il))));
        end
        
        %     ix = randi([1, size(XL,1)],instnum,1);
        %     iy = randi([1, size(YL,2)],instnum,1);
        xd = XL(ix,1);
        yd = YL(1,iy);
        stdmet = 1;
        linemet = 1;
        if length(unique(round(xd))) < 2
            display('all in one line')
            linemet = 0;
        end
        %     if std(xd)<1.4 || std(yd)<3.5
        %         display('not spread enough');
        %         stdmet = 0;
        %     end
        if stdmet == 1 && linemet == 1
            for i = 1:instnum
                zt1    = squeeze(z(ix(i),iy(i),:));
                T      = 1:1:length(zt1);
                j      = find(isnan(zt1));
                nanmet(i) = 1;
                if 100*length(j)/length(zt1)<threshold
                    zt2    = zt1;
                    T2     = T;
                    zt2(j) = [];
                    T2(j)  = [];
                    zt1    = interp1(T2,zt2,T);
                    zt1     = fillmissing(zt1,'linear');
                    zt1    = zt1-nanmean(zt1);
                else
                    display('Too many NaNs');
                    nanmet(i) = 0;
                end
                zd(i,:) = zt1;
                %             hd(i) = hL(ix(i),iy(i));
            end
        end
    end
        if min(nanmet) == 1
            ID              = IDdef;
%             ID.depth     = hd;
            ID.data       = zd'; % 6000 x 14
            N             = length(xd);
            ID.layout     = [xd';yd; zeros(N,1)']; % 3 x 14
            options = {'PLOTTYPE',0};
            [Smout,EPout]   = dirspec(ID,SM,EP,options);
            
            % Ouput in spectral poower density for the ith frequency and
            % jth direction (i.e., per unit Hz.deg)
            df              = diff(Smout.freqs(1:2));
            dtheta          = diff(Smout.dirs(1:2));
            Sout            = Smout.S;%sqrt(2*Smout.S*dtheta*df);
            Spec(:,:,count) = Smout.S;
            
            [~,ith]=nanmin(abs(Smout.dirs-179)); 
%             Sd              = nansum(abs(Sout),1)*df;
%             Sd              = [Sd(ith+1:end) flipud(Sd(1:ith))];
%             for fi = 1:size(Sout,1)
%             	Sf(fi)  = trapz(Smout.S(fi,:),Smout.dirs); 
%                 Sf1(fi)  = trapz(Sout(fi,:),Smout.dirs);
%             end
            Sf = nansum(Sout,2)*dtheta;
            
            for fi = 1:size(Sout,2)
%             	Sd(fi)  = trapz(Smout.freqs,Smout.S(:,fi)); 
                Sd(fi)  = trapz(Smout.freqs,Sout(:,fi));
            end
%             Sd = Sd';
%             Sd1 = Sd1';
%             Sd2 = nansum(Sout,1)*df;
            Sd              = [Sd(ith+1:end) flipud(Sd(1:ith))];
%             Sd1              = [Sd1(ith+1:end) flipud(Sd1(1:ith))];
%             Sd2              = [Sd2(ith+1:end) flipud(Sd2(1:ith))];
            
            subplot(211)
            semilogy(Smout.freqs,Sf,'m','linewidth',1);
%             hold on
%             semilogy(Smout.freqs,Sf1,'r','linewidth',1);
%             semilogy(Smout.freqs,Sf2,'b','linewidth',1);
            
            subplot(212)
            semilogy(Smout.dirs-180,Sd,'m','linewidth',1);
%             hold on
%             semilogy(Smout.dirs-180,Sd1,'r','linewidth',1);
%             semilogy(Smout.dirs-180,Sd2,'b','linewidth',1);
            
            count = count + 1;
            count
            writeVideo(v,getframe(gcf))
        end
    if selrand == 1
        if count == iternum+1
            break
        end
    end
end
close(v)

% testS = sqrt(2*Smout.S*dtheta*df);
Savg = mean(Spec,3);%sqrt(2*mean(Spec,3)*dtheta*df);
dirs = Smout.dirs;
freqs = Smout.freqs;

% figure
% [~,ith]=nanmin(abs(dirs-179));
% for i = 1:6
%     tempspec = Spec(:,:,i);
%     Sd            = nansum(abs(tempspec),1)*df;
%     Sd = [Sd(ith+1:end) flipud(Sd(1:ith))];
%     semilogy(theta,Sd,'linewidth',1);
%     hold on
%     pause(1)
% end

% theta=dirs;
% EfthetaSwift = Savg;
% [~,ith]=nanmin(abs(theta-179));
% % theta = theta-180
% Sd            = nansum(abs(EfthetaSwift),1)*df;
% Sd = [Sd(ith+1:end) flipud(Sd(1:ith))];
% semilogy(theta,Sd,'k','linewidth',3);

% theta=90-theta; 
% theta(theta<0)=theta(theta<0)+360;
% EfthetaSwift(EfthetaSwift==0)=NaN; 
% theta(theta > 180) = theta(theta > 180) - 360; 
% [theta,I]=unique(theta);
% EfthetaSwift=EfthetaSwift(:,I); 
% [theta, dsort] = sort(theta);
% EfthetaSwift = EfthetaSwift(:,dsort); 
% indNan = sum(isnan(EfthetaSwift),2)>0;
% Sd            = nansum(EfthetaSwift,1)*df;
% figure
% plot(theta-180,Sd*10,'k','linewidth',3);
% hold on
% plot(W.theta*(180/pi),nansum(W.S,2)*dfwafo)


specsave = ['E:\data\processed\dirspec_wafo\','simulated_Hmo',num2str(Hmo*100),'_Tp',num2str(floor(Tp)),'_s',num2str(s),'_Sfdestimates_test',num2str(testnum)];
eval(['save -v7.3 ',specsave,' dirs',' freqs',' Savg',' Spec',' W',' Hmo',' Tp',' s'])

%%
fact= 8;


dthetawafo    = (W.theta(2)-W.theta(1))*(180/pi); 
dfwafo        = (W.w(2)-W.w(1))/(2*pi);

dtheta        = abs(dirs(2)-dirs(1));
df            = abs(freqs(2)-freqs(1));
Sf            = nansum(Savg,2)*dtheta;
Sd            = nansum(Savg,1)*df;

[Spp,Spfreq] = pwelch(squeeze(z(5,5,:)),256,256/2,[],1/dt);
hs =  4*sqrt(trapz(Spfreq,Spp))
hs2 =  4*sqrt(trapz(Smout.freqs,abs(Sf)))
hs3 = 4*sqrt(trapz(W.w/(2*pi),nansum(W.S,1)*dthetawafo))

figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
tempspec = sum(Spec,2)*dtheta;
tempspec = squeeze(tempspec(:,1,:));
for i = 1:size(Spec,3)
    semilogy(freqs,tempspec(:,i),'y','linewidth',1);
    hold on
end
semilogy(freqs,Sf,'k','linewidth',3);
hold on
semilogy(W.w/(2*pi),nansum(W.S,1)*dthetawafo/fact,'r','linewidth',3)
% semilogy(W.w/(2*pi),(nansum(W.S,1)*dthetawafo),'g','linewidth',3)
% semilogy(W.w/(2*pi),((nansum(W.S,1)*dthetawafo).^2)/(2*dthetawafo*dfwafo),'m','linewidth',3)
semilogy(Spfreq,Spp,'b','linewidth',3)
semilogy(fsmall,fsmall4,'b','linewidth',2);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
ylim([10^-4.5 10^-1])
xlim([0 3]);
h1=gca;
set(h1,'fontsize',15);

[~,ith]=nanmin(abs(dirs-179));
Sd = [Sd(ith+1:end) flipud(Sd(1:ith))];

subplot(212)
tempspec = sum(Spec,1)*df;
tempspec = squeeze(tempspec(1,:,:))';
for i = 1:size(Spec,3)
    temps = [tempspec(ith+1:end,i); fliplr(tempspec(1:ith,i))];
    semilogy(dirs-180,temps,'y','linewidth',1);
    hold on
end
semilogy(dirs-180,Sd,'k','linewidth',3);
hold on
semilogy(W.theta*(180/pi),nansum(W.S,2)*dfwafo/fact,'r','linewidth',3)
hold on
grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([-180 180]);
ylim([10^-7 10^-3.8])
h1=gca;
set(h1,'fontsize',15);

Sname1 = ['E:\data\processed\dirspec_wafo\','simulated_Hmo',num2str(Hmo*100),'_Tp',num2str(floor(Tp)),'_s',num2str(s),'_test',num2str(testnum),'_semilog'];
print(Sname1,'-dpng')

figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
plot(Smout.freqs,Sf,'k','linewidth',3);
hold on
plot(W.w/(2*pi),nansum(W.S,1)*dthetawafo/fact,'r','linewidth',3)
semilogy(Spfreq,Spp,'b','linewidth',3)
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
ylim([0 0.03])
xlim([0 3]);
h1=gca;
set(h1,'fontsize',15);

subplot(212)
plot(Smout.dirs-180,Sd,'k','linewidth',3);
hold on
plot(W.theta*(180/pi),nansum(W.S,2)*dfwafo/fact,'r','linewidth',3)
hold on
grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([-180 180]);
ylim([0 0.0001])
h1=gca;
set(h1,'fontsize',15);

Sname1 = ['E:\data\processed\dirspec_wafo\','simulated_Hmo',num2str(Hmo*100),'_Tp',num2str(floor(Tp)),'_s',num2str(s),'_polarplot_test',num2str(testnum)];
print(Sname1,'-dpng')

%%

figure('units','inches','position',[1 1 17 5],'color','w');

conf = [0.9631 1.03156];

ax0 = axes('Position',[0.03 -2 0.28 0.96]);
hcb = pcolor(Savg)
caxis([10^(-8) 10^(-3.5)])
set(gca,'ColorScale','log')
colormap(jet)
cb = colorbar('Position',[0.035 0.095 0.022 0.65])
cb.Ruler.MinorTick = 'on';
h1=gca;
set(h1,'fontsize',18)

% ax05 = axes('Position',[0.03 0.9 0.9 0.2]);
% text(10^(-5),180,'(a) $\splot(T40cam.Sft(50,:))igma_{\theta}~\mathrm{(^{\circ})}$','interpreter','latex','color','k')

ax1 = axes('Position',[0.09 0.08 0.24 0.75]);
% Savg(:,size(Savg,2)+1) = NaN*size(Savg,1);
[h,c] = polarPcolor(freqs,[dirs],Savg,'ncolor',10,'colormap','jet','Ncircles',4,'Nspokes',9,'colBar',0);%,'RtickLabel',['0' '0.5' '1.0' '1.5' '2.0'])%,'typeRose','default')
% h1=gcf;
% set(h1,'fontsize',20)
% set(h1,'ZScale','log')
caxis([10^(-8) 10^(-3.5)])
% ylabel(c,' radial wind speed (m/s)');

ax2 = axes('Position',[0.39 0.08 0.24 0.75]);
sp = Savg;
sp(:,size(sp,2)+1) = NaN*size(sp,1);
fp = W.f(1):diff(W.f(1:2)):diff(W.f(1:2))*(size(sp,1)-2);
sp = zeros(length(W.dir),length(fp));
sp(:,1:length(W.f)) = W.S;
[h,c] = polarPcolor(fp,W.dir',sp/fact,'ncolor',10,'colormap','jet','Ncircles',4,'Nspokes',9,'colBar',0);%,'RtickLabel',['0' '1' '2' '3'])%,'typeRose','default')
% h1=gcf;
% set(h1,'fontsize',20)
% set(h1,'ZScale','log')
caxis([10^(-8) 10^(-3.5)])

ax3 = axes('Position',[0.72 0.15 0.26 0.83]);
% ADD THE CONFIDENCE INTERVAL
semilogy(dirs-180,Sd,'k','linewidth',3);
hold on
semilogy(W.theta*(180/pi),nansum(W.S,2)*dfwafo/fact,'r','linewidth',3)
% plot([190 190-180],[10^(-4)*conf(1) 10^(-6)*conf(2)],'k','LineWidth',2.5)
% plot([185 195],[10^(-4)*conf(2)+0.00005 10^(-6)*conf(2)+0.00005],'k','LineWidth',1)
% plot([185 195],[10^(-4)*conf(1)-0.00005 10^(-6)*conf(1)-0.00005],'k','LineWidth',1)
% grid
xlabel('$\theta~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{\eta\eta}(\theta)~(\mathrm{m^{2}/^{\circ}})$','interpreter','latex')
xlim([-180 180]);
ylim([10^(-8) 10^(-3.5)])
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'YScale','log')
set(h1,'fontsize',15);
h2 = legend('Estimated','Simulated')
set(h2,'interpreter','latex','fontsize',18,'orientation','vertical','Location','northeast');
text(-340-700,10^(-3.8),'(a) Estimated','interpreter','latex','color','k','fontsize',22)
text(-130-530,10^(-3.8),'(b) Simulated','interpreter','latex','color','k','fontsize',22)
text(95-270,10^(-3.8),'(c)','interpreter','latex','color','k','fontsize',22)
text(-410-750,10^(-4.5),'$S_{\eta\eta}(f,\theta)$','interpreter','latex','color','k','fontsize',22)
% text(-390,10^(-4.58),'$\left(\frac{\mathrm{m}^2}{^{\circ}f}\right)$','interpreter','latex','color','k','fontsize',16)
% text(-240,10^(-5.0),'$f$ (Hz)','interpreter','latex','color','w','fontsize',18)
% text(-30,10^(-5.0),'$f$ (Hz)','interpreter','latex','color','w','fontsize',18)
% text(-290,10^(-4.494),'$\theta=$','interpreter','latex','color','k','fontsize',18)
% text(-82,10^(-4.494),'$\theta=$','interpreter','latex','color','k','fontsize',18)


Sname1 = ['E:\data\processed\dirspec_wafo\','simulated_Sfd_Hmo',num2str(Hmo*100),'_Tp',num2str(floor(Tp)),'_s',num2str(s),'_test',num2str(testnum)];
print(Sname1,'-dpng')
