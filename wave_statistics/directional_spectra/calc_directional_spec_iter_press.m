%This code will read the gridded camera data and then use that infromation
%to determine the two-dimensional frequency directional spectrum

%% STEP 1: Clear All
close all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\codes\insitu'))
addpath(genpath('E:\codes\trc_lab_experiment\toolbox')) 

%% Input trial details

Tinfo.spread = 30;
Tinfo.Hs = 0.25;
Tinfo.Tp = 2;
Tinfo.tide = 1.07;
Tinfo.filt = 1;
instnum = 5;
iternum = 200;
inHz = 100;

%% SM defaults

SM.dirs       = [0:1:359];
SM.freqs      = [0.01:0.01:3];
SM.S          = zeros(300,360);
EP.nfft       = inHz*(2^5); %  32 sec window
EP.dres       = 360;
EP.method     = 'IMLM';
IDdef.datatypes  = {'elev' 'elev' 'elev' 'elev' 'elev'};
IDdef.fs         = inHz;
% IDdef.depth      = 0.1664;

%% Get trail details

Tinfo = trial_files(Tinfo);

% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);
[camera.time] = cam_time(Tinfo);

dx = 0.5;
dy = 0.5;
xlid    = [27.5:dx:33];
ylid    = [-10:dy:4];

%% STEP 2: Load Bathymetry
A        = load('E:/data/processed/lidar/Riegl/TRC_bathymetry_estimate_line.mat');
% dy       = dy;
% yp       = ylid;
% h        = Tinfo.tide-A.h;
X = A.xp;
% [X,Y]    = meshgrid(A.xp,yp);
% X        = X';
% Y        = Y';
% H        = repmat(h,[1 length(yp)]);
% [~,ihmin] = min(abs(xlid(1)-X));
% [~,ihmax] = min(abs(xlid(2)-X));
% H = h(ihmin:ihmax);
% IDdef.depth = mean(H);

%% Load lidar data

[press,wg] = load_insitu(Tinfo);
IDdef.depth = mean(Tinfo.tide-press.xyz(3,:));
z        = press.eta-mean(press.eta,1);
XL       = press.xyz(1,:);
YL       = press.xyz(2,:);

%% STEP 4: Select Points

% for it = 1:300
%     ix = randi([1, size(XL,2)],instnum,1);
%     xd = XL(ix);
%     yd = YL(ix);
%     testx(it) = nanstd(xd);
%     testy(it) = nanstd(yd);
% end
count = 1;

for it = 1:600
    ix = randperm(size(XL,2),instnum);
    xd = XL(ix);
    yd = YL(ix);
    stdmet = 1;
    if std(yd)<3
        display('not spread enough');
        stdmet = 0;
    end
    if stdmet == 1
        ID              = IDdef;
        ID.data       = z(:,ix); % 6000 x 14
        N             = size(xd,2);
        ID.layout     = [xd;yd; zeros(1,N)]; % 3 x 14
        [Smout,EPout] = dirspec(ID,SM,EP);
        
        Spec(:,:,count) = Smout.S;
        
        count = count + 1;
        count
        if count == iternum+1
            break
        end
    end
end


Savg = mean(Spec,3);
dirs = Smout.dirs;
freqs = Smout.freqs;
if Tinfo.filt == 1
    specsave = [Tinfo.savefolder,'dirspec_press_regx',num2str(round(xlid(1))),'_',num2str(round(xlid(end))),'_regyneg',num2str(abs(round(ylid(1)))),'_',num2str(round(ylid(end))),'_dx',num2str(dx*100),'_dy',num2str(dy*100),'_filt_inst',num2str(instnum),'_iter',num2str(iternum)];
end
eval(['save -v7.3 ',specsave,' dirs',' freqs',' Savg',' Spec'])

%% Try DIWASP


dtheta        = abs(dirs(2)-dirs(1));
df            = abs(freqs(2)-freqs(1));
Sf            = nansum(Savg,2)*dtheta;
Sd            = nansum(Savg,1)*df;
fsmall        = [1:0.01:3];
fsmall4       = 10^-2*fsmall.^-4;

figure('units','inches','position',[1 1 10 8],'color','w');
subplot(211)
semilogy(Smout.freqs,Sf,'k','linewidth',3);
hold on
semilogy(fsmall,fsmall4,'b','linewidth',2);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$S_{f}~(\mathrm{m^{2}/Hz})$','interpreter','latex')
ylim([10^-4.5 10^-1.5])
xlim([0 3]);
h1=gca;
set(h1,'fontsize',15);

subplot(212)
semilogy(Smout.dirs,Sd,'k','linewidth',3);
hold on
grid
xlabel('$\mathrm{Dir.}~\mathrm{(^{\circ})}$','interpreter','latex')
ylabel('$S_{d}~(\mathrm{m^{2}/degree})$','interpreter','latex')
xlim([0 360]);
h1=gca;
set(h1,'fontsize',15);



