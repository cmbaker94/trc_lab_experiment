function estimate_spread(Tinfo,rsinst,rsiter,pginst,pgiter,frange,matorpy)
% estimate directional spectra properties

addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\codes\insitu'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rinst = 14;
% riter = 200;

% frange = [0.3 0.8];%[0.25 1.25];
dirrange = [0 359];

Tinfo = trial_files(Tinfo);

% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

if Tinfo.filt == 1
    subname = '_filt';
else
    subname = '';
end

%% STEP 2: Define figure folders

% STEP 1: Load data

if matorpy ==0
    geninfo = 'regx28-33_regy-10-4_dx50_dy50';
    Sftheta = csvread([Tinfo.savefolder,'spec_cam_',geninfo,'_grid_rand_',num2str(rsinst),'inst_',num2str(rsiter),'iter',subname,'.csv']);
    freq = csvread([Tinfo.savefolder,'freq_cam_',geninfo,'_grid_rand_',num2str(rsinst),'inst_',num2str(rsiter),'iter',subname,'.csv'])';
    dir = csvread([Tinfo.savefolder,'dirs_cam_',geninfo,'_grid_rand_',num2str(rsinst),'inst_',num2str(rsiter),'iter',subname,'.csv'])';
elseif matorpy == 1
    geninfo = 'regx28_33_regyneg10_4_dx50_dy50';
    S = load([Tinfo.savefolder,'dirspec_cam_',geninfo,'_',num2str(rsinst),'inst_',num2str(rsiter),'iter',subname,'.mat']);
    % below is temp
    iS = mean(abs(S.Spec),[1,2]);
    Sft = abs(S.Spec(:,:,squeeze(iS(1,1,:))>10^(-6)));
    Sftheta = mean(Sft,3);
    display(['size ',num2str(size(Sft,3))])
%     Sftheta = abs(S.Savg);
    freq = S.freqs;
    dir = S.dirs;
end
cam = calc_dir_prop(Sftheta,freq,dir,frange,dirrange);
[sprd_cam,S_cam,dir_cam] = fit_cos2s(Tinfo,cam.dir,cam.Sd)
% sprd_cam
% dir_cam
clear Sftheta freq dir

%% lidar

% geninfo = 'regx28-33_regy-10-4_dx50_dy50';
if matorpy == 0
    geninfo = 'regx28-33_regy-10-4_dx50_dy50';
    Sftheta = csvread([Tinfo.savefolder,'spec_lid_',geninfo,'_grid_rand_',num2str(rsinst),'inst_',num2str(rsiter),'iter',subname,'.csv']);
    freq = csvread([Tinfo.savefolder,'freq_lid_',geninfo,'_grid_rand_',num2str(rsinst),'inst_',num2str(rsiter),'iter',subname,'.csv'])';
    dir = csvread([Tinfo.savefolder,'dirs_lid_',geninfo,'_grid_rand_',num2str(rsinst),'inst_',num2str(rsiter),'iter',subname,'.csv'])';
elseif matorpy == 1
    geninfo = 'regx28_33_regyneg10_4_dx50_dy50';
    S = load([Tinfo.savefolder,'dirspec_lid_',geninfo,'_',num2str(rsinst),'inst_',num2str(rsiter),'iter',subname,'.mat']);
    % below is temp
    iS = mean(abs(S.Spec),[1,2]);
    Sft = abs(S.Spec(:,:,squeeze(iS(1,1,:))>10^(-6)));
    Sftheta = mean(Sft,3);
    display(['size ',num2str(size(Sft,3))])
%     Sftheta = abs(S.Savg);
    freq = S.freqs;
    dir = S.dirs;
end
lid = calc_dir_prop(Sftheta,freq,dir,frange,dirrange);
[sprd_lid,S_lid,dir_lid] = fit_cos2s(Tinfo,lid.dir,lid.Sd);
% sprd_cam
% dir_lid
clear Sftheta freq dir

%%
if matorpy == 0
    Sftheta = csvread([Tinfo.savefolder,'spec_press_press_sz_',num2str(pginst),'inst_',num2str(pgiter),'iter.csv']);
    freq = csvread([Tinfo.savefolder,'freq_press_press_sz_',num2str(pginst),'inst_',num2str(pgiter),'iter.csv'])';
    dir = csvread([Tinfo.savefolder,'dirs_press_press_sz_',num2str(pginst),'inst_',num2str(pgiter),'iter.csv'])';
elseif matorpy == 1
    S = load([Tinfo.savefolder,'dirspec_press_sz_rand_',num2str(pginst),'inst_',num2str(pgiter),'iter_hyd.mat']);
    % below is temp
    iS = mean(abs(S.Spec),[1,2]);
    Sft = abs(S.Spec(:,:,squeeze(iS(1,1,:))>10^(-6)));
    Sftheta = mean(Sft,3);
    display(['size ',num2str(size(Sft,3))])
    display('using hydrostatic version')
%     Sftheta = abs(S.Savg);
    freq = S.freqs;
    dir = S.dirs;
end
pg = calc_dir_prop(Sftheta,freq,dir,frange,dirrange);
clear Sftheta
[sprd_pg,S_pg,dir_pg] = fit_cos2s(Tinfo,pg.dir,pg.Sd);
% sprd_pg
% dir_pg

%%

% puvfile     = ['E:\data\processed\insitu\',Tinfo.sz.tdate,'\PUV_p06.csv'];
% puv06 = read_puv_csv(puvfile);
% puvfile     = ['E:\data\processed\insitu\',Tinfo.sz.tdate,'\PUV_p11.csv'];
% puv11 = read_puv_csv(puvfile);
% freqrange = [0.2 1.25];
% spec = T0puv06.th;
% [emean] = calc_energy_weighted_mean(spec,freq,freqrange)

if matorpy == 0
    psname = [Tinfo.savefolder,'sprd_cos2sdist',subname,'.mat'];
elseif matorpy == 1
    psname = [Tinfo.savefolder,'sprd_cos2sdist',subname,'_matver.mat'];
end
eval(['save -v7.3 ',psname,' sprd_*',' S_*',' dir_*']);

end
