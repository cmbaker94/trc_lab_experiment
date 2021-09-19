function estimate_spread(Tinfo,rsinst,rsiter,pginst,pgiter,frange)
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

geninfo = 'regx28-33_regy-10-4_dx50_dy50';
Sftheta = csvread([Tinfo.savefolder,'spec_cam_',geninfo,'_grid_rand_',num2str(rsinst),'inst_',num2str(rsiter),'iter',subname,'.csv']);
freq = csvread([Tinfo.savefolder,'freq_cam_',geninfo,'_grid_rand_',num2str(rsinst),'inst_',num2str(rsiter),'iter',subname,'.csv'])';
dir = csvread([Tinfo.savefolder,'dirs_cam_',geninfo,'_grid_rand_',num2str(rsinst),'inst_',num2str(rsiter),'iter',subname,'.csv'])';
cam = calc_dir_prop(Sftheta,freq,dir,frange,dirrange);
[sprd_cam,S_cam,dir_cam] = fit_cos2s(cam.dir,cam.Sd)
sprd_cam
dir_cam
clear Sftheta freq dir

%% lidar

geninfo = 'regx28-33_regy-10-4_dx50_dy50';
Sftheta = csvread([Tinfo.savefolder,'spec_lid_',geninfo,'_grid_rand_',num2str(rsinst),'inst_',num2str(rsiter),'iter',subname,'.csv']);
freq = csvread([Tinfo.savefolder,'freq_lid_',geninfo,'_grid_rand_',num2str(rsinst),'inst_',num2str(rsiter),'iter',subname,'.csv'])';
dir = csvread([Tinfo.savefolder,'dirs_lid_',geninfo,'_grid_rand_',num2str(rsinst),'inst_',num2str(rsiter),'iter',subname,'.csv'])';
lid = calc_dir_prop(Sftheta,freq,dir,frange,dirrange);
clear Sftheta freq dir
[sprd_lid,S_lid,dir_lid] = fit_cos2s(lid.dir,lid.Sd);
sprd_cam
dir_lid

%%

Sftheta = csvread([Tinfo.savefolder,'spec_press_press_sz_',num2str(pginst),'inst_',num2str(pgiter),'iter.csv']);
freq = csvread([Tinfo.savefolder,'freq_press_press_sz_',num2str(pginst),'inst_',num2str(pgiter),'iter.csv'])';
dir = csvread([Tinfo.savefolder,'dirs_press_press_sz_',num2str(pginst),'inst_',num2str(pgiter),'iter.csv'])';
pg = calc_dir_prop(Sftheta,freq,dir,frange,dirrange);
clear Sftheta
[sprd_pg,S_pg,dir_pg] = fit_cos2s(pg.dir,pg.Sd);
sprd_pg
dir_pg

%%

% puvfile     = ['E:\data\processed\insitu\',Tinfo.sz.tdate,'\PUV_p06.csv'];
% puv06 = read_puv_csv(puvfile);
% puvfile     = ['E:\data\processed\insitu\',Tinfo.sz.tdate,'\PUV_p11.csv'];
% puv11 = read_puv_csv(puvfile);
% freqrange = [0.2 1.25];
% spec = T0puv06.th;
% [emean] = calc_energy_weighted_mean(spec,freq,freqrange)


psname = [Tinfo.savefolder,'sprd_cos2sdist',subname,'.mat'];
eval(['save -v7.3 ',psname,' sprd_*',' S_*',' dir_*']);

end
