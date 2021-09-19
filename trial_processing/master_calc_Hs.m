function master_calc_Hs(Tinfo)
% Compute sea-surface elevation spectra, Hs, and Tp for in situ, stereo,
% and lidar

% Set up paths and clear workspace
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\insitu'))
addpath(genpath('E:\code\cameras'))
addpath(genpath('E:\code\trc_lab_experiment\toolbox'))

%%%%%%%%%%%%%%%%%%%%% USER MANIPULATED SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% spectral analysis details
sdet.WL              = 2^5;
sdet.OL              = 2^4;
sdet.frange     = [0.25, 1.2];
sdet.nancutoff  = 0.1;

%% STEP 1: Create paths, files and naming

% general path and names
datapath    = 'E:\';

Tinfo = trial_files(Tinfo);

% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

% % Bathymetry
% addpath([datapath,'codes\analytic'])
% % load trc bathy
bathymetry = load('E:\data\processed\lidar\Riegl\TRC_bathymetry_estimate_line.mat');
% SWAN0 = load('\Users\cmbaker9\Documents\Research\Lab_Experiments\codes\analytic\SWAN_output_0_degree_spread.mat');
% SWAN40 = load('\Users\cmbaker9\Documents\Research\Lab_Experiments\codes\analytic\SWAN_output_40_degree_spread.mat');

%% STEP 2: Stereo Reconstruction + VeloLidar: compute spectra
sdet.plotspec = 0;

% stereo
sdet.Hz     = 8;
cam         = calc_Hs_stereo(Tinfo,sdet);
if Tinfo.filt == 1
    sname       = 'stereo_Hs_spectra_filt';
else
    sname       = 'stereo_Hs_spectra';
end
eval(['save -v7.3 ',Tinfo.savefolder,sname,' cam Tinfo sdet'])

%%
% lidar
sdet.Hz     = 10;
lidar       = calc_Hs_lidar(Tinfo,sdet);
if Tinfo.filt == 1
    sname       = 'velodyne_Hs_spectra_filt';
else
    sname       = 'velodyne_Hs_spectra';
end
eval(['save -v7.3 ',Tinfo.savefolder,sname,' lidar Tinfo sdet'])

%% Insitu

sdet.Hz         = 100;
sdet.maxfac     = 1.2;

time            = [Tinfo.cam.timevec(1) Tinfo.cam.timevec(end)+datenum(0,0,0,0,0,1/Tinfo.cam.Hz)]-datenum(Tinfo.cam.tstart(1:end-3),'mm-dd-yyyy-HHMM');

szvsis          = 0; % 0 for sz array
[sz.Hs,sz.Tp,sz.See,sz.Seec,sz.f,sz.inst,sz.xyz,sz.Tinfo,sz.data] = calc_Hs_spectra_insitu(Tinfo,time,sdet.WL,sdet.OL,sdet.frange,sdet.maxfac,bathymetry,szvsis); % ADD BACK IN TIME

%%
szvsis          = 1; % 1 for is array and want to grab time range from sz array
[is.Hs,is.Tp,is.See,is.Seec,is.f,is.inst,is.xyz,is.Tinfo,sz.data] = calc_Hs_spectra_insitu(Tinfo,time,sdet.WL,sdet.OL,sdet.frange,sdet.maxfac,bathymetry,szvsis);

% time            = [Tinfo.cam.timevec(1) Tinfo.cam.timevec(end)+datenum(0,0,0,0,0,1/Tinfo.cam.Hz)];
% 
% szvsis          = 0; % 0 for sz array
% [sz.Hs,sz.Tp,sz.See,sz.Seec,sz.f,sz.inst,sz.xyz,sz.Tinfo,sz.data] = calc_Hs_spectra_insitu(Tinfo.sz.tdate,time,sdet.WL,sdet.OL,sdet.frange,sdet.maxfac,bathymetry,szvsis,[]); % ADD BACK IN TIME
% 
% %%
% szvsis          = 1; % 1 for is array and want to grab time range from sz array
% [is.Hs,is.Tp,is.See,is.Seec,is.f,is.inst,is.xyz,is.Tinfo,sz.data] = calc_Hs_spectra_insitu(Tinfo.is.tdate,time,sdet.WL,sdet.OL,sdet.frange,sdet.maxfac,bathymetry,szvsis,Tinfo.sz.tdate);

% create summary table
sz.sumtable = table(sz.inst',sz.xyz,sz.Hs',sz.Tp');
is.sumtable = table(is.inst',is.xyz,is.Hs',is.Tp');
sz.sumtable.Properties.VariableNames = {'Instrument' 'xyz' 'Hs' 'Tp'};
is.sumtable.Properties.VariableNames = {'Instrument' 'xyz' 'Hs' 'Tp'};

% sz.Hs(find(sz.Hs<0.05)) = NaN;
% is.Hs(find(is.Hs<0.05)) = NaN;

sname       = 'insitu_Hs_spectra';
eval(['save -v7.3 ',Tinfo.savefolder,sname,' sz',' is',' sdet'])

end

