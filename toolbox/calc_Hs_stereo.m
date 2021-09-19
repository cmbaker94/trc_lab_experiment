function [cam] = calc_Hs_lidar(Tinfo,sdet)
% This function will compute the sea-surface elevation spectra and compute
% the significant wave height from stereo measurments over a specified time
% range. 
% INPUT:
% Tcam          = Trial information about the camera setup
% sdet          = structure with defined parameters for computing the
%                 sea-surface elevation spectra and Hs
%       WL      = window length, default = []
%       OL      = overlap length, default = []
%       Hz      = measurement frequency, default = 8 Hz
%       frange  = frequency range Hs integrated over, 2x1 vector: 
%                 [start, end], default = [0.25, 1] Hz
%       nancutoff = percent of nan's allowable in timeseries when computing
%                   the spectra, default = 10%
%       plotspec = flag, plot each spectra as computed, where 1 = plot
%                  and 0 = no plots, default = 0 
% OUTPUT:
% lidar         = structure with spectra (m^2/Hz), Hs (m), frequency (Hz)


% set defaults if not specified in input
if ~isfield(sdet,'WL')
    sdet.WL = [];
    sdet.OL = [];
end
if ~isfield(sdet,'Hz')
    sdet.Hz = 8;
end
if ~isfield(sdet,'frange')
    sdet.frange = [0.25, 1.25];
end
if ~isfield(sdet,'nancutoff')
    sdet.nancutoff = 0.1;
end
if ~isfield(sdet,'plotspec')
    sdet.plotspec = 0;
end

% load file
if Tinfo.filt == 1
    load([Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_res',num2str((Tinfo.cam.dx)*100),'cm_filtered.mat']); 
else
    load([Tinfo.cam.datafolder,'dem_region_x',num2str(Tinfo.cam.regx(1)),'to',num2str(Tinfo.cam.regx(end)),'m_yneg',num2str(abs(Tinfo.cam.regy(1))),'to',num2str(Tinfo.cam.regy(end)),'m_res',num2str((Tinfo.cam.dx)*100),'cm.mat']);
end

% Stereo Bulk Stats
% Hs_4std = 4*nanstd(cam.z,[],3);
% zmean = nanmean(cam.z,3);

[cam.Hs,cam.Tp,cam.See,cam.f,cam.Seec,cam.nanratio] = calc_Hs_spectra_remotesensing(z,squeeze(x(1,:)),squeeze(y(:,1)),sdet.WL,sdet.OL,sdet.frange,sdet.Hz,sdet.nancutoff,sdet.plotspec);
cam.x = squeeze(x(1,:));
cam.y = squeeze(y(:,1))';
end