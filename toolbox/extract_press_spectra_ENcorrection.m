function spec = extract_press_spectra_ENcorrection(insitu,gages)
% Extract the spectra for insitu, lidar, and camera
% INPUT:
% insitu= processes insitu gage statistics
% cam   = processes stereo statistics
% lidar = processes lidar statistics
% gages = names of insitu gages to process data from (wg, press)
% grab_tseries = grab timeseries, if 0 = NO or 1 = YES
% OUTPUT:
% spec  = structure with the following data:
%   wg    = spectral parameters for wg
%   press = spectral parameters for pressure gages
%   cam   = spectral parameters for stereo
%   lidar = spectral parameters for lidar, velodyne

% pressure gage
press2find = gages{2}(6:end);
id = find(contains(insitu.inst,['eta',press2find]));
spec.press.xyz = insitu.xyz(id,:);
spec.press.Hs = insitu.Hs(:,id);
spec.press.Tp = insitu.Tp(:,id);
spec.press.See = insitu.See(id,:)';
spec.press.f = insitu.f(id,:)';
spec.press.Seec = squeeze(insitu.Seec(id,:,:));
idx = round(insitu.xyz(id,1));
ind = round(insitu.xyz(:,1))==idx;
spec.press.See_yavg = nanmean(insitu.See(ind,:),1);
spec.press.See_ystd = nanstd(insitu.See(ind,:),1);
spec.press.Hs_yavg = nanmean(insitu.Hs(ind));
spec.press.Hs_ystd = nanstd(insitu.Hs(ind));

end