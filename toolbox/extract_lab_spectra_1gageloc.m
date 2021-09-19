function [spec] = extract_lab_spectra_1gageloc(insitu,cam,lidar,gage,grab_tseries);
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

% wave gage
if gage{1}(1:2) == 'wg'
    wg2find = gage{1};
    idtemp = find(contains(insitu.inst,wg2find));
    for i = 1:length(idtemp)
        tempcheck = insitu.inst(idtemp(i));
        if tempcheck{1} == gage{1}
            id = idtemp(1);
            break
        end
    end
    spec.wg.xyz = insitu.xyz(id,:);
    spec.wg.Hs = insitu.Hs(:,id);
    spec.wg.Tp = insitu.Tp(:,id);
    spec.wg.See = insitu.See(id,:)';
    spec.wg.f = insitu.f(id,:)';
    spec.wg.Seec = squeeze(insitu.Seec(id,:,:));
    idx = round(insitu.xyz(id,1));
    ind = round(insitu.xyz(:,1))==idx;
    spec.wg.See_yavg = nanmean(insitu.See(ind,:),1);
    spec.wg.See_ystd = nanstd(insitu.See(ind,:),1);
    spec.wg.Hs_yavg = nanmean(insitu.Hs(ind));
    spec.wg.Hs_ystd = nanstd(insitu.Hs(ind));
    clear wg
    
elseif gage{1}(1:2) == 'pr'
    % pressure gage
    press2find = gage{1};
    idtemp = find(contains(insitu.inst,press2find));
    for i = 1:length(idtemp)
        tempcheck = insitu.inst(idtemp(i));
        if tempcheck{1} == gage{1}
            id = idtemp(i);
            break
        end
    end
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
    
    % Cameras
    [temp,ix] = nanmin(abs(cam.x-spec.press.xyz(1,1)));
    [temp,iy] = nanmin(abs(cam.y-spec.press.xyz(1,2)));
    spec.cam.x = cam.x(ix);
    spec.cam.y = cam.y(iy);
    if grab_tseries == 1
        spec.cam.z  = squeeze(cam.z(iy,ix,:));
        spec.cam.eta = squeeze(cam.eta(iy,ix,:));
    end
    spec.cam.See = squeeze(cam.See(iy,ix,:));
    spec.cam.See_yavg = squeeze(nanmean(cam.See(:,ix,:),1));
    spec.cam.See_ystd = squeeze(nanstd(cam.See(:,ix,:),1));
    spec.cam.f = squeeze(cam.f(iy,ix,:));
    if isnan(mean(spec.cam.f))
        spec.cam.f = squeeze(nanmean(cam.f(:,ix,:),1));
    end
    spec.cam.Seec = squeeze(cam.Seec(iy,ix,:,:));
    spec.cam.Hs = squeeze(cam.Hs(iy,ix));
    spec.cam.Hs_yavg = squeeze(nanmean(cam.Hs(:,ix),1));
    spec.cam.Hs_ystd = squeeze(nanstd(cam.Hs(:,ix),1));
    spec.cam.Tp = squeeze(cam.Tp(iy,ix));
    
    % Lidar
    [temp,ix] = nanmin(abs(lidar.x-spec.press.xyz(1,1)));
    [temp,iy] = nanmin(abs(lidar.y-spec.press.xyz(1,2)));
    spec.lidar.x = lidar.x(ix);
    spec.lidar.y = lidar.y(iy);
    if grab_tseries == 1
        spec.lidar.z  = squeeze(lidar.z(iy,ix,:));
        spec.lidar.eta = squeeze(lidar.eta(iy,ix,:));
    end
    spec.lidar.See = squeeze(lidar.See(iy,ix,:));
    spec.lidar.See_yavg = squeeze(nanmean(lidar.See(:,ix,:),1));
    spec.lidar.See_ystd = squeeze(nanstd(lidar.See(:,ix,:),1));
    spec.lidar.f = squeeze(lidar.f(iy,ix,:));
    if isnan(mean(spec.lidar.f))
        spec.lidar.f = squeeze(nanmean(lidar.f(:,ix,:),1));
    end
    spec.lidar.Seec = squeeze(lidar.Seec(iy,ix,:,:));
    spec.lidar.Hs = squeeze(lidar.Hs(iy,ix));
    spec.lidar.Hs_yavg = squeeze(nanmean(lidar.Hs(:,ix),1));
    spec.lidar.Hs_ystd = squeeze(nanstd(lidar.Hs(:,ix),1));
    spec.lidar.Tp = squeeze(lidar.Tp(iy,ix));
end
end