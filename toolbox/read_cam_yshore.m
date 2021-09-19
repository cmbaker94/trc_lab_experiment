function cam = read_cam_yshore(Tinfo,xshore,xavg,tavg)
% extract alongshore transects
% tavg = number of bins to average over in time
% xavg = cross-shore average of stored transects

load([Tinfo.cam.datafolder,'dem_transect_x',num2str(xshore),'m_xavg',num2str(xavg*100),'cm.mat']);
cam.time    = Tinfo.cam.timevec;

eval(['z',num2str(xshore),'(y',num2str(xshore),'(:,1)>13,:,:) = NaN;'])
eval(['z',num2str(xshore),'(y',num2str(xshore),'(:,1)<-13,:,:) = NaN;'])
eval(['cam.z = z',num2str(xshore),';'])
eval(['cam.y = y',num2str(xshore),';'])
cam.stillwat = nanmean(cam.z,2);

if tavg>1
    for j = 1:length(cam.z(1,:))-tavg
        ztemp(:,j) = nanmean(cam.z(:,j:j+tavg-1),2);
    end
    cam.z = ztemp;
end

