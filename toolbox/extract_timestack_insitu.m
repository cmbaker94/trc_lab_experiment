function [inst] = extract_timestack_insitu(time,camera,press,lidar,inst)
% extract a chunk of the timeseries for plotting
% INPUT:
% time = time necessary to extract
% camera = camera data
% press = insitu data
% lidar = lidar data
% instname = name of instruments, {onshore, offshore}, string format,
%            example: instname = {'6','11'};
% OUTPUT:
% inst: structure timeseries at the offshore and onshore instrament

% sensorloc = {'on','off'};
g = 9.81;
rho = 1000;

% Identify time to extract
% in situ
[temp,iI(1)] = nanmin(abs(time(1)-press.time));
[temp,iI(2)] = nanmin(abs(time(end)-press.time));
inst.It          = press.time(iI(1):iI(2));
inst.Ito         = (inst.It-camera.time(1))*24*3600;

% camera
[temp,iC(1)] = nanmin(abs(time(1)-camera.time));
[temp,iC(2)] = nanmin(abs(time(end)-camera.time));
inst.Ct = camera.time(iC(1):iC(2));
inst.Cto = (inst.Ct-camera.time(1))*24*3600;
% lidar
[temp,iL(1)] = nanmin(abs(time(1)-lidar.time));
[temp,iL(2)] = nanmin(abs(time(end)-lidar.time));
inst.Lt = lidar.time(iL(1):iL(2));
inst.Lto = (inst.Lt-camera.time(1))*24*3600;


for i = 1:length(inst.name)
    instnum     = find(contains(press.name,join(['press',inst.name(i)],'')));
    x           = press.xyz(1,instnum);
    y           = press.xyz(2,instnum);
    if round(x)<29
        sensorloc = 'off';
    elseif round(x)>29
        sensorloc = 'on';
    end
    h.press     = nanmean(press.press(:,instnum));
    h.camera    = nanmean(camera.z(:,instnum));
    h.lidar     = nanmean(lidar.z(:,instnum));
    loc.Iz      = (press.press(iI(1):iI(2),instnum)-h.press)/(g*rho);
    loc.Ifz     = press.eta(iI(1):iI(2),instnum);
    loc.Cz      = camera.z(iC(1):iC(2),instnum)-h.camera;
    loc.Lz      = lidar.z(iL(1):iL(2),instnum)-h.lidar;
    eval(['inst.',sensorloc,' = loc;']);
    eval(['inst.',sensorloc,'.h = h;'])
    eval(['inst.',sensorloc,'.x = x;'])
    eval(['inst.',sensorloc,'.y = y;'])
end

