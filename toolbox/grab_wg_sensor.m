function [x,y,eta] = grab_wg_sensor(wg,sensor)
% select sensor and extract timeseries and location

instnum     = find(contains(wg.names,join(['wg',sensor],'')));
x           = wg.xyz(1,instnum);
y           = wg.xyz(2,instnum);
eta         = wg.z(:,instnum);