function [x,y,eta] = grab_press_sensor(press,sensor)
% select sensor and extract timeseries and location

instnum     = find(contains(press.name,join(['press',sensor],'')));
x           = press.xyz(1,instnum);
y           = press.xyz(2,instnum);
eta         = press.eta(:,instnum);