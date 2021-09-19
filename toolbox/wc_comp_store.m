function [Tinfo] = wc_comp_store(Tinfo)
% Create the file naming structure and folders based on the wave conditions
% INPUT: Tinfo
% OUTPUT: file structures

Tinfo.comp = ['Hs',num2str(Tinfo.Hs*100),'_Tp',num2str(Tinfo.Tp),'_tide',num2str(Tinfo.tide*100),'_spread',num2str(Tinfo.spread)];
datarange = [datestr(Tinfo.cam.timevec(1)); datestr(Tinfo.cam.timevec(end))];
Tinfo.datarange = ['time_', datarange(1,13:14), datarange(1,16:17), '-' , datarange(2,13:14), datarange(2,16:17)];
Tinfo.procpath = [Tinfo.datapath,'data\processed\conditions\',Tinfo.comp];

eval(['mkdir ',Tinfo.datapath,'data\processed\conditions\'])
eval(['mkdir ',Tinfo.procpath])
eval(['mkdir ',Tinfo.procpath,'\',Tinfo.cam.tstart])
eval(['mkdir ',Tinfo.procpath,'\',Tinfo.cam.tstart,'\',Tinfo.datarange])
Tinfo.savefolder = [Tinfo.procpath,'\',Tinfo.cam.tstart,'\',Tinfo.datarange,'\'];

% figure folder
fssubfolder = datestr(date,'yy-mm-dd');
Tinfo.figfolder   = [Tinfo.datapath,'figures\meas_comp\',Tinfo.cam.trialname,'\',Tinfo.cam.trimname,fssubfolder,'\'];

% make figure folders
eval(['mkdir ',Tinfo.datapath,'figures\meas_comp\',Tinfo.cam.trialname]);
eval(['mkdir ',Tinfo.datapath,'figures\meas_comp\',Tinfo.cam.trialname,'\',Tinfo.cam.trimname]);
eval(['mkdir ',Tinfo.figfolder])