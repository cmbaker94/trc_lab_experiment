% processing code prior to stereo reconstructions
clear all
close all
clc

trial = 'TRM-08-29-2018-2349UTC';

% copy files after 
frames = 00000:07199;

filestart = ['H:\TRC_Fall_Experiment\',trial,'\',trial,'_Scene1_JPEG\'];%_JPEG
fileend =  ['H:\TRC_Fall_Experiment\',trial,'\',trial,'_Scene1_JPEG_15min\'];

subfolder = ['frames_',num2str(frames(1),'%05.f'),'-',num2str(frames(end),'%05.f'),'\'];
eval(['!mkdir ',fileend])


cam = 2;
for i = 1:length(cam)
    camtemp = ['c',num2str(cam(i))];
%     eval(['!mkdir ',fileend,subfolder,'images\',camtemp])
    for j = 1:length(frames) 
        file2copy = [filestart,'\Movie1_Scene1_',camtemp,'_',num2str(frames(j),'%05.f')];
         file2rename = ['Movie1_Scene1_',camtemp,'_',num2str(frames(j),'%05.f')];
%         file2copy = [filestart,'Movie1_Scene1_',camtemp,'_',num2str(frames(j),'%05.f')];
        eval(['!copy ',file2copy,'* ',fileend])
        
        fullrename = dir([fileend,file2rename,'*.jpg']);
        eval(['cd ',fileend])
        eval(['!ren ',fullrename.name,' ',file2rename,'.jpg'])
    end
end
