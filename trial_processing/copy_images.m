function copy_images(trial,frames)

% filestart = ['D:\TRC_Fall_Experiment\',trial,'\',trial,'_Scene1_JPEG\'];
filestart = ['H:\TRC_Fall_Experiment\',trial,'\',trial,'_Scene1_JPEG\'];%_JPEG
fileend =  ['F:\Metashape\',trial,'_Scene1\'];



subfolder = ['frames_',num2str(frames(1),'%05.f'),'-',num2str(frames(end),'%05.f'),'\'];
eval(['!mkdir ',fileend])
eval(['!mkdir ',fileend,subfolder])
eval(['!mkdir ',fileend,subfolder,'images\'])
eval(['!mkdir ',fileend,subfolder,'dems\'])
eval(['!mkdir ',fileend,subfolder,'pointclouds\'])

cam = 1:3;
for i = 1:length(cam)
    camtemp = ['c',num2str(cam(i))];
    eval(['!mkdir ',fileend,subfolder,'images\',camtemp])
    for j = 1:length(frames) 
%         file2copy = [filestart,camtemp,'\Movie1_Scene1_',camtemp,'_',num2str(frames(j),'%05.f')];
        file2copy = [filestart,'Movie1_Scene1_',camtemp,'_',num2str(frames(j),'%05.f')];
        eval(['!copy ',file2copy,'* ',fileend,subfolder,'images\',camtemp,'\'])
    end
end
end
