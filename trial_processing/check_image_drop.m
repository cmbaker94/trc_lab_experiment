function check_image_drop(trial)

addpath('E:\code\cameras\reconstruction_processing')

% folder = ['D:\TRC_Fall_Experiment\',trial,'\',trial,'_Scene1_JPEG\'];
folder = ['H:\TRC_Fall_Experiment\',trial,'\',trial,'_Scene1_JPEG\'];

cam = 1:3;
for i = 1:length(cam)
    camtemp = ['c',num2str(cam(i))];
    jpgtemp = ['Movie1_Scene1_',camtemp,'_'];
%     fileList = dir(fullfile([folder,camtemp,'\'], [jpgtemp,'*']));
    fileList = dir(fullfile(folder, [jpgtemp,'*']));
    for j = 1:length(fileList) 
        tempname = fileList(j,1).name;
        timestreams = str2double(tempname(24:end-4));
        timeepoc = streams2Epoch(timestreams);
        timestamp(i,j) = timeepoc;
    end
    eval(['diff',num2str(i),' = diff(squeeze(timestamp(i,:))'');']);
end

scatter(1:length(diff1),diff1)
display(['max c1: ',num2str(max(diff1))])
display(['max c1: ',num2str(max(diff2))])
display(['max c1: ',num2str(max(diff3))])

display(['min c1: ',num2str(min(diff1))])
display(['min c1: ',num2str(min(diff2))])
display(['min c1: ',num2str(min(diff3))])

%%
tvec = 0:1/8:(length(diff2)/8)-1/8;

figure('units','inches','position',[1 1 10 5],'Color','w');
scatter(tvec,diff2,10,'k')
hold on
title(trial,'interpreter','latex','fontsize',20)
box on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
set(h1,'ticklength',1*get(h1,'ticklength'));
set(h1,'fontsize',15);
xlabel('$t$ (s)','interpreter','latex','fontsize',20);
ylabel('d$t$ (s)','interpreter','latex','fontsize',20);
% ylim([10^-4.1 10^-2])
xlim([tvec(1) tvec(end)])
sname = ['timestep_check_imagedrop'];
% print(['H:\TRC_Fall_Experiment\',trial,'\',sname],'-dpng')
print(['H:\TRC_Fall_Experiment\',trial,'\',sname],'-dpng')
end