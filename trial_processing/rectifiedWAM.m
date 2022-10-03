function rectifiedWAM(Tinfo,ixlim,iylim,idxdy,awin,ovr)
% function to create a wave-averaged movie. credits to B. Bruder for code
% INPUT:
% Tinfo             trial information
% ixlim, iylim      x,y limitis (m)
% idxdy             resolution (m)
% awin              averaging window (s)
% ovr               overlap window (s)
% OUTPUT:
% saved mat and video vile of wave averaged movie

%% STEP 1: Create paths, files and naming

Tinfo = trial_files(Tinfo);

% Stereo Reconstructions
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

%% Load File Names + Times
subname = '';
geoname = ['x',num2str(ixlim(1)),'to',num2str(round(ixlim(2))),'_y',num2str(round(iylim(2))),'_res',num2str(idxdy*100),'cm',subname];
odir = [Tinfo.savefolder(1:74),'orthos\',geoname];

cnames = ls([odir,'\c2_*.tiff']);
L = cnames(:,4:8);
% L = L(1:360,:);
% display('only selecting subsection')
t = str2num(L)/Tinfo.cam.Hz;
dt=mode(diff(t));
ts=(t-t(1));

%% Create Time Stamps and Windows

Tbin= (awin/2):ovr: (ts(end)-awin/2); % Centered On
Tlow=Tbin-awin/2; % Upper Limit
Thigh=Tbin+awin/2; % Lower Limit

%% Create Structure to Hold Images
Io=imread(fullfile(odir,['c2_',L(1,:),'.tiff']));
[r c co]=size(Io);
Iwam=zeros(r,c,length(Tbin));
Ncounter=zeros(1,length(Tbin));
clear('Io')


%% Load Images

for k=1:length(L)
    
% Load Image (Gray- Matlab cannot handle RGB array size)
I=double(rgb2gray(imread(fullfile(odir,['c2_',L(k,:),'.tiff']))));

% Find Bins Image Belongs into
tcheck=ts(k);    
lcheck=Tlow-tcheck;
hcheck=Thigh-tcheck;
ind=find(lcheck<=0 & hcheck>0);

% Add Image to Average
Iwam(:,:,ind)=Iwam(:,:,ind)+I;

% Add Counter for average at end
Ncounter(ind)=Ncounter(ind)+1;

k
end

disp('Loaded Images')

%% Take Average
for k=1:length(Tbin)
    Iwam(:,:,k)=Iwam(:,:,k)./Ncounter(k);
end
Iwam=uint8(Iwam);

TbinS=Tbin;
Tbin=Tbin+t(1);
disp('Averaged Images')

%% Save File (MAT)

oname=['c2_',L(1,:),'_',L(end,:),'_win',num2str(awin),'_ovr',num2str(ovr)];
sname = fullfile(Tinfo.savefolder(1:74),'WAMs',geoname,oname);
eval(['!mkdir ',sname])
save(fullfile(sname,[oname 'R.mat']),'Iwam','Ncounter','TbinS','awin','ovr','TbinS')
disp('saved Mat File')

%% Make A Movie
    
v = VideoWriter(fullfile(sname,[oname,'R.avi']));
v.FrameRate=1; %?????
open(v)

for k=1:length(Tbin)
writeVideo(v,Iwam(:,:,k))

k
end  
   close(v) 
   disp('saved movie')

