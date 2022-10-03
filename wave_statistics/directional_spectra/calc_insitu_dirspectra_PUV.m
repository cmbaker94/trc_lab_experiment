% This code calculates the directional spectra using in situ pressure and
% velocities with the PUV function.

clear all
close all
clc

%% STEP 0: Input trial information

% Trial info
Tinfo.Hs        = 0.25;
Tinfo.Tp        = 2;
Tinfo.tide      = 1.07;
Tinfo.spread    = 40;
time_range      = [5 35];
Tinfo.insituA   = 1; % 1 = surf zone, 2 = inner shelf

% spectral analysis details
units = 1/9807; % convert Pa to m freshwater
sampfreq = 100; % Hz
window_length = 180; 
merge = 5; 
plt = 0;

% fill in the rest of the trial information based on wave conditions
% FIX THIS IN THE FUTURE
Tinfo = trial_files(Tinfo);

% general path and names
Tinfo.cam = TRC_camera_info(Tinfo.cam);

% Data and figure storage
[Tinfo] = wc_comp_store(Tinfo);

% STEP 1: Load Stereo Reconstruction Data

[camera.time] = cam_time(Tinfo);

%% Step 1: Load data

datapath = 'E:/data/processed/insitu/';
if Tinfo.insituA == 1
    insitudata = Tinfo.sz.tdate;
elseif Tinfo.insituA == 2
    insitudata = Tinfo.is.tdate;
end

pstruct = load([datapath,insitudata,'/',insitudata,'-insitu.mat']); % 40 deg sprd, sz array

%% STEP 2: compute directional spectra

% get file pressure gauge names
pnames  = fieldnames(pstruct.press.xyzd);
it      = sampfreq*60*time_range(1):sampfreq*60*time_range(end)-1;
var = {'freq', 'SSE', 'k', 'PP', 'PU', 'PV', 'UU', 'VV', 'UV', ...
        'dir1', 'spread1', 'dir2', 'spread2', 'spread2alt', 'dir2_swell', 'spread2_swell', ...
        'Hsig_swell', 'Hsig_ig', 'centroid_swell', 'centroid_ig', 'in', 'out', ...
        'PsqP', 'PsqU', 'Depth', 'correction', 'a1', 'b1', 'a2', 'b2'};
% compute bulk statistics
for i = 1:length(pnames)
    % load pressure data
    instnum(i) = str2double(pnames{i}(6:end));
    eval(['p = pstruct.press.press',pnames{i}(6:end),'(it)*units;'])
    eval(['u = pstruct.vel.u',pnames{i}(6:end),'(it);'])
    eval(['v = pstruct.vel.v',pnames{i}(6:end),'(it);'])
    eval(['xyz = pstruct.press.xyzd.press',pnames{i}(6:end),'(1,1:3);']) 
    
    zsens = mean(p);
    Depth = zsens + .05; % assume sensor is about 5cm above the bed
    p = p - mean(p);

    [ freq, SSE, k, PP, PU, PV, UU, VV, UV, ...
        dir1, spread1, dir2, spread2, spread2alt, dir2_swell, spread2_swell, ...
        Hsig_swell, Hsig_ig, centroid_swell, centroid_ig, in, out, ...
        PsqP, PsqU, Depth, correction, ...
        a1, b1, a2, b2]  =  ...
        my_PUVspectra( p*100, u*100, v*100 , Depth*100, ...
        window_length, merge, plt );
    
    for j = 1:length(var)
        eval(['p',num2str(instnum(i),'%02d'),'.',var{j},'=',var{j},';'])
    end
    
%     figure
    subplot(3,1,1)
    semilogy(freq,SSE,'m-','LineWidth',2); hold on
    xlim([0 1.1])
    ylabel('spec')
    grid on
    subplot(3,1,2)
    plot(freq,dir2,'m','LineWidth',2); hold on
    xlim([0 1.1])
    ylim([-20 20])
    ylabel('dir')
    grid on
    subplot(3,1,3)
    plot(freq,spread2*180/pi,'m','LineWidth',2); hold on
    %plot(freq,spread2alt*180/pi,'y','LineWidth',2); hold on
    xlim([0 1.1])
%     ylim([0 40])
    ylabel('spread')
    grid on
    xlabel('freq (hz)')
end
% clear vel



%% STEP 3: Save data

% for i = 1:length(pnames)
%     eval(['savenames(i) = "p',pnames{i}(6:end),'";'])
% end

sname = [datapath,insitudata,'/',insitudata,'-PUV.mat'];
eval(['save -v7.3 ',sname,' p01',' p02',' p03',' p04',' p05',' p06',' p07',' p08',' p09',' p10',' p11',' p12'])%' eta'

for i = 1:length(instnum)
    eval(['M(:,1) = p',num2str(instnum(i),'%02d'),'.freq;'])
    eval(['M(:,2) = p',num2str(instnum(i),'%02d'),'.SSE;'])
    eval(['M(:,3) = p',num2str(instnum(i),'%02d'),'.dir2;'])
    eval(['M(:,4) = p',num2str(instnum(i),'%02d'),'.spread2;'])
    fname = [datapath,insitudata,'/PUV_p',num2str(instnum(i),'%02d'),'.csv'];
    writematrix(M,fname)
end