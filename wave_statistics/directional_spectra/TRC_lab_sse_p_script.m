clear all
close all
clc

datapath = '/Users/cmbaker9/Documents/Research/Lab_Experiments/data/processed/insitu/';
run = '08-30-2018-2129UTC';
%pstruct = load('09-01-2018-2213UTC-insitu.mat'); %0 deg sprd
%pstruct = load('09-06-2018-1841UTC-insitu.mat'); %40 deg sprd, is array
pstruct = load([datapath,run,'/',run,'-insitu.mat']); % 40 deg sprd, sz array

units = 1/9807; % convert Pa to m freshwater
sampfreq = 100; % Hz

p = pstruct.press.press6*units;
p = p(sampfreq*60*5:sampfreq*60*35-1); % clip to middle 35 mins
zsens = mean(p);
p = p - mean(p);

plt = 1;

u = pstruct.vel.u6;
v = pstruct.vel.v6;

u = u(sampfreq*60*5:sampfreq*60*35-1); % clip to middle 35 mins
v = v(sampfreq*60*5:sampfreq*60*35-1); % clip to middle 35 mins

Depth = zsens + .05; % assume sensor is about 5cm above the bed
window_length = 180; 
merge = 5; 

[ freq, SSE, k, PP, PU, PV, UU, VV, UV, ...
    dir1, spread1, dir2, spread2, spread2alt, dir2_swell, spread2_swell, ...
    Hsig_swell, Hsig_ig, centroid_swell, centroid_ig, in, out, ...
    PsqP, PsqU, Depth, correction, ...
    a1, b1, a2, b2]  =  ...
    my_PUVspectra( p*100, u*100, v*100 , Depth*100, ...
    window_length, merge, plt );


% %figure
% subplot(2,1,1)
% semilogy(freq,SSE,'r-','LineWidth',2); hold on
% xlim([0 1.1])
% subplot(2,1,2)
% errorbar(freq,dir2,spread2alt*180/pi/2,'r','LineWidth',2); hold on
% xlim([0 1.1])

figure
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
ylim([0 30])
ylabel('spread')
grid on
xlabel('freq (hz)')
