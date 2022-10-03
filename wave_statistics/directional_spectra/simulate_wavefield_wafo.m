% Test directional spectra method with WAFO toolbox

% create synthetic spectra
close all
clear all
clc
addpath(genpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS'))
addpath(genpath('E:\code\wafo'))


%% Define variables

Hmo     = 0.25;
Tp      = 2;
s       = 170;% maximum: 170

% Define simulation domain
dx  = 0.2;
dt  = 1/8;
L   = 30;
T   = 60;

% calc method for  remote sensing in lab  
instnum = 14;
iternum = 120;
inHz = dt;

%% Create spectra

sdata   = [Hmo Tp];
wc      = 33/Tp;%(1/Tp)*2*pi;
w       = linspace(0,wc,257);
S       = jonswap(w,sdata);

data    = [s, nan, wc];
D       = spreading(linspace(-pi,pi,101),'cos2s',0,data,S.w,0);
W       = mkdspec(S,D);
S.note  = ['Demospec: ',S.note,', truncated at 2*wp'];

%% Perform Simulation

Nx  = round(L/dx);
Ny  = Nx+1;
Nt  = round(T/dt);
x   = linspace(0,L-dx,Nx); % m
y   = linspace(0,L,Ny); % m
t   = linspace(0,T-dt,Nt); % s
dy  = dx;

% Call WAFO codes to simulate surface trace
% if fftdim = 2, 2D ifft at each time.  if fftdim = 1, 1D ifft at each space
fftdim =1;
[Y,~] = seasim(W,Nx,Ny,Nt,dx,dy,dt,fftdim,1);


% [Y,S_Measured]=simsurfCB(W);

figure
for i = 1:10%length(Y.t)
    pcolor(Y.x,Y.y,Y.Z(:,:,i))
    shading interp
    colorbar
    colormap(cmocean('balance'));
    caxis([-0.25 0.25])
    pause(0.5)
end

figure
pcolor(W.w/(2*pi),W.theta*(180/pi),W.S)
shading interp
colorbar

figure
plot(W.w/(2*pi),nansum(W.S,1))

figure
plot(W.theta*(180/pi),nansum(W.S,2))
