function [x,y,h,z]=domain2Dlab(Tinfo)
% create the theoretical domain for the lab
% INPUT:
% Tinfo  = trial information, needed for tide
% OUTPUT:
% x - cross-shore
% y - alongshore
% h - empty facility
% z - filled still water facility

bathymetry = load('E:\data\processed\lidar\Riegl\TRC_bathymetry_estimate_line.mat');
x = bathymetry.xp;
y = [-14:diff(x(1:2)):14.01];
h = repmat(bathymetry.h',length(y),1);
ih(1) = find(round(y,2)==-13.30);
ih(2) = find(round(y,2)==13.30);
h(1:ih(1),:)= max(h,[],'all')+zeros(size(h(1:ih(1),:)));
h(ih(2):end,:) = max(h,[],'all')+zeros(size(h(ih(2):end,:)));
iz =  find(round(bathymetry.h,2)==Tinfo.tide);
z=h;
z(ih(1)+1:ih(2)-1,1:iz(1)) = Tinfo.tide+zeros(size(h(ih(1)+1:ih(2)-1,1:iz(1))));

[x,y]=meshgrid(x,y);