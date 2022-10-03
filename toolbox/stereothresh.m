function [Ioz] = stereothresh(Zim,thresh,minarea)
% Function to smooth and threshold  stereo 


Zim(Zim>0.75) = NaN;
% Zim = smoothdata(Zim,1,'movmedian',[3 3],'omitnan');%[4 4]
% Zim = smoothdata(Zim,2,'movmedian',[1 1],'omitnan');% [1 1]

% THIS  ONE
Zim = smoothdata(Zim,1,'movmean',[6 6],'omitnan');%[5 5]%  just this one!
Zim = smoothdata(Zim,2,'movmean',[1 1],'omitnan');% [3 3]

% Zim = movmax(Zim,4,1,'omitnan');
% Zim = movmax(Zim,4,2,'omitnan');

Iz   = 0*ones(size(Zim));
Iz(Zim>thresh) = 1;
% Ioz = bwareaopen(Iz,minarea,8);
Ioz = Iz;

