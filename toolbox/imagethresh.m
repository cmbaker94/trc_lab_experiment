function [Iob] = imagethresh(IM,thresh,minarea,idxmin,idxmax)
% Function to smooth and threshold  images 

% Irbw = smoothdata(IM,1,'movmedian',[3 3],'omitnan');
% Irbw = smoothdata(Irbw,2,'movmedian',[3 3],'omitnan'); %[1,1]

% THIS  ONE
Irbw = smoothdata(IM,1,'movmean',[6 6],'omitnan'); %[2  2]
Irbw = smoothdata(Irbw,2,'movmean',[1 1],'omitnan'); %[3,3]

% Irbw = movmax(IM,4,1,'omitnan');
% Irbw = movmax(Irbw,1,2,'omitnan');

Irbw = Irbw(:,idxmin:idxmax);
% apply the threshold to
Ib   = imbinarize(Irbw,thresh);
% Iob = bwareaopen(Ib,minarea,8);
Iob = Ib;
