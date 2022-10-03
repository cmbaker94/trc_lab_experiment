function [stereo] = combthresh_prep_stereo(F1,pltflag,ixsel)
% prepare sterio to for combined thresholding
if pltflag == 1
    stereo.orig = load(F1,'x','y','z');
    
    [~,idxmin] = min(abs(ixsel(1)-stereo.orig.x(1,:)));
    [~,idxmax] = min(abs(ixsel(2)-stereo.orig.x(1,:)));
    
    stereo.x = stereo.orig.x(:,idxmin:idxmax);
    stereo.y = stereo.orig.y(:,idxmin:idxmax);
    stereo.z = squeeze(stereo.orig.z(:,idxmin:idxmax,:));
    stereo.orig.msse = nanmean(stereo.orig.z,3);
else
    stereo = load(F1,'x','y','z');
    
    [~,idxmin] = min(abs(ixsel(1)-stereo.x(1,:)));
    [~,idxmax] = min(abs(ixsel(2)-stereo.x(1,:)));
    
    stereo.x = stereo.x(:,idxmin:idxmax);
    stereo.y = stereo.y(:,idxmin:idxmax);
    stereo.z = squeeze(stereo.z(:,idxmin:idxmax,:));
end

% mssetemp = nanmean(stereo.z,3);
% mssetemp(mssetemp<Tinfo.tide) = NaN;
mssetemp = movmean(stereo.z,5*8*2,3,'omitnan');% 10 waves
[iy(1)] = find(stereo.y(:,1)==-13.05);
[iy(2)] = find(stereo.y(:,1)==13.05);

for i = 1:size(mssetemp,3)
    mssetemp(iy(1):iy(2),:,i) = smoothdata(mssetemp(iy(1):iy(2),:,i),1,'movmean',[6 6],'omitnan');%[5 5]%  just this one!
    mssetemp(iy(1):iy(2),:,i) = smoothdata(mssetemp(iy(1):iy(2),:,i),2,'movmean',[6 6],'omitnan');% [3 3]
end
for i = 1:20
    mssetemp(:,:,i) = mssetemp(:,:,20);
    mssetemp(:,:,end-i+1) = mssetemp(:,:,end-20);
end
stereo.msse = mssetemp;