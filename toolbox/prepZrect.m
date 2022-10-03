function  [Z]  = prepZrect(Tinfo,stereox,stereoy,stereoz,X,Y)
% Code to iterpolate and fill missing data in stereo to use for rectifying
% images
% INPUT:
% Tinfo: trial information
% stereox,y,z: x,y,z information at the same instance in time as the image
% X,Y: XY position to interpolate to
% OUTPUT:
% Z elevation gridded to X,Y with no nans
stereoz(stereoz<0.75) = NaN;
Z = interp2(stereox,stereoy,stereoz,X,Y,'linear');

% create 2-d bathymetry
[x,y,~,z]=domain2Dlab(Tinfo);
z = interp2(x,y,z,X,Y);
z = smoothdata(z,1,'movmean',[6 6],'omitnan');%[5 5]
z = smoothdata(z,2,'movmean',[5 5],'omitnan');

xtemp = X(1,:);
ytemp = Y(:,1);
if max(xtemp) > 34.5
    xr = [27 30.4 34.5 34.5 30 27 27];
    yr = [-13.2 -13.2 -11 10.8 13.2  13.2 -13.2];
elseif max(xtemp) == 33.5
    xr = [27 30.4 33.5 33.5 30 27 27];
    yr = [-13.2 -13.2 -11.8 11.6 13.2  13.2 -13.2];
end
for i = 1:length(xr)
    iX(i)  = find(round(xtemp,2)==xr(i));
    iY(i) = find(round(ytemp,2)==yr(i));
end
stereoregion = poly2mask(iX,iY,length(ytemp),length(xtemp));

Z(~stereoregion) = z(~stereoregion);
% znan = ~isnan(Z);
% Z1  = griddata(X(znan),Y(znan),Z(znan),X,Y);
Z = fillmissing(Z,'linear');

Z = smoothdata(Z,1,'movmean',[6 6],'omitnan');%[5 5]
Z = smoothdata(Z,2,'movmean',[5 5],'omitnan');% [3 3]