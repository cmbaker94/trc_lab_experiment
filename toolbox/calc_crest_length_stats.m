function [output] = calc_crest_length_stats(input,Wtank,Cnorm,dx,xrange)
% Fuction to compute the wave crest statistics 
% INPUT:
% input structure read in from either an image or stereo ylen file
% Wtank: Width of tank
% Cnorm: number of frames to normalize by
% dx: number of bins for histogram
% xrange: cross-shore extent
% OUPUT:
% structure of crest length statistics

nearedge = 26;
maxclen = sqrt(5^2+Wtank^2);

% format info
yalen = abs(input.yalen(:,1)-input.yalen(:,2));
yclen = input.yclen;
ytlen = input.yctlen;
cang = input.cang;

% remove crests outside of region
yalen(input.crestends<2)=NaN;
yclen(input.crestends<2)=NaN;
ytlen(input.crestends<2)=NaN;
cang(input.crestends==0)=NaN;
% assess lengths > Wtank
il = input.yalen(:,1)<-nearedge/2;
ir = input.yalen(:,2)>nearedge/2;
% cw = cang(yalen>Wtank);
cw = cang(yclen>maxclen);
yclen(yclen>maxclen)= Wtank*cosd(abs(cw)-90);
% yclen(yalen>Wtank)= Wtank*cosd(abs(cw)-90); % if yalen>Wtank, then take then overwrite with Wtank *(cos(abs(cang))-90)
% yalen(yalen>Wtank)=Wtank;

% store angle output
cang(cang<0) = cang(cang<0)+180;
output.cang     = cang;

% rename and prep
% yalen(cang>150) = NaN;
% yalen(cang<30) = NaN;
ya = yalen;
% yclen(cang>150) = NaN;
% yclen(cang<30) = NaN;
yc = yclen;
%     ya = yalen(max([il'; ir'])==0);
%     ya(ya>Wtank)=Wtank;
yt = ytlen;
% cang(cang>150) = NaN;
% cang(cang<30) = NaN;
output.cangstd = nanstd(cang);
display(output.cangstd)

% count crest numbers
Cno = sum(input.crestends)-sum(il)-sum(ir);
% Cno = 2*length(yalen)-sum(il)-sum(ir);
output.ends   = Cno/Cnorm;

% compute Nmin/Nmax
Nmin = Cno/26.5;
Nmax = (2*length(input.yim)-sum(input.crestends>1))/xrange;
Nratio = Nmin/Nmax;
% Nmin = Cno;
% Nmax = (2*length(input.yim)-sum(input.crestends>1));
% Nratio = Nmin/Nmax;
% count  = 1;
% for iim =  min(input.yim):max(input.yim)
%     iy = find(input.yim==iim);
%     ilim = input.yalen(iy,1)<-nearedge/2;
%     irim = input.yalen(iy,2)>nearedge/2;
%     Cnoim = sum(input.crestends(iy))-sum(ilim)-sum(irim);
%     Nmin = Cnoim/26.5;
%     Nmax = (2*length(input.yim(iy))-sum(input.crestends(iy)>1))/xrange;
% %     Nratioim(count) = Nmin/Nmax;
%     if Nmax > 0
%         Nratioim(count) = Nmin/Nmax;
%         count = count+1;
%     end
% end
% [xbins,nums]     = bin_data(Nratioim,0.01);
% Nratio = mode(xbins);

display(Nratio)
output.Nratio = Nratio;

% alongshore component
output.ya        = yalen;
output.meana     = nanmean(ya);
output.meda      = nanmedian(ya);
output.modea     = mode(ya);
[xbins,nums]     = bin_data(yalen,dx);
output.numbina   = nums;
output.bina      = xbins;

% along-crest eliptical component
output.yc        = yclen;
output.meanc     = nanmean(yc);
output.medc      = nanmedian(yc);
output.modec     = mode(yc);
[xbins,nums]     = bin_data(yclen,dx);
output.numbinc   = nums;
output.binc      = xbins;

% along-crest following contour component
output.yt        = ytlen;
output.meant     = nanmean(yt);
output.medt      = nanmedian(yt);
output.modet     = mode(yt);
ytlen(ytlen>Wtank) =  Wtank;
display('Changing ytlen to equal to tank width')
[xbins,nums]     = bin_data(ytlen,dx);
output.numbint   = nums;
output.bint      = xbins;
dang = 3;
[xbins,nums]     = bin_data(cang,dang);
output.numbinang   = nums;
output.binang      = xbins;

output.yim = input.yim;

% 2dpdf
% output.pdf = mvnrnd(yclen,cang,dang);

