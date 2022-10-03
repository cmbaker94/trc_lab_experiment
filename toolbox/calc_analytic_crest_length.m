function [yc,Cno,Lsz] = calc_analytic_crest_length(Tinfo,spread,xreg,sprdcon,T,calcmethod)
% This code will calculate the crest lenths and number of crest ends.

% CREST LENGTH:
% Assuming 2 wave trains with wave angles of spread/2 and -spread/2. This 
% is assumed to be constant in the cross-shore. Therefore, the offshore properties 
% are used. 

% NUMBER OF CREST ENDS:
% The number of crests is computed as 2*(width of tank/crest length) with 
% one crest end added with (width of tank/crest length) is a fraction
% (Note, this accounts for waves that have one edge of the crest touching
% the wall). This code assumes ~2 waves can fit in the region at one time. 
% This is based on the region it is computed over (for example, x=25-31.5 m = 6.5 m) and the
% average wave length within this range (for example, Lavg = 3.27 m).

% INPUT:
% Tinfo = trial information
% xreg  = cross-shore distance computed over
% sprdcon = spreading constant to multiply spread by to get theta angle
% OUPUT:
% yc = crest length
% Cno = number of crest ends

% defaults
Wtank   = 26.5; % width of tank
g       = 9.81; % gravity
load('E:/data/processed/lidar/Riegl/TRC_bathymetry_estimate_line.mat');
h = -(h-Tinfo.tide);
% T = 1.5434;
% Tm = Tp/1.2958;

% Calc yc
theta = spread*sprdcon;
omega   = 2*pi/T;
cons    = omega^2*h(1)/g;
kh      = dispersi(cons);
k       = kh/h(1);
L      = 2*pi/k;
if calcmethod == 0
    yp      = L/(sind(theta)-sind(-theta)); % Dalrymple
elseif calcmethod == 1
    yp = L/deg2rad(theta); % L-H
end
yc      = 0.5*yp;

% number of crest ends per cross-shore location
% Cno = floor(2*Wtank/yp);
Cno = 2*Wtank/yp;
% if Wtank/yp > floor(Wtank/yp)
%     Cno = Cno+1;
% end

% Calc Lavg
for i = 2:length(h)
    cons   = omega^2*h(i)/g;
    kh = dispersi(cons);
    k(i) = kh/h(i);
end
Len      = 2*pi./k; % crest length as a function of cross-shore location
% Compute the average wave length
[~,ioff] = nanmin(abs(xp-xreg(1)));
[~,ion] = nanmin(abs(xp-xreg(2)));
Lavg = nanmean(Len(ioff:ion)); % average crest length in region
% number of crest lengths per region
nosz = (xreg(2)-xreg(1))/Lavg; % number of wave in surf zone

% Total crest ends per region
Cno = (Cno*nosz);

% Remove values for 0 deg case
if isinf(yc)
    yc = NaN;
    Cno = NaN;
end

Lsz = Lavg;




