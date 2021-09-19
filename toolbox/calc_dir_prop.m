function [dirprop] = calc_dir_prop(Sftheta,freq,dir,frange,dirrange)
% compute directional properties from 2d directional spectra (Sftheta),
% with frequencies (freq) and directions (dir). The max frequency is given
% as fmax.
rad2deg  = 180/pi;
[val, ifmin] = min(abs(freq-frange(1)));
[val, ifmax] = min(abs(freq-frange(2)));

[val, idmin] = min(abs(dir-dirrange(1)));
[val, idmax] = min(abs(dir-dirrange(2)));

% Sft = Sftheta;
% f = freq(ifmin:ifmax);

dtheta        = abs(dir(2)-dir(1));
df            = abs(freq(2)-freq(1));
Sf            = nansum(Sftheta,2)*dtheta;
Sd            = nansum(Sftheta(ifmin:ifmax,:),1)*df;

Sft           = Sftheta(:,idmin:idmax);
theta         = dir(idmin:idmax);
a1F           = (Sft*(cosd(theta).*dtheta)')./(nansum(Sft,2)*dtheta);
b1F           = (Sft*(sind(theta).*dtheta)')./(nansum(Sft,2)*dtheta);
a2F           = (Sft*(cosd(2*theta).*dtheta)')./(nansum(Sft,2)*dtheta);
b2F           = (Sft*(sind(2*theta).*dtheta)')./(nansum(Sft,2)*dtheta);

th_1          = atan2(b1F,a1F)*rad2deg; 
th_2          = 0.5*atan2(b2F,a2F)*rad2deg; 
sig_1         = (rad2deg)*sqrt(abs(2*(1-(a1F.*cosd(th_1)+b1F.*sind(th_1)))));
sig_2         = (rad2deg)*sqrt(abs(0.5*(1-(a2F.*cosd(2*th_2)+b2F.*sind(2*th_2)))));

dirprop.sig_1 = sig_1;
dirprop.th_1 = th_1;
dirprop.sig_2 = sig_2;
dirprop.th_2 = th_2;
dirprop.Sf = Sf;
dirprop.Sd = Sd;
dirprop.Sftheta = Sftheta;
dirprop.dir = dir;
dirprop.freq = freq;