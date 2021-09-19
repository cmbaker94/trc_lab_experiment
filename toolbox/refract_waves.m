function sig_sz = refract_waves(Tinfo,region)

load('E:/data/processed/lidar/Riegl/TRC_bathymetry_estimate_line.mat');
g = 9.81;
h = -(h-Tinfo.tide);


theta0 = Tinfo.spread;
omega   = 2*pi/Tinfo.Tp;
cons    = omega^2*h(1)/g;
kh      = dispersi(cons);
k       = kh./h(1);
k0      = k(1);
theta   = asind(sind(theta0)*(k0./k));

for i = 2:length(h)
    k0 = k(i-1);
    theta0 = theta(i-1);
    cons   = omega^2*h(i)/g;
    kh = dispersi(cons);
    k(i) = kh/h(i);
    theta(i)   = asind(sind(theta0)*(k0./k(i)));
end

[~,id1] = min(abs(xp-region(1)));
[~,id2] = min(abs(xp-region(2)));
sig_sz = nanmean(theta(id1:id2));
end