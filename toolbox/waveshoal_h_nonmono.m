function wave = waveshoal_h_nonmono(T, H0, theta0, gamma, h)
% shoal a wave plus depth-limited breaking, Melissa Moulton

% Constants
g=9.81; % m/s^2

% Calculate wavelengths in deep water at depths h

% Deep water wavelength:
Ldeep = g*T.^2/(2*pi); 

% Wavelength, Ole Madsen approx:
L = Ldeep.*(1-exp(-(2.*pi.*h./Ldeep).^1.25)).^0.4;

% Calculate group and phase speeds at depths h
c = L./T; % Phase speed
k = 2.*pi./L; % Wavenumber
cg = (L./(2.*T)).*(1+2.*(k).*h./sinh(2.*(k).*h)); % Group velocity

% Calculate group and phase speeds at depth h0
c0 = c(1); % Phase speed at depth h0
cg0 = cg(1); % Phase speed at depth h0
   
% Compute wave height and angle at depths h
theta = asind(c./c0.*sind(theta0));
H = H0*sqrt(cg0./cg).*sqrt(cosd(abs(theta0))./cosd(abs(theta)));

% Calculate breaking variables
breaking_index = find(H./h>gamma);
breaking_index = breaking_index(1);
breaking_depth = h(breaking_index);
breaking_height = H(breaking_index);
breaking_angle = theta(breaking_index);

% NaN values above shoreline
theta(h<0) = NaN;
H(h<0) = NaN;
L(h<0) = NaN;
c(h<0) = NaN;
cg(h<0) = NaN;

% Store variables
wave.h = h;
wave.L = L;
wave.Ldeep = Ldeep;
wave.H = H;
wave.c = c;
wave.cg = cg;
wave.theta = theta;
wave.breaking_depth = breaking_depth;
wave.breaking_height = breaking_height;
wave.breaking_angle = breaking_angle;
wave.breaking_index = breaking_index;

% Compute profile onshore of breaking

H(breaking_index) = h(breaking_index)*gamma;
binds = breaking_index:(length(h)-1);

for ii=binds
    if h(ii+1)<h(ii)
        % Depth-limited breaking
        H(ii+1) = h(ii+1)*gamma;
    else
        % Re-shoaling)
        h_2 = h(ii:end);
        H0_2 = H(ii);
        theta0_2 = theta(ii);
        
        wave_2 = waveshoal_subf(T, H0_2, theta0_2, gamma, h_2);
        
        H_2 = wave_2.H;
        theta_2 = wave_2.theta;
        breaking_index2 = wave_2.breaking_index;
        
        % Fill in values
        H(ii:end) = H_2;
        theta(ii:end) = theta_2;
        
        % Second breaking region
        H(ii+breaking_index2:end) = h(ii+breaking_index2:end)*gamma;        
        
        break
    end
end

% NaN values above shoreline
theta(h<0) = NaN;
H(h<0) = NaN;

wave.H = H;
wave.theta = theta;

end

function wave_sub = waveshoal_subf(T, H0, theta0, gamma, h)
% subfunction for re-shoaling

% Constants
g=9.81; % m/s^2

% Calculate wavelengths in deep water at depths h

% Deep water wavelength:
Ldeep = g*T.^2/(2*pi); 

% Wavelength, Ole Madsen approx:
L = Ldeep.*(1-exp(-(2.*pi.*h./Ldeep).^1.25)).^0.4;

% Calculate group and phase speeds at depths h
c = L./T; % Phase speed
k = 2.*pi./L; % Wavenumber
cg = (L./(2.*T)).*(1+2.*(k).*h./sinh(2.*(k).*h)); % Group velocity

% Calculate group and phase speeds at depth h0
c0 = c(1); % Phase speed at depth h0
cg0 = cg(1); % Phase speed at depth h0

% Compute wave height and angle at depths h
theta = asind(c./c0.*sind(theta0));
H = H0*sqrt(cg0./cg).*sqrt(cosd(abs(theta0))./cosd(abs(theta)));

% NaN values above shoreline
theta(h<0) = NaN;
H(h<0) = NaN;
L(h<0) = NaN;
c(h<0) = NaN;
cg(h<0) = NaN;

% Calculate breaking variables
breaking_index = find(H./h>gamma);
breaking_index = breaking_index(1);
breaking_depth = h(breaking_index);
breaking_height = H(breaking_index);
breaking_angle = theta(breaking_index);

% Store variables
wave_sub.h = h;
wave_sub.L = L;
wave_sub.Ldeep = Ldeep;
wave_sub.H = H;
wave_sub.c = c;
wave_sub.cg = cg;
wave_sub.theta = theta;
wave_sub.breaking_depth = breaking_depth;
wave_sub.breaking_height = breaking_height;
wave_sub.breaking_angle = breaking_angle;
wave_sub.breaking_index = breaking_index;

end