%This code is a simple methodology to estimate sea-surface elevation
%spectra from pressure measurements

%% STEP 1: Clear Workspace
clear all %#ok<CLALL>

%% STEP 2: Load Pressure
fname   = 'press2.txt';
eval(['A =load(','''',fname,'''',');']);
rho     = 1000;
g       = 9.81;

%% STEP 3: Estimate eta
eta     = A/(rho*g);         % Convert pressure to meters
mp      = mean(eta);         % Get the mean depth
eta     = eta - mp;
hp      = 0.05;              % Get the height of pressure sensor above bottom
mp      = mp+hp;             % Total depth

%% STEP 4: Spectra
[Spp,F,Sppc] = pwelch(eta,4096,2048,[],100,'ConfidenceLevel',0.95);
omega        = 2*pi*F;            % Angular frequency
omega2       = omega.^2;
kh           = dispersion((omega2.*mp)/g); % Dispersion equation solution
K            = kh/mp;                    % Wave Number

%% STEP 5: Depth Correction
Kp         = (cosh(K.*hp)./cosh(kh)).^2;
Kpp        = ones(size(Kp));
maxfac     = 0.8; %3.1; %0.80;
ifreq      = find(F<maxfac);
Kpp(ifreq) = Kp(ifreq);
Sp         = Spp./Kpp;
Spc        = Sppc./Kpp;

figure
semilogy(F,Sp,'linewidth',2);
xlim([0 3]);
ylim([10^-4 10^-1.5]);
hold on
semilogy(F,Spc(:,1),'linewidth',2);
semilogy(F,Spc(:,2),'linewidth',2);
grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$S_{\eta \eta}~\mathrm{(m^{2}/Hz)}$','interpreter','latex')

%% STEP 6: Get Hs
dF          = abs(F(2)-F(1));
fmax        = maxfac;
j           = find(F<=fmax);
Hs          = 4*sqrt(nansum(Sp(j))*dF);
Hs1         = 4*sqrt(nansum(Spc(j,1))*dF);
Hs2         = 4*sqrt(nansum(Spc(j,2))*dF);

fig1 = figure('Color','w');
subplot(2,1,1)
semilogy(F,Sp,'c--','linewidth',2);
xlim([0 1.5]);
ylim([10^-4 10^-1.5]);
hold on
semilogy(F,Spp,'r-','linewidth',2);
semilogy(F,Spc(:,1),'linewidth',2);
semilogy(F,Spc(:,2),'linewidth',2);
% semilogy(freq,SSE/(100^2),'Color',[1 1 1]/2,'linewidth',2); % MM comparing...
% semilogy(freq,PP/(100^2),'k-','linewidth',2); % MM comparing...
set(gca,'FontSize',16); hold on;

grid
xlabel('$\mathrm{Freq.}~\mathrm{(Hz)}$','interpreter','latex')
ylabel('$S_{\eta \eta}~\mathrm{(m^{2}/Hz)}$','interpreter','latex')

subplot(2,1,2)
plot(F,1./Kp,'c','linewidth',2);
hold on
% plot(freq,correction.^2,'k--','linewidth',2) % MM comparing
xlim([0 1.5]);
ylim([0 25]);
ylabel('correction')
grid
set(gca,'FontSize',16); hold on;

