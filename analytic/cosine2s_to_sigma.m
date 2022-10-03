% calculating std of cosine2s distrubtion
% Code created by C.M. Baker, Sept 2021

rad = deg2rad([-90:1:90]);
S = 1:0.05:171; % full range of s distributions in increments of 1
TH = 0;
pltflg = 0;
for i = 1:length(S)
    sdist = gamma(S(i)+1)/2/sqrt(pi)./gamma(S(i)+1/2).*cos((rad-TH)/2).^(2*S(i));
    left = trapz(rad,sdist.*cos(rad))^2;
    right = trapz(rad,sdist.*sin(rad))^2;
    sprd(i) = rad2deg(sqrt(2*(1-sqrt(left+right))));
    sigma(i) = rad2deg(std(rad,sdist));
end

figure('Color','w')
plot(sprd,sigma,'LineWidth',2,'Color','r')
hold on
plot([0:60],[0:60],'LineStyle','-.','LineWidth',2,'Color',[0.5 0.5 0.5])

xlabel('$\sigma_{\theta}~(^{\circ})$','interpreter','latex')
ylabel('std($\cos-2s)$','interpreter','latex')
% xlim([-180 180]);
% ylim([10^(-6) 10^(-2)])
h1=gca;
set(h1,'tickdir','in','xminortick','on','yminortick','on');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
% set(h1,'YScale','log')
set(h1,'fontsize',15);