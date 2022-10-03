% total variance vs crest length
clear all
clc
close all
datafolder = 'E:\data\processed\trial_comparison\';
load([datafolder,'wave_crest_Hs25.mat'])
load([datafolder,'total_variance.mat'])
%%
sprd = [0, 10, 20, 30, 40];
szavg = nanmean([szout szin],2);
cm = cmocean('thermal');
cmo = cm(1:floor(255/4):256,:);

figure('units','inches','position',[1 1 8 5],'color','w')
ax1 = axes('Position',[0.18 0.15 0.65 0.75]);
    plot(Hs25.stend,szavg,'k','LineWidth',2)
%     plot(Hs25.stend,nanmean(szout,2),'k','LineWidth',2)
    hold on
%     plot(Hs25.stend,szavg,'r','LineWidth',2)
%     errorbar(sprd,nanmean(isout,2),nanstd(isout,[],2),'Color','k','LineWidth',1.5,'LineStyle','None')
    scatter(Hs25.stend,szavg,200,sprd,'fill','MarkerEdgeColor','k','LineWidth',1)
%     scatter(Hs25.stend,nanmean(szout,2),200,sprd,'fill','MarkerEdgeColor','k','LineWidth',1)
    
%     scatter(Hs25.stend,szavg,70,'r','fill','MarkerEdgeColor','k','LineWidth',1)
%     grid on
    h1=gca;
    % set(h1, 'XScale', 'log')
    set(h1,'tickdir','in','xminortick','on','yminortick','on');
    set(h1,'ticklength',1.5*get(h1,'ticklength'));
    set(h1,'fontsize',20);
    box on
    colormap(cmo)
    caxis([-4 44])
    hc = colorbar('Position', [0.86 0.15 0.05 0.48]);
    
%         colormap(ax3,flipud(cmocean('deep')));
%     hc = colorbar(ax3,'Position', [0.76 0.095 0.05 0.74]);
%     caxis(ax3,[-0.15001 0.15001])
%     xlim([-2 42])
%     ylim([0 4.2]*10^-3)
    xlabel('$\#$ Crest Ends','interpreter','latex','fontsize',22)
%     ylabel([{'Surfzone Avg.,'},{'$\sigma_{vel}$ (m$^2$/s$^2$)'}],'interpreter','latex','fontsize',22)
    ylabel('$\langle \sigma_{vel} \rangle _{\mathrm{sz}}$ (m$^2$/s$^2$)','interpreter','latex','fontsize',22)
    text(15.5, 3.15*10^-3,'$\sigma_{\theta}~(^{\circ})$','interpreter','latex','fontsize',24)
    % set(h1,'yticklabel',[],'xticklabel',[]);
    % text(ax3,0.05, 10^(-1.3),'(c)','interpreter','latex','fontsize',24);
    % title(ax3,'$\sigma_{\theta}=40^{\circ}$','interpreter','latex','fontsize',28);
    % errorbar(T40.specon.wg.f(67),10^(-1.3),10^-1.85,'Color','k','LineWidth',2.5)%wg.See(3)-wg.Seec(3,1))
%     h2 = legend('outer surf zone','surfzone avg','interpreter','latex','fontsize',22)
%     set(h2,'orientation','vertical','Location','northeast')
   sname = ['E:\figures\meas_comp\velocity\vel_power_vs_numcrest'];
    print([sname],'-dpng')
