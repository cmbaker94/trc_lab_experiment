% Create scatter plot of mean sprd and dir spec

clear all
close all
clc

% 
% 0 deg
% angle
% ﻿{'PUV06': 0.09153443931243539,
%  'PUV11': -0.9363417420730632,
%  'wg': 0.38215473363135954,
%  'press_sz': -4.182026134360725,
%  'press_is': -0.7364310247007092,
%  'cam_sz': 2.310283594804259}
% 
% spread 
% ﻿{'PUV06': 12.483828771805785,
%  'PUV11': 12.60794971653462,
%  'wg': 15.94559308743408,
%  'press_sz': 19.32217790409021,
%  'press_is': 27.90347611894105,
%  'cam_sz': 22.90062162660438}
% 
% 40 deg
% angle
% ﻿{'PUV06': -1.0184143419250185,
%  'PUV11': -0.6216192973906992,
%  'wg': 0.1435695503004702,
%  'press_sz': 0.29081632780111705,
%  'press_is': -0.06941372628285486,
%  'cam_sz': 10.362776432496593}
% 
% spread
% ﻿{'PUV06': 16.38726105908312,
%  'PUV11': 17.72156842208246,
%  'wg': 31.73809226472663,
%  'press_sz': 36.6329927914275,
%  'press_is': 35.498922042936485,
%  'cam_sz': 34.535301712447044}

%'PUV06', 'PUV11', 'wg', 'press_sz', 'press_is', 'cam_sz'

deg0(1,:) = [0.0915, -0.9363, 0.3821, -4.182, -0.736, 2.3102];
deg0(2,:) = [12.4838, 12.6079, 15.94559, 19.322, 27.90348, 22.9006];

deg40(1,:) = [-1.01841, -0.621619, 0.143569, 0.290816, -0.06941, 10.36278];
deg40(2,:) = [16.38726, 17.7216, 31.73809, 36.633, 35.4989, 34.5353];

%%

figure('units','inches','position',[1 1 10 6],'Color','w');
% ax1 = axes('Position',[0.12 0.15 0.8 0.77]);
subplot(1,2,1)
plot([-45, 45],[-45, 45],'Color',[0.8 0.8 0.8],'LineStyle','-.','LineWidth',2)
hold on
plot([-45, 45],[0, 0],'Color',[0.8 0.8 0.8],'LineStyle','-','LineWidth',.5)
plot([0, 0],[-45, 45],'Color',[0.8 0.8 0.8],'LineStyle','-','LineWidth',.5)
scatter(0,deg0(1,1),200,'r','^','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(0,deg0(1,2),200,'r','^','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(0,deg0(1,3),200,'y','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(0,deg0(1,4),200,'r','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(0,deg0(1,5),200,'m','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(0,deg0(1,6),200,'b','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(0,deg40(1,1),200,'r','^','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(0,deg40(1,2),200,'r','^','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(0,deg40(1,3),200,'y','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(0,deg40(1,4),200,'r','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(0,deg40(1,5),200,'m','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(0,deg40(1,6),200,'b','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)

axis equal
shading flat
% grid on
ylim([-11 11])
xlim([-11 11])
%     text(27,17,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
    h1=gca;
    set(h1,'tickdir','out','xminortick','on','yminortick','on');
    set(h1,'ticklength',2*get(h1,'ticklength'));
%     set(h1,'ydir','reverse');
    set(h1,'fontsize',26);
    xlabel('Wavemaker','interpreter','latex','fontsize',20);
    ylabel('Measured','interpreter','latex','fontsize',20);
%     title('Angle of Incidence, $\theta$ ($^{\circ}$)','interpreter','latex','fontsize',20);
%     box on
% set(gca,'color','none')
    	set(h1,'ytick',[-10:10:10],'yticklabel',{'-10' '0' '10'});
        set(h1,'xtick',[-10:10:10],'yticklabel',{'-10' '0' '10'});
    
    subplot(1,2,2)
    plot([-45, 45],[-45, 45],'Color',[0.8 0.8 0.8],'LineStyle','-.','LineWidth',2)
hold on
plot([-45, 45],[0, 0],'Color',[0.8 0.8 0.8],'LineStyle','-','LineWidth',.5)
plot([0, 0],[-45, 45],'Color',[0.8 0.8 0.8],'LineStyle','-','LineWidth',.5)
plot([-45, 45],[20, 20],'Color',[0.8 0.8 0.8],'LineStyle','-','LineWidth',.5)
plot([20, 20],[-45, 45],'Color',[0.8 0.8 0.8],'LineStyle','-','LineWidth',.5)
plot([-45, 45],[40, 40],'Color',[0.8 0.8 0.8],'LineStyle','-','LineWidth',.5)
plot([40, 40],[-45, 45],'Color',[0.8 0.8 0.8],'LineStyle','-','LineWidth',.5)

    scatter(0,deg0(2,1),100,'r','^','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)

scatter(0,deg0(2,2),200,'r','^','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(0,deg0(2,3),200,'y','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(0,deg0(2,4),200,'r','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(0,deg0(2,5),200,'m','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(0,deg0(2,6),200,'b','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(40,deg40(2,1),200,'r','^','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(40,deg40(2,2),200,'r','^','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(40,deg40(2,3),200,'y','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(40,deg40(2,4),200,'r','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(40,deg40(2,5),200,'m','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
scatter(40,deg40(2,6),200,'b','fill','MarkerEdgeColor','k','MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)

axis equal
shading flat
% grid on
ylim([-5 45])
xlim([-5 45])
%     text(27,17,'$H_s~(\mathrm{m})$','interpreter','latex','fontsize',20);
    h1=gca;
    set(h1,'tickdir','out','xminortick','on','yminortick','on');
    set(h1,'ticklength',2*get(h1,'ticklength'));
%     set(h1,'ydir','reverse');
    set(h1,'fontsize',26);
    xlabel('Wavemaker','interpreter','latex','fontsize',20);
%     ylabel('Measured','interpreter','latex','fontsize',20);
%     title('Directional Spread, $\sigma_{\theta}$ ($^{\circ}$)','interpreter','latex','fontsize',20);
       	set(h1,'ytick',[0:20:40],'yticklabel',{'0' '20' '40'});
        set(h1,'xtick',[0:20:40],'yticklabel',{'0' '20' '40'});
%     box on
% set(gca,'color','none')
    
    figfolder ='/Users/cmbaker9/Documents/Research/Lab_Experiments/figures/wave_statistic/directional/'
    Sname1 = [figfolder,'scatter_dir_sprd'];
print(Sname1,'-dpng')


