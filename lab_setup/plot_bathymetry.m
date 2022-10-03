%% Load bathymetry and analytic solution

for bathyanalysis = 1
    
    analytic.h = sz.Tinfo.stillwat-bathymetry.h;
    analytic.x = bathymetry.xp;
    
    % wave conditions, gamma
    analytic.T = 2; % period, s
    analytic.H0 = 0.3; % offshore wave height, m
    analytic.theta0 = 0; % offshore wave direction
    analytic.gamma = 0.8; % breaking wave gamma
    
    % run the wave shoaling/breaking code
    wave = waveshoal_h_nonmono(analytic.T, analytic.H0, analytic.theta0, analytic.gamma, analytic.h);
    
    % plot bathymetry
    figure('units','inches','position',[1 1 10 3],'color','w');
    plot(analytic.x,-analytic.h,'k','LineWidth',2)
    hold on
    plot(analytic.x,analytic.h*0,'k--','LineWidth',2);
    fill([analytic.x; flipud(analytic.x)],[-analytic.h; -2*ones(size(analytic.h))],[0.7 0.7 0.7])
    xlim([20 34])
    ylim([-1.1 0.1])
    xlabel('$x$ (m)','interpreter','latex','fontsize',18);
    ylabel('$h$ (m)','interpreter','latex','fontsize',18);
    grid on
    h1=gca;
    set(h1,'tickdir','in','xminortick','on','yminortick','on');
    set(h1,'ticklength',1*get(h1,'ticklength'));
    set(h1,'fontsize',20);
    Sname2 = ['/Users/cmbaker9/Documents/Research/Lab_Experiments/figures/LIDAR/Riegl/bathymetry/TRC_bathymetry'];
    print(Sname2,'-dpng')
    Sname2 = [figfolder,'TRC_bathymetry'];
    print(Sname2,'-dpng')
    
end
