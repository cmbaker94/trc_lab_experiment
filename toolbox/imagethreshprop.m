function [IMc] = imagethreshprop(IM,IMc,Ithresh,X,Y,thresh,idxmin,idxmax,minarea,xreg,ixlim,iylim,resy,pltflag,imageno,i,Tinfo)
% OUTPUT:
% yalen, yim, resx, resy, xreg, yclen, cang, yctlen, elipgeom, xycent, X,
% Y, Ithresh

Irbw = smoothdata(IM,1,'movmedian',[3 3],'omitnan');
Irbw = smoothdata(Irbw,2,'movmedian',[1 1],'omitnan');
Irbw = Irbw(:,idxmin:idxmax);
% apply the threshold to
Ib   = imbinarize(Irbw,thresh);
Iob = bwareaopen(Ib,minarea,8);
Ithresh(:,:,i) = Iob;
%     Iob = bwareafilt(Iob,5);
[L,n] = bwlabel(Iob);

shapeheight = [];
for ni = 1:n
    obtemp = 0.*Iob;
    obtemp(L==ni)=1;
    props = regionprops(obtemp, 'BoundingBox', 'MajorAxisLength', 'Orientation', 'Centroid', 'MinorAxisLength');
    aedge(1)    = (props.BoundingBox(2)*resy)+Y(1);
    aedge(2)    = aedge(1)+(props.BoundingBox(4)*resy);
    IMc.yalen        = [IMc.yalen; aedge];
    IMc.yclen       = [IMc.yclen; props.MajorAxisLength*resy];
    IMc.cang        = [IMc.cang; props.Orientation];
    IMc.yim         = [IMc.yim imageno(i)];
    IMc.elipgeom.major_ax    = [IMc.elipgeom.major_ax props.MajorAxisLength*resy];
    IMc.elipgeom.minor_ax	= [IMc.elipgeom.minor_ax props.MinorAxisLength*resy];
    IMc.elipgeom.orientation = [IMc.elipgeom.orientation props.Orientation];
    ctemp       = props.Centroid*resy+[xreg(1) iylim(1)];
    IMc.elipgeom.centroid    = [IMc.elipgeom.centroid; ctemp];
    
    count = 0;
    for ro = 1 : size(obtemp, 1)
        if sum(obtemp(ro,:))>0
            count = count+1;
            B(1) = (find(obtemp(ro, :), 1, 'first')*resy)+xreg(1);
            B(2) = (find(obtemp(ro, :), 1, 'last')*resy)+xreg(1);
            C(count, 1) = mean(B);
            C(count, 2) = (ro*resy)+iylim(1);
        end
        clear B
    end
    IMc.xycent{length(IMc.yim)} = C(1:floor(length(C)/5):end,:);
    d = diff(C);
    totlen = sum(sqrt(sum(d.^2,2)));
    IMc.yctlen = [IMc.yctlen; totlen];
    clear C
end
if pltflag == 1
    %         figure('units','inches','position',[1 1 5 8],'color','w')
    % plot
    ax2 = axes('Position',[0.15 0.1 0.4 0.8]);
    pcolor(X,Y,Irbw)%(1:end-1,:))
    %     imagesc(nu,nv,Irbw)
    hold on
    props = regionprops(Iob, 'BoundingBox', 'MajorAxisLength', 'Orientation', 'Centroid', 'MinorAxisLength');
    
    t = linspace(0,2*pi,50);
    for ni = 1:length(props)
        a = (props(ni).MajorAxisLength*resy)/2;
        b = (props(ni).MinorAxisLength*resy)/2;
        Xc = props(ni).Centroid(1)*resy+xreg(1);
        Yc = props(ni).Centroid(2)*resy+iylim(1);
        phi = deg2rad(-props(ni).Orientation);
        x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
        y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
        plot(x,y,'r','Linewidth',5)
    end
    shading interp
    colormap('gray');
    %     hc = colorbar('Position', [0.88 0.1 0.03 0.73]);
    caxis([0 1])
    axis equal
    box on
    xlim(ixlim)
    ylim(iylim)
    h1=gca;
    %         set(h1,'ydir','reverse')
    set(h1,'tickdir','out','xminortick','on','yminortick','on');
    set(h1,'ticklength',1*get(h1,'ticklength'));
    set(h1,'fontsize',16);
    set(h1,'ytick',[-10:5:10],'yticklabel',{'-10' '-5' '0' '5' '10'});
    ylabel('$y$ (m)','interpreter','latex','fontsize',20);
    xlabel('$x$ (m)','interpreter','latex','fontsize',20);
    
    ax3 = axes('Position',[0.46 0.1 0.45 0.8]);
    pcolor(X,Y,double(Iob))%(1:end-1,:)))
    hold on
    shading faceted
    colormap('gray');
    % hc = colorbar('Position', [0.575 0.1 0.02 0.3]);
    caxis([0 1])
    axis equal
    %     grid on
    box on
    shading interp
    xlim(xreg)
    ylim(iylim)
    h1=gca;
    %         set(h1,'ydir','reverse')
    set(h1,'tickdir','out','xminortick','on','yminortick','on');
    set(h1,'ticklength',1*get(h1,'ticklength'));
    set(h1,'fontsize',16);
    set(h1,'ytick',[-10:5:10],'yticklabel',{'' '' '' '' ''});
    text(ax2,13.95,26.25,'$\eta$ (m)','interpreter','latex','fontsize',20);
    xlabel('$x$ (m)','interpreter','latex','fontsize',20);
    text(ax2,25,17.5,['Threshold~ = ',num2str(thresh)],'interpreter','latex','fontsize',15);
    text(ax2,25,16.5,['Min Area~~ = ',num2str(minarea),' pix'],'interpreter','latex','fontsize',15);
    text(ax2,25,15.5,['Frame ~~~~~ = ',num2str(imageno(i))],'interpreter','latex','fontsize',15);
    
    sname = ['thresh_crest_',num2str(i,'%04.f')];
    print([Tinfo.figfolder,'images/',sname],'-dpng')
    pause(0.1)
    clf
end