
function [Iprop] = threshprop_combends(Icomb,Iprop,Y,xreg,iylim,resy,imageno,i)
% OUTPUT:
% yalen, yim, resx, resy, xreg, yclen, cang, yctlen, elipgeom, xycent, X,
% Y, Ithresh

[L,n] = bwlabel(Icomb);

for ni = 1:n
    obtemp = 0.*Icomb;
    obtemp(L==ni)=1;
    props = regionprops(obtemp, 'BoundingBox', 'MajorAxisLength', 'Orientation', 'Centroid', 'MinorAxisLength');
    aedge(1)            = (props.BoundingBox(2)*resy)+Y(1);
    aedge(2)            = aedge(1)+(props.BoundingBox(4)*resy);
    Iedge(ni)           = aedge;
    Iyclen(ni)          = props.MajorAxisLength*resy;
    Icang(ni)           = props.Orientation;
    Iyim(ni)            = imageno(i);
    Imajor_ax(ni)       = props.MajorAxisLength*resy;
    Iminor_ax(ni)       = props.MinorAxisLength*resy;
    Iorientation(ni)    = props.Orientation;
    Icentroid(ni)       = props.Centroid*resy+[xreg(1) iylim(1)];
    
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
    D(:,2) = C(:,2);
    D(:,1) = movmean(movmedian(C(:,1),14),14);
    Iprop.xycent{length(Iprop.yim)} = D(1:floor(length(D)/5):end,:);
    d = diff(D);
    totlen = sum(sqrt(sum(d.^2,2)));
    if size(C,1)<2
        totlen=idxdy;
    end
    Iprop.yctlen = [Iprop.yctlen; totlen];
    
    % crest ends
    Dx = D([1 end],1);
    Dy = D([1 end],2);
    idx=find(Dx<xreg(1)+0.2 | Dx>xreg(2)-0.2);
%     idy=find(Dy<-12.8 | Dy>12.8);
    crestno =  length(Dx)-length(idx);
    Ice_x(ni) = D([1 end],1);
    Ice_y(ni) = Dy = D([1 end],2);
    Icrestends(ni) = [Iprop.crestends; crestno];
    if crestno == 2
        Iexitregflag(ni) = 0;
    elseif crestno < 2
        Iexitregflag(ni) = 1;
    end
    clear C D
    
    % elipse
%     t = linspace(0,2*pi,50);
%     a = (props.MajorAxisLength*resy)/2;
%     b = (props.MinorAxisLength*resy)/2;
%     Xc = props.Centroid(1)*resy+xreg(1);
%     Yc = props.Centroid(2)*resy+min(Y,[],'all');
%     phi = deg2rad(-props.Orientation);
%     xe = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
%     ye = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
%     [~,iymin]=min(ye);
%     [~,iymax]=max(ye);
%     Xe =  [xe(iymin) xe(iymax)];
%     Ye = [ye(iymin) ye(iymax)];
%     Xe = Xc+[a*cos(phi) -a*cos(phi)];
%     Ye = Yc+[a*sin(phi) -a*sin(phi)];
%     idx=find(Xe>xreg(1) & Xe<xreg(2));
%     Iprop.crestends = [Iprop.crestends; length(idx)]; % using elliptical
%     to count  crest ends
end

for ni = 1:n
    if 
end

Iprop.yalen        = [Iprop.yalen; Iaedge];
Iprop.yclen       = [Iprop.yclen; Iyclen];
Iprop.cang        = [Iprop.cang; Icang];
Iprop.yim         = [Iprop.yim Iyim];
Iprop.elipgeom.major_ax    = [Iprop.elipgeom.major_ax; Imajor_ax];
Iprop.elipgeom.minor_ax	= [Iprop.elipgeom.minor_ax; Iminor_ax];
Iprop.elipgeom.orientation = [Iprop.elipgeom.orientation; Iorientation];
Iprop.elipgeom.centroid    = [Iprop.elipgeom.centroid; Icentroid];

Iprop.exitregflag = [Iprop.exitregflag; Iexitflag];
Iprop.crestends = [Iprop.crestends; Icrestends]; % using elliptical

% Iprop.crestends = crestends;
end