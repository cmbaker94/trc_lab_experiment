
function [Iprop] = threshprop(Icomb,Iprop,Y,xreg,iylim,resy,imageno,i,distconnect)
% OUTPUT:
% yalen, yim, resx, resy, xreg, yclen, cang, yctlen, elipgeom, xycent, X,
% Y, Ithresh

Icombconnect  = bwdist(Icomb)<=distconnect;
[Lconnect,n] = bwlabel(Icombconnect);
% [L,n] = bwlabel(Icomb);
L = Lconnect.*Icomb; 

for ni = 1:n
    obtemp = 0.*Icomb;
    obtemp(L==ni)=1;
%     propstemp = regionprops(obtemp, 'BoundingBox', 'MajorAxisLength', 'Orientation', 'Centroid', 'MinorAxisLength');
%     edgetemp = round(propstemp.BoundingBox([2 4]));
%     obtemp(1:edgetemp(1)+distconnect,:) = 0;
%     obtemp(edgetemp(2)-distconnect:end,:) = 0;
    
    props = regionprops(obtemp, 'BoundingBox', 'MajorAxisLength', 'Orientation', 'Centroid', 'MinorAxisLength');
    aedge(1)    = (props.BoundingBox(2)*resy)+Y(1);
    aedge(2)    = aedge(1)+(props.BoundingBox(4)*resy);
    Iprop.yalen        = [Iprop.yalen; aedge];
    Iprop.yclen       = [Iprop.yclen; props.MajorAxisLength*resy];
    Iprop.cang        = [Iprop.cang; props.Orientation];
    Iprop.yim         = [Iprop.yim imageno(i)];
    Iprop.elipgeom.major_ax    = [Iprop.elipgeom.major_ax props.MajorAxisLength*resy];
    Iprop.elipgeom.minor_ax	= [Iprop.elipgeom.minor_ax props.MinorAxisLength*resy];
    Iprop.elipgeom.orientation = [Iprop.elipgeom.orientation props.Orientation];
    ctemp       = props.Centroid*resy+[xreg(1) iylim(1)];
    Iprop.elipgeom.centroid    = [Iprop.elipgeom.centroid; ctemp];
    
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
    Iprop.crestends = [Iprop.crestends; crestno]; % using elliptical
    if crestno == 2
        Iprop.exitregflag = [Iprop.exitregflag; 0];
    elseif crestno < 2
        Iprop.exitregflag = [Iprop.exitregflag; 1];
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
% Iprop.crestends = crestends;
end