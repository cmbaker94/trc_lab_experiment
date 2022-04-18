function [r2,bias,rmse,pfit,slope,r] = r2_bias_stats(X,Y)

notnan      = find(~isnan(Y));
x           = X(notnan);
y           = Y(notnan);
p = polyfit(x,y,1);
yfit        = polyval(p,x);
yresid      = y - yfit;
SSresid     = sum(yresid.^2);
SStotal     = (length(y)-1) * var(y);
r2         = 1 - SSresid/SStotal;
bias  = sum(y-x)/length(x);

cc = corrcoef(x,y);
r = cc(1,2);

rmse = sqrt(sum((y-x).^2)/length(x));
pfit.x = x;
pfit.y=yfit;
if ~isnan(rmse)
    slope = (yfit(end)-yfit(1))/(x(end)-x(1));
else
    slope = NaN;
end

end
