function [r2,bias,rmse] = r2_bias_stats(X,Y)

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

rmse = sqrt(sum((y-x).^2)/length(x));
end
