function [xbins,nums] = bin_data(x,dx)

xbins = [0:dx:max(x)];
xin = round(x/dx)*dx;
% if upvdown == 0 %down
%     xin = floor(x) + floor( (x-floor(x))/dx) * dx;
% elseif upvdown == 1 %up
%     xin = floor(x) + ceil( (x-floor(x))/dx) * dx
% end

nums = zeros(size(xbins));
for i = 1:length(xbins)
    nums(i) = sum(xbins(i)==xin);
end