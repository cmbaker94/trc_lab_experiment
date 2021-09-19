function [registrationError] =  fit1dregerror(x1,y1,x2,y2,pltflg)
% function to find cloest line
% from here: https://www.mathworks.com/matlabcentral/answers/484637-from-a-set-of-curves-how-to-find-the-one-that-the-most-similar-to-the-reference-curve

addpath('C:\Users\cmbaker9\Documents\MATLAB\MTOOLS')

% close all
% x1=linspace(0,11,100);
% y1=sin(x1);
%     A=1.5;  a=1.1; t=0.5;  %stretch parameters
% x2=a*(linspace(3,12,150)+t);  y2=A*sin(x2);
% x3=x1; y3=x3/20;
% reg1d(x1,y1,x2,y2)
% reg1d(x1,y1,x3,y3)


F=griddedInterpolant(x1,y1,'cubic');
fun=@(p,xdata) F(xdata/p(1)-p(2));
flist={fun};
a0=(max(x2)-min(x2))/(max(x1)-min(x1));
t0=mean(x2/a0)-mean(x1);
[pr,Ar]=fminspleas(flist,[a0,t0],x2,y2);
efun=@(x) Ar*fun(pr,x);

registrationError=norm(efun(x2)-y2)

if pltflg == 1
    figure;
    subplot(2,1,1)
    plot(x1,y1,'x',x2,y2);
    title 'INITIAL';legend('reference','test');
    set(gca,'FontSize',15);
    subplot(2,1,2)
    plot(x2,efun(x2),'x', x2,y2);
    title("ALIGNED: Error ="+registrationError);
    legend('deformed reference','test');
    set(gca,'FontSize',15);
end

end