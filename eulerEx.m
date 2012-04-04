%Euler method to do numerical integration
%dy/dt = 3*exp(-4*t)-2*y(t), y(0) = 1
% function [xx, yy] = eulerEx
% 
% clc
% close all
% clear all
% 
% 
% h = 0.1;  %h is the time step
% t = 0:h:10; %initialize time variable
% 
% y = 1.0
% 
% for i = 1:(length(t)-1)
%     k1 = EulerM(y(i), t(i));
%     %k1 = 3*exp(-4*t(i)) - 2*y(i);
%     y(i+1) = y(i) + h*k1;
% end
% 
% yexact = 2.5*exp(-2*t)-1.5*exp(-4*t);
% 
% plot(t, y, 'b--', t,yexact,'r-');
% legend('Approximate', 'exact');
% title('Euler Approxiamtion, h = 0.1');
% xlabel('Time');
% ylabel('y*(t), y(t)');
% 
% 
% function [slope] = EulerM(x, t)
% slope = 3*exp(-4*t) - 2*x

clc 
clear all
close all

h = 0.1;
t = 0:h:5;

c1 = 2;  %capacitor C1
c2 = 2;  %capacitor C2
L = 1; %inductor I
k1 = 3;   %for resistor R1 = 1 + k1*f2, R1(0) = 1
k2 = 0.5;  %for resistor R2 = 1 + k2*f6, R2(0) = 1
k3 = 3;  %for resistor R3 = 1 + k3*f11, R3(0) = 1
source = 5;  %voltage source

x=[0;0;0];
xArray = x;

for i = 1:(length(t)-1)
    r3 = 1 + k3*x(3)
    r2 = (1 + sqrt(1 + 4*k2*(x(1) - x(2))))/2
    r1 = (1 + sqrt(1 + 4*k1*(source - x(1))))/2

    slope(3,1) = (x(2) - r3*x(3))/L;
    slope(2,1) = (x(1) - x(2))/(c1*r2) - x(3)/c2;
    slope(1,1) = -(1/(c1*r1) + 1/(c1*r2))*x(1) + x(2)/(c1*r2) + source/(c1*r1)

    
    x = x + h*slope
    xArray = [xArray x];
end


    







