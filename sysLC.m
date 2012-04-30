%This is the system to equation to simulate electrical LC system
function [xk] = sysLC(k, sysxkm1, uk)

%discretize the continuous system using Euler Method
%uk is the system noise
%fourth order Runge Kutta integration for the electronic system
%global c1 c2 L k1 k2 k3 source

c1 = 2;  %capacitor C1
c2 = 2;  %capacitor C2
L = 3; %inductor I
k1 = 3;   %for resistor R1 = 1 + k1*f2, R1(0) = 1
k2 = 3;  %for resistor R2 = 1 + k2*f6, R2(0) = 1
k3 = 3;  %for resistor R3 = 1 + k3*f11, R3(0) = 1
source = 5;  %voltage source
h = 0.15;   %step size for the euler method

r3 = 1 + k3*sysxkm1(3);
r2 = (1 + sqrt(1 + 4*k2*(sysxkm1(1) - sysxkm1(2))))/2;
r1 = (1 + sqrt(1 + 4*k1*(source - sysxkm1(1))))/2;
% if k > 200
%     r1 = 100
% else
%     r1 = (1 + sqrt(1 + 4*k1*(source - sysxkm1(1))))/2;
% end

slope(3,1) = (sysxkm1(2) - r3*sysxkm1(3))/L;
slope(2,1) = (sysxkm1(1) - sysxkm1(2))/(c1*r2) - sysxkm1(3)/c2;
slope(1,1) = -(1/(c1*r1) + 1/(c1*r2))*sysxkm1(1) + sysxkm1(2)/(c1*r2) + source/(c1*r1);

xk = sysxkm1 + h * slope + uk;   %add the noise
