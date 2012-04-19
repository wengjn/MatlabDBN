%Euler method to do numerical integration
%x(1)---e4; x(2)---e8; x(3)---f10
%
function [slope] = EulerMethod(x)
%fourth order Runge Kutta integration for the electronic system
%global c1 c2 L k1 k2 k3 source
c1 = 2;  %capacitor C1
c2 = 2;  %capacitor C2
L = 3; %inductor I
k1 = 3;   %for resistor R1 = 1 + k1*f2, R1(0) = 1
k2 = 3;  %for resistor R2 = 1 + k2*f6, R2(0) = 1
k3 = 3;  %for resistor R3 = 1 + k3*f11, R3(0) = 1
source = 5;  %voltage source

r3 = 1 + k3*x(3);
r2 = (1 + sqrt(1 + 4*k2*(x(1) - x(2))))/2;
r1 = (1 + sqrt(1 + 4*k1*(source - x(1))))/2;

slope(3,1) = (x(2) - r3*x(3))/L;
slope(2,1) = (x(1) - x(2))/(c1*r2) - x(3)/c2;
slope(1,1) = -(1/(c1*r1) + 1/(c1*r2))*x(1) + x(2)/(c1*r2) + source/(c1*r1);
