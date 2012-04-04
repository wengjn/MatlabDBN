function [dx1, dx2, dx3, dx4, dx5] = RungeKutta(x)
%fourth order Runge Kutta integration for the electronic system
%global c1 c2 L k1 k2 k3 source
c1 = 2e-6;  %capacitor C1
c2 = 10e-6;  %capacitor C2
L = 1e-3; %inductor I
k1 = 3;   %for resistor R1 = 1 + k1*f2, R1(0) = 1
k2 = 0.5;  %for resistor R2 = 1 + k2*f6, R2(0) = 1
k3 = 3;  %for resistor R3 = 1 + k3*f11, R3(0) = 1
source = 5;  %voltage source
dt = 0.5;

dx1(3,1) = x(2)/L - (1+k3*x(3))*x(3)/L;
r2 = (1+sqrt(1+4*k2*(x(1)-x(2))))/2;
r1 = (1+sqrt(1+4*k1*(source-x(1))))/2;
dx1(2,1) = x(1)/(c2*r2) - x(2)/(c2*r2) - x(3)/c2;
dx1(1,1) = -x(1)*(1/(c1*r1)+1/(c1*r2)) + x(2)/(c1*r2) + source/(c1*r1);

dx1 = dx1*dt;
xtemp = x + dx1 / 2;
dx2(3,1)= xtemp(2)/L - (1+k3*xtemp(3))*xtemp(3)/L;
%r2 = (1+sqrt(1+4*k2*(xtemp(1)-xtemp(2))))/2;
%r1 = (1+sqrt(1+4*k1*(source-xtemp(1))))/2;
dx2(2,1) = xtemp(1)/(c2*r2) - xtemp(2)/(c2*r2) - xtemp(3)/c2;
dx2(1,1) = -xtemp(1)*(1/(c1*r1)+1/(c1*r2)) + xtemp(2)/(c1*r2) + source/(c1*r1);

dx2 = dx2 * dt;
xtemp = x + dx2 / 2;
dx3(3,1)= xtemp(2)/L - (1+k3*xtemp(3))*xtemp(3)/L;
%r2 = (1+sqrt(1+4*k2*(xtemp(1)-xtemp(2))))/2;
%r1 = (1+sqrt(1+4*k1*(source-xtemp(1))))/2;
dx3(2,1) = xtemp(1)/(c2*r2) - xtemp(2)/(c2*r2) - xtemp(3)/c2;
dx3(1,1) = -xtemp(1)*(1/(c1*r1)+1/(c1*r2)) + xtemp(2)/(c1*r2) + source/(c1*r1);

dx3 = dx3 * dt;
xtemp = x + dx3;
dx4(3,1)= xtemp(2)/L - (1+k3*xtemp(3))*xtemp(3)/L;
%r2 = (1+sqrt(1+4*k2*(xtemp(1)-xtemp(2))))/2;
%r1 = (1+sqrt(1+4*k1*(source-xtemp(1))))/2;
dx4(2,1) = xtemp(1)/(c2*r2) - xtemp(2)/(c2*r2) - xtemp(3)/c2;
dx4(1,1) = -xtemp(1)*(1/(c1*r1)+1/(c1*r2)) + xtemp(2)/(c1*r2) + source/(c1*r1);

dx4 = dx4 * dt;
return;