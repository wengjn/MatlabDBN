%This is the system equation to simulate electrical LC system
function [xk] = sysLC(k, sysxkm1, uk)

xk = 0.5*sysxkm1 + 25*sysxkm1/(1+sysxkm1.^2) + 8*cos(1.2*(k-1)) + uk;


