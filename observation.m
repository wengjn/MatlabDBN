%observation equation to simulate the electrical LC system

function [yk] = observation(k, xk, vk)

%linearize the states variable in this system
%vk is the observation noise
yk = xk(1) + 2*xk(2) + 3 * xk(3) + vk;